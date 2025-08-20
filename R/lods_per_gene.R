# utilizing the qtl2rest container
# this will load the data

source('/app/qtl2rest/load.R')

#' Converts a named list map to a tibble with marker_id, chr, and pos
map_list_to_tbl <- function(map) {
    # Converts named list map to a tibble with marker_id, chr, and pos
    df_list <- lapply(names(map), function(chr) {
        data.frame(
            marker_id = names(map[[chr]]),
            chr = chr,
            pos = as.numeric(map[[chr]]),
            stringsAsFactors = FALSE
        )
    })
    tibble::as_tibble(do.call(rbind, df_list))
}


# the main worker - FIXED VERSION
find_qtls_per_gene <- function(scan1_obj, map, annotations, window_size = 2.5) {

    # Step 0: Convert map to tibble
    map_tbl <- map_list_to_tbl(map)
    markers <- map_tbl |>
        dplyr::filter(!is.na(pos)) |>
        janitor::clean_names() |>
        dplyr::select(marker_id, chr, pos)

    # fix if marker `pos` is in bp and convert to Mb
    if (all(markers$pos > 2e5, na.rm = TRUE)) {
        markers$pos <- markers$pos / 1e6
    }

    # Step 1: Convert scan1 to long format
    # marker_id, gene_id, lod
    lod_long <- scan1_obj |>
        # columns are gene ids
        as.data.frame() |>
        # add marker_id as column from rowname
        tibble::rownames_to_column('marker_id') |>
        # pivot to marker_id, gene_id, lod
        tidyr::pivot_longer(
            cols = -marker_id,
            names_to = 'gene_id',
            values_to = 'lod',
            values_transform = list(lod = as.numeric)
        )

    # Step 2: Join LODs with marker info
    # marker_id, gene_id, lod, chr, pos
    lod_joined <- dplyr::inner_join(lod_long, markers, by = 'marker_id')

    # Step 3: Clean annotations - FIXED: Don't filter by chromosomes in map
    # Only keep genes that are actually in the scan1 object
    scan1_genes <- unique(lod_long$gene_id)
    annotations_cleaned <- annotations |>
        dplyr::filter(gene_id %in% scan1_genes) |>
        dplyr::distinct(gene_id, symbol, chr_gene = chr, start)

    # fix if annotation `start` is in bp and convert to Mb
    if (all(annotations_cleaned$start > 2e5, na.rm = TRUE)) {
        annotations_cleaned$start <- annotations_cleaned$start / 1e6
    }

    # Step 4: Add nearest marker info - FIXED: Add error handling for qtl2::find_marker
    nearest_tbl <- annotations_cleaned |>
        dplyr::mutate(
            nearest_id = purrr::map2_chr(chr_gene, start, function(chr, pos) {
                tryCatch({
                    qtl2::find_marker(map, chr, pos)
                }, error = function(e) {
                    NA_character_  # Return NA if chromosome not in map or other error
                })
            })
        )

    # Step 5: Join nearest marker LOD - FIXED: Handle NA nearest_id
    nearest_info <- nearest_tbl |>
        dplyr::left_join(
            lod_joined,
            by = c('gene_id', 'nearest_id' = 'marker_id')
        ) |>
        dplyr::rename(
            nearest_id = nearest_id,
            nearest_chr = chr,
            nearest_pos = pos,
            nearest_lod = lod
        )

    # Step 6: Get global max LOD
    max_info <- lod_joined |>
        dplyr::group_by(gene_id) |>
        dplyr::slice_max(order_by = lod, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::rename(
            max_id = marker_id,
            max_chr = chr,
            max_pos = pos,
            max_lod = lod
        )

    # Step 7: Add max in-region and out-of-region LODs
    region_lods <- annotations_cleaned |>
        dplyr::inner_join(lod_joined, by = 'gene_id') |>
        dplyr::mutate(
            lower = pmax(0, start - window_size / 2.0),
            upper = start + window_size / 2.0,
            in_region = (chr == chr_gene) & (pos >= lower) & (pos <= upper)
        ) |>
        dplyr::group_by(gene_id) |>
        dplyr::summarize(
            # max within region
            local_id = if (any(in_region)) marker_id[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            local_chr = if (any(in_region)) chr[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            local_pos = if (any(in_region)) pos[which.max(replace(lod, !in_region, -Inf))] else NA_real_,
            local_lod = if (any(in_region)) max(lod[in_region]) else NA_real_,

            # max outside region
            distal_id = if (any(!in_region)) marker_id[which.max(replace(lod, in_region, -Inf))] else NA_character_,
            distal_chr = if (any(!in_region)) chr[which.max(replace(lod, in_region, -Inf))] else NA_character_,
            distal_pos = if (any(!in_region)) pos[which.max(replace(lod, in_region, -Inf))] else NA_real_,
            distal_lod = if (any(!in_region)) max(lod[!in_region]) else NA_real_
        ) |>
        dplyr::ungroup()

    # Step 8: Combine everything
    # Final column ordering
    result <- nearest_info |>
        dplyr::left_join(max_info, by = "gene_id") |>
        dplyr::left_join(region_lods, by = "gene_id") |>
        dplyr::select(
            gene_id, symbol, chr_gene, start,
            nearest_id, nearest_chr, nearest_pos, nearest_lod,
            max_id,     max_chr,     max_pos,     max_lod,
            local_id,   local_chr,   local_pos,   local_lod,
            distal_id,  distal_chr,  distal_pos,  distal_lod
        )
}



find_qtls_per_gene_ext <- function(scan1_obj, map, annotations, window_size = 2.5) {

    # Step 0: Convert map to tibble
    cat('Converting map to tibble\n')
    map_tbl <- map_list_to_tbl(map)
    markers <- map_tbl |>
        dplyr::filter(!is.na(pos)) |>
        janitor::clean_names() |>
        dplyr::select(marker_id, chr, pos)

    # fix if marker `pos` is in bp and convert to Mb
    if (all(markers$pos > 2e5, na.rm = TRUE)) {
        markers$pos <- markers$pos / 1e6
    }

    # Step 1: Convert scan1 to long format
    # marker_id, gene_id, lod
    cat('Converting scan1 object to long format\n')
    lod_long <- scan1_obj |>
        # columns are gene ids
        as.data.frame() |>
        # add marker_id as column from rowname
        tibble::rownames_to_column('marker_id') |>
        # pivot to marker_id, gene_id, lod
        tidyr::pivot_longer(
            cols = -marker_id,
            names_to = 'gene_id',
            values_to = 'lod',
            values_transform = list(lod = as.numeric)
        )

    # Step 2: Join LODs with marker info
    # marker_id, gene_id, lod, chr, pos
    cat('Joining LODs with marker info\n')
    lod_joined <- dplyr::inner_join(lod_long, markers, by = 'marker_id')

    # Step 3: Clean annotations - FIXED: Don't filter by chromosomes in map
    # Only keep genes that are actually in the scan1 object
    cat('Cleaning annotations\n')
    scan1_genes <- unique(lod_long$gene_id)
    annotations_cleaned <- annotations |>
        dplyr::filter(gene_id %in% scan1_genes) |>
        dplyr::distinct()

    # fix if annotation `start` is in bp and convert to Mb
    if (all(annotations_cleaned$start > 2e5, na.rm = TRUE)) {
        annotations_cleaned$start <- annotations_cleaned$start / 1e6
    }

    # Step 4: Add nearest marker info - FIXED: Add error handling for qtl2::find_marker
    cat('Finding nearest markers\n')
    nearest_tbl <- annotations_cleaned |>
        dplyr::mutate(
            nearest_id = purrr::map2_chr(chr_gene, start, function(chr, pos) {
                tryCatch({
                    qtl2::find_marker(map, chr, pos)
                }, error = function(e) {
                    NA_character_  # Return NA if chromosome not in map or other error
                })
            })
        )

    # Step 5: Join nearest marker LOD - FIXED: Handle NA nearest_id
    cat('Joining nearest marker LODs\n')
    nearest_info <- nearest_tbl |>
        dplyr::left_join(
            lod_joined,
            by = c('gene_id', 'nearest_id' = 'marker_id')
        ) |>
        dplyr::rename(
            nearest_id = nearest_id,
            nearest_chr = chr,
            nearest_pos = pos,
            nearest_lod = lod
        )

    # Step 6: Get global max LOD
    cat('Finding global max LODs\n')
    max_info <- lod_joined |>
        dplyr::group_by(gene_id) |>
        dplyr::slice_max(order_by = lod, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::rename(
            max_id = marker_id,
            max_chr = chr,
            max_pos = pos,
            max_lod = lod
        )

    # Step 7: Add max in-region and out-of-region LODs
    cat('Adding max in-region and out-of-region LODs\n')
    region_lods <- annotations_cleaned |>
        dplyr::inner_join(lod_joined, by = 'gene_id') |>
        dplyr::mutate(
            lower = pmax(0, start - window_size / 2.0),
            upper = start + window_size / 2.0,
            in_region = (chr == chr_gene) & (pos >= lower) & (pos <= upper)
        ) |>
        dplyr::group_by(gene_id) |>
        dplyr::summarize(
            # max within region
            local_id = if (any(in_region)) marker_id[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            local_chr = if (any(in_region)) chr[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            local_pos = if (any(in_region)) pos[which.max(replace(lod, !in_region, -Inf))] else NA_real_,
            local_lod = if (any(in_region)) max(lod[in_region]) else NA_real_,

            # max outside region
            distal_id = if (any(!in_region)) marker_id[which.max(replace(lod, in_region, -Inf))] else NA_character_,
            distal_chr = if (any(!in_region)) chr[which.max(replace(lod, in_region, -Inf))] else NA_character_,
            distal_pos = if (any(!in_region)) pos[which.max(replace(lod, in_region, -Inf))] else NA_real_,
            distal_lod = if (any(!in_region)) max(lod[!in_region]) else NA_real_
        ) |>
        dplyr::ungroup()

    cat('Finding max_ext LODs\n')
    region_lods_max_ext <- annotations_cleaned |>
        dplyr::inner_join(lod_joined, by = 'gene_id') |>
        dplyr::mutate(
            lower = pmax(0, ext_max_pos - window_size / 2.0),
            upper = ext_max_pos + window_size / 2.0,
            in_region = (chr == ext_max_chr) & (pos >= lower) & (pos <= upper)
        ) |>
        dplyr::group_by(gene_id) |>
        dplyr::summarize(
            # max within region
            max_lod_ext = if (any(in_region)) max(lod[in_region]) else NA_real_,
            max_pos_ext = if (any(in_region)) pos[which.max(replace(lod, !in_region, -Inf))] else NA_real_,
            max_chr_ext = if (any(in_region)) chr[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            max_id_ext = if (any(in_region)) marker_id[which.max(replace(lod, !in_region, -Inf))] else NA_character_
        ) |>
        dplyr::ungroup()

    cat('Finding distal_ext LODs\n')
    region_lods_distal_ext <- annotations_cleaned |>
        dplyr::inner_join(lod_joined, by = 'gene_id') |>
        dplyr::mutate(
            lower = pmax(0, ext_max_pos - window_size / 2.0),
            upper = ext_max_pos + window_size / 2.0,
            in_region = (chr == ext_max_chr) & (pos >= lower) & (pos <= upper)
        ) |>
        dplyr::group_by(gene_id) |>
        dplyr::summarize(
            # max within region
            distal_lod_ext = if (any(in_region)) max(lod[in_region]) else NA_real_,
            distal_pos_ext = if (any(in_region)) pos[which.max(replace(lod, !in_region, -Inf))] else NA_real_,
            distal_chr_ext = if (any(in_region)) chr[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            distal_id_ext = if (any(in_region)) marker_id[which.max(replace(lod, !in_region, -Inf))] else NA_character_
        ) |>
        dplyr::ungroup()

    # Combine max and distal ext LODs
    cat('Combining max and distal ext LODs\n')
    region_lods_ext <- region_lods_max_ext |>
        dplyr::left_join(region_lods_distal_ext, by = "gene_id")


    # Step 8: Combine everything
    # Final column ordering
    cat('Combining all results\n')
    result <- nearest_info |>
        dplyr::left_join(max_info, by = "gene_id") |>
        dplyr::left_join(region_lods, by = "gene_id") |>
        dplyr::left_join(region_lods_ext, by = "gene_id") |>
        dplyr::select(
            gene_id, symbol, chr_gene, start,
            nearest_id, nearest_chr, nearest_pos, nearest_lod,
            max_id,        max_chr,        max_pos,        max_lod,
            local_id,      local_chr,      local_pos,      local_lod,
            distal_id,     distal_chr,     distal_pos,     distal_lod,
            max_id_ext,    max_chr_ext,    max_pos_ext,    max_lod_ext,
            distal_id_ext, distal_chr_ext, distal_pos_ext, distal_lod_ext
        )
}


retrieve_lods <- function(dataset_id, start_id, end_id) {

    # use qtl2api helper functions
    dataset <- qtl2api::get_dataset_by_id(dataset_id)
    dataset <- qtl2api::synchronize_dataset(dataset)

    ds_addcovar <- qtl2api:::get_covar_matrix(dataset)
    ds_addcovar <- ds_addcovar$covar_matrix

    # Map a bunch of genes
    ds_expr = dataset$data[, start_id:end_id]

    cat('Performing scan1\n')
    scan1_obj <- qtl2::scan1(genoprobs = genoprobs, pheno = ds_expr, kinship = K, addcovar = ds_addcovar, cores = 20)

    cat('Done scan, finding lods\n')
    new_data <- find_qtls_per_gene(scan1_obj[,,drop=FALSE], map, dataset$annot_mrna)

    cat('Done find_qtls_per_gene\n')

    return(new_data)
}


retrieve_lods_ext <- function(dataset_id, ext_file_name, start_id, end_id) {

    # use qtl2api helper functions
    dataset <- qtl2api::get_dataset_by_id(dataset_id)
    dataset <- qtl2api::synchronize_dataset(dataset)

    ds_addcovar <- qtl2api:::get_covar_matrix(dataset)
    ds_addcovar <- ds_addcovar$covar_matrix

    # Map a bunch of genes
    ds_expr = dataset$data[, start_id:end_id]

    cat('Performing scan1\n')
    scan1_obj <- qtl2::scan1(genoprobs = genoprobs, pheno = ds_expr, kinship = K, addcovar = ds_addcovar, cores = 20)

    cat('Done scan1, reading external data:', ext_file_name, '\n')
    ext_data <- readr::read_csv(ext_file_name, show_col_types = FALSE)
    ext_data <- ext_data |> dplyr::select(
        gene_id,
        symbol,
        chr_gene,
        start,
        ext_max_id = max_id,
        ext_max_chr = max_chr,
        ext_max_pos = max_pos,
        ext_max_lod = max_lod,
        ext_max_id_outside = max_id_outside,
        ext_max_chr_outside = max_chr_outside,
        ext_max_pos_outside = max_pos_outside,
        ext_max_lod_outside = max_lod_outside
    )

    cat('Finding lods\n')
    new_data <- find_qtls_per_gene_ext(scan1_obj[,,drop=FALSE], map, ext_data)

    cat('Done find_qtls_per_gene\n')

    return(new_data)
}


write_lods <- function(dataset_id, start_id, end_id) {
    new_data <- retrieve_lods(dataset_id, start_id, end_id)

    file_name <- make.names(paste0(dataset_id, "_", start_id, "_", end_id, ".csv"))

    cat('Dumping to file: ', file_name, '\n')
    readr::write_csv(new_data, file = file_name, col_names = TRUE)
    cat('Done\n')
}


write_lods_ext <- function(dataset_id, ext_file_name, start_id, end_id) {
    new_data <- retrieve_lods_ext(dataset_id, ext_file_name, start_id, end_id)

    file_name <- make.names(paste0(dataset_id, "_", start_id, "_", end_id, ".ext.csv"))

    cat('Dumping to file: ', file_name, '\n')
    readr::write_csv(new_data, file = file_name, col_names = TRUE)
    cat('Done\n')
}


# parse the args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 3) {
    # call it
    write_lods(args[1], as.integer(args[2]) + 1, as.integer(args[3]) + 1)
} else if (length(args) == 4) {
    write_lods_ext(args[1], args[2], as.integer(args[3]) + 1, as.integer(args[4]) + 1)
} else {
    cat('Usage: Rscript lods_per_gene.R <dataset_id> <start_id> <end_id>\n')
    cat('or\n')
    cat('Usage: Rscript lods_per_gene.R <dataset_id> <ext_file.csv> <start_id> <end_id>\n')
    quit(status = 1)
}


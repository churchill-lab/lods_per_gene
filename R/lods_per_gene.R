# utilizing the qtl2rest container
# this will load the data

source('/app/qtl2rest/load.R')


# helper function
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

# the main worker
find_qtls_per_gene <- function(scan1_obj, map, annotations, window_size = 5.0) {

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
    lod_long <- scan1_obj |>
        as.data.frame() |>
        tibble::rownames_to_column('marker_id') |>
        tidyr::pivot_longer(
            cols = -marker_id,
            names_to = 'gene_id',
            values_to = 'lod',
            values_transform = list(lod = as.numeric)
        )

    # Step 2: Join LODs with marker info
    lod_joined <- dplyr::inner_join(lod_long, markers, by = 'marker_id')

    # Step 3: Clean annotations
    annotations_cleaned <- annotations |>
        dplyr::filter(chr %in% names(map)) |>
        dplyr::distinct(gene_id, symbol, chr_gene = chr, start)

    # fix if annotation `start` is in bp and convert to Mb
    if (all(annotations_cleaned$start > 2e5, na.rm = TRUE)) {
        annotations_cleaned$start <- annotations_cleaned$start / 1e6
    }

    # Step 4: Add nearest marker info
    nearest_tbl <- annotations_cleaned |>
        dplyr::mutate(
            nearest_id = purrr::map2_chr(chr_gene, start, ~ qtl2::find_marker(map, .x, .y))
        )

    # Step 5: Join nearest marker LOD
    nearest_info <- dplyr::inner_join(
        nearest_tbl,
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
            max_lod_within = if (any(in_region)) max(lod[in_region]) else NA_real_,
            max_pos_within = if (any(in_region)) pos[which.max(replace(lod, !in_region, -Inf))] else NA_real_,
            max_chr_within = if (any(in_region)) chr[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            max_id_within = if (any(in_region)) marker_id[which.max(replace(lod, !in_region, -Inf))] else NA_character_,

            # max outside region
            max_lod_outside = if (any(!in_region)) max(lod[!in_region]) else NA_real_,
            max_pos_outside = if (any(!in_region)) pos[which.max(replace(lod, in_region, -Inf))] else NA_real_,
            max_chr_outside = if (any(!in_region)) chr[which.max(replace(lod, in_region, -Inf))] else NA_character_,
            max_id_outside = if (any(!in_region)) marker_id[which.max(replace(lod, in_region, -Inf))] else NA_character_
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
            max_id_within, max_chr_within, max_pos_within, max_lod_within,
            max_id_outside, max_chr_outside, max_pos_outside, max_lod_outside
        )
}


# the main worker - FIXED VERSION
find_qtls_per_gene_y <- function(scan1_obj, map, annotations, window_size = 5.0) {

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
    lod_long <- scan1_obj |>
        as.data.frame() |>
        tibble::rownames_to_column('marker_id') |>
        tidyr::pivot_longer(
            cols = -marker_id,
            names_to = 'gene_id',
            values_to = 'lod',
            values_transform = list(lod = as.numeric)
        )

    # Step 2: Join LODs with marker info
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
            max_lod_within = if (any(in_region)) max(lod[in_region]) else NA_real_,
            max_pos_within = if (any(in_region)) pos[which.max(replace(lod, !in_region, -Inf))] else NA_real_,
            max_chr_within = if (any(in_region)) chr[which.max(replace(lod, !in_region, -Inf))] else NA_character_,
            max_id_within = if (any(in_region)) marker_id[which.max(replace(lod, !in_region, -Inf))] else NA_character_,

            # max outside region
            max_lod_outside = if (any(!in_region)) max(lod[!in_region]) else NA_real_,
            max_pos_outside = if (any(!in_region)) pos[which.max(replace(lod, in_region, -Inf))] else NA_real_,
            max_chr_outside = if (any(!in_region)) chr[which.max(replace(lod, in_region, -Inf))] else NA_character_,
            max_id_outside = if (any(!in_region)) marker_id[which.max(replace(lod, in_region, -Inf))] else NA_character_
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
            max_id_within, max_chr_within, max_pos_within, max_lod_within,
            max_id_outside, max_chr_outside, max_pos_outside, max_lod_outside
        )
}

write_lods <- function(dataset_id, start_id, end_id) {

    # use qtl2api helper functions
    dataset <- qtl2api::get_dataset_by_id(dataset_id)
    dataset <- qtl2api::synchronize_dataset(dataset)

    ds_addcovar <- qtl2api:::get_covar_matrix(dataset)
    ds_addcovar <- ds_addcovar$covar_matrix

    # Map a bunch of genes
    ds_expr = dataset$data[, start_id:end_id]

    print('Performing scan1')
    scan1_obj <- qtl2::scan1(genoprobs = genoprobs, pheno = ds_expr, kinship = K, addcovar = ds_addcovar, cores = 20)

    print('Done scan, finding lods')
    new_data <- find_qtls_per_gene_y(scan1_obj[,,drop=FALSE], map, dataset$annot_mrna)

    print('Done finding local and max')
    file_name <- make.names(paste0(dataset_id, "_", start_id, "_", end_id, ".csv"))

    print(paste0('Dumping to file: ', file_name))
    readr::write_csv(new_data, file = file_name, col_names = TRUE)

    print('Done')
}

# parse the args
args <- commandArgs(trailingOnly = TRUE)

# call it
write_lods(args[1], as.integer(args[2]) + 1, as.integer(args[3]) + 1)



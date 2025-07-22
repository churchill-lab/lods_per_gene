#!/bin/bash

# Usage: ./submit.sh <container_file> <r_script_file> <data_dir> <db_file> <dataset> <batch_size> <total_rows_in_data>

# Example:
# /pat/to/bin/submit.sh \ # THIS SCRIPT (IN REPO)
# /path/to/singularity/lods_per_gene.sif \ # SINGULARITY FILE, (CAN BE BUILT FROM lods_per_gene.def)
# /path/to/R/lods_per_gene.R \ MAIN SCRIPT (IN REPO)
# /datapath/to/rdata \ # DATA DIRECTORY, SHOULD HAVEL ALL RData AND rds FILES 
# /datapath/snps/fv.2021.snps.db3 \ # QTL2API USES FOR SNPS
# dataset.test.v3 \ # NAME OF DATASET
# 100 \ # BATCH SIZE (NUMBER OF GENES PER BATCH)
# 22117  # TOTAL # OF GENES IN DATASET

# === Input Arguments ===
container=$1
r_script=$2
data_dir=$3
db_file=$4
dataset=$5
batch_size=$6
total_rows=$7

# === Validate Input ===
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <container_file> <r_script_file> <data_dir> <db_file> <dataset> <batch_size> <total_rows>"
    exit 1
fi

# === Compute number of jobs ===
num_jobs=$(( (total_rows + batch_size - 1) / batch_size ))  # rounds up
max_array_index=$((num_jobs - 1))

echo "Submitting SLURM array job with $num_jobs tasks to process $total_rows rows in batches of $batch_size"

# === Submit Job ===
sbatch --array=0-${max_array_index} \
       --export=ALL,CONTAINER=$container,R_SCRIPT=$r_script,DATA_DIR=$data_dir,DB_FILE=$db_file,DATASET=$dataset,BATCH_SIZE=$batch_size,TOTAL_ROWS=$total_rows \
       /flashscratch/mvincent/gbrs_qtl_lods/lods_per_gene/bin/lods_per_gene.sh



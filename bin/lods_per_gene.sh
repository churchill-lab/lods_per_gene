#!/bin/bash

#SBATCH --job-name=lods_per_gene         # Job name
#SBATCH --mail-type=NONE                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=matt.vincent@jax.org # Where to send mail
#SBATCH --nodes=1                        # How many nodes
#SBATCH --cpus-per-task=20               # CPUs per task
#SBATCH --mem=64gb                       # Memory
#SBATCH --time=01:00:00                  # Time limit hrs:min:sec
#SBATCH --output=lods_per_gene_%A-%a.out # Standard output and error log

pwd; hostname; date

module load singularity

# Read from environment variables exported by submit_peaks.sh
container=$CONTAINER
r_script=$R_SCRIPT
data_dir=$DATA_DIR
db_file=$DB_FILE
dataset=$DATASET
batch_size=$BATCH_SIZE
total_rows=$TOTAL_ROWS

# Calculate start and end indices
start_idx=$((SLURM_ARRAY_TASK_ID * batch_size))
end_idx=$((start_idx + batch_size - 1))
if [ $end_idx -ge $total_rows ]; then
    end_idx=$((total_rows - 1))
fi

echo "[$SLURM_ARRAY_TASK_ID] Processing rows $start_idx to $end_idx"

echo "container=$container"
echo "r_script=$r_script"
echo "data_dir=$data_dir"
echo "db_file=$db_file"
echo "dataset=$dataset"
echo "start_idx=$start_idx"
echo "end_idx=$end_idx"

singularity exec --bind $data_dir:/app/qtl2rest/data/rdata --bind $db_file:/app/qtl2rest/data/sqlite/ccfounders.sqlite $container Rscript $r_script $dataset $start_idx $end_idx

date


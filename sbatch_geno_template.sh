#!/bin/bash

#SBATCH --array=1-22
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=6-00:00:00
#SBATCH --partition=batch
#SBATCH --job-name=amr_1000_geno
#SBATCH --output=R-%x.%j.out

n=$SLURM_ARRAY_TASK_ID

cd /syntetic_data/

CONFIG=data/amr/config_1000 # prefix for config file

# generate a config for each chromosome
cp ${CONFIG}.yaml ${CONFIG}$n.yaml
sed -i 's/${chr}'"/$n/g" ${CONFIG}$n.yaml

mkdir -p /syntetic_data/logs

TIMEFILE="/syntetic_data/logs/time_geno_amr_1000_chr${n}.txt"

echo "=== GENO amr 1000 ===" > "$TIMEFILE"
echo "Start: $(date)" >> "$TIMEFILE"

/usr/bin/time -v -a -o "$TIMEFILE" singularity exec \
  --bind /syntetic_data/data:/data \
  --env JULIA_DEPOT_PATH=/data/.julia \
 /syntetic_data/synthetic-genetic-data.sif \
  generate_geno 8 "${CONFIG}${n}.yaml"

echo "End: $(date)" >> "$TIMEFILE"

#!/bin/bash
#SBATCH --job-name=amr_1000_pheno
#SBATCH --output=R-%x.%j.out
#SBATCH --time=6-00:00:00
#SBATCH --mem=256G
#SBATCH --partition=batch

cd /syntetic_data/

CONFIG=data/amr/config_1000 # prefix for config file

# replace the chromosome wildcard with all
cp ${CONFIG}.yaml ${CONFIG}_pheno.yaml
sed -i 's/${chr}'"/all/g" ${CONFIG}_pheno.yaml

  echo "=== PHENO amr 1000 ==="
  echo "Start: $(date)"
  echo "--- generate_pheno ---"
  
  singularity exec \
  --bind /syntetic_data/data:/data \
  --env JULIA_DEPOT_PATH=/data/.julia \
  /syntetic_data/synthetic-genetic-data.sif \
  generate_pheno ${CONFIG}_pheno.yaml

  echo "--- convert ---"
  singularity exec \
  --bind /syntetic_data/data:/data \
  --env JULIA_DEPOT_PATH=/data/.julia \
  /syntetic_data/synthetic-genetic-data.sif \
  convert ${CONFIG}_pheno.yaml

echo "End: $(date)"

#!/bin/bash
#SBATCH --job-name=amr_1000_eval
#SBATCH --output=R-%x.%j.out
#SBATCH --time=6-00:00:00
#SBATCH --mem=128G
#SBATCH --partition=batch

cd /syntetic_data/

CONFIG=data/amr/config_1000 # prefix for config file

SIZE=$(basename "${CONFIG}" | sed 's/^config_//')

mkdir -p /syntetic_data/logs

run_eval () {
  local chr_value=$1       
  local suffix=$2           
  local aats=$3 kinship=$4 ld_corr=$5 ld_decay=$6 maf=$7 pca=$8 gwas=$9 mia=${10}

  local yaml="${CONFIG}${suffix}_evaluation.yaml"

  cp "${CONFIG}.yaml" "$yaml"
  sed -i 's/${chr}'"/${chr_value}/g" "$yaml"

  sed -i \
    -e "s/aats: \(true\|false\)/aats: ${aats}/" \
    -e "s/kinship: \(true\|false\)/kinship: ${kinship}/" \
    -e "s/ld_corr: \(true\|false\)/ld_corr: ${ld_corr}/" \
    -e "s/ld_decay: \(true\|false\)/ld_decay: ${ld_decay}/" \
    -e "s/maf: \(true\|false\)/maf: ${maf}/" \
    -e "s/pca: \(true\|false\)/pca: ${pca}/" \
    -e "s/gwas: \(true\|false\)/gwas: ${gwas}/" \
    -e "s/mia: \(true\|false\)/mia: ${mia}/" \
    "$yaml"

  
    echo "=== Evaluation amr 1000${suffix} (chr=${chr_value}) ==="
    echo "Start: $(date)"
    echo "--- validate ---"
    singularity exec \
      --bind /syntetic_data/data:/data \
      --env JULIA_DEPOT_PATH=/data/.julia \
      /syntetic_data/synthetic-genetic-data.sif \
      validate "$yaml"
    echo "End: $(date)"
  
}

case "$SIZE" in
  1000)
    # all chr: only pca, maf, ld_corr, ld_decay, gwas
    run_eval "all" ""       false false true  true  true  true  true  false
    # chr 21: also aats, kinship, mia
    run_eval "21"  "_chr21" true  true  true  true  true  true  true  true
    ;;

  10000)
    # all true, all chr
    run_eval "all" "" true true true true true true true true
    ;;

  100000|1000000|2000000)
    # only pca, maf, ld_corr, ld_decay, gwas
    run_eval "all" "" false false true true true true true false
    ;;

  *)
    exit 1
    ;;
esac
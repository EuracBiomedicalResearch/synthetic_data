# How to use this tool?

For a more comprehensive guide, see the **README** file.

---

## Step 1: Pull the container

Pull the container (Docker or Singularity) using this command:

```bash
singularity pull docker://benonsen/synthetic-genetic-data
```

If some singularity out of memory issues appear, set the cache and tmp dir before pulling the container.
```bash
export SINGULARITY_TMPDIR=<path-to-tmp>/tmp
export SINGULARITY_CACHEDIR=<path-to-tmp>/tmp/cache
```

---

## Step 2: Initialization

Run `init` and install all the needed packages using this command. Also, bind the data folder to make it persistent; move the `config.yaml` file into `data/config.yaml` beforehand.

If you encounter any Julia-related errors, add this after the bind parameter to the command:

```
--env JULIA_DEPOT_PATH=/data/.julia
```

Command:

```bash
singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif init
```

---

## Step 3: Fetch reference data

Fetch all the reference data. This will download the reference dataset:

```bash
singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif fetch
```

---

## Step 4: Generate data

To generate data, run:

```bash
singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif generate_geno <number-of-threads> data/config.yaml
```

In practice, we do not use this single command; instead, we use a script that automatically submits all the jobs with the correct dependencies. I will move the files to the correct directory.

---

## File structure

**What files do we have?**

* `/data/inputs/processed/1KG+HGDP/`:
  Reference dataset from 1KG. This path contains 5 folders: `EUR`, `CSA`, `EAS`, `AMR`, `AFR`.
  Each folder contains the `CausalList` file for a trait.

  * Trait 1 corresponds to `CausalList1`
  * Trait 2 corresponds to `CausalList2`
  * etc.
    First, specify the SNP and then the effect size.
    I provided a script (see create_causalList.py) to convert GWAS Catalog associations into a `CausalList` file automatically. As input provide the downloaded file from GWASCatalog and the output file.

* `/data/`:
  Also contains 5 folders (`eur`, `csa`, `eas`, `amr`, `afr`) with the corresponding `config.yaml` files.

  * For each sample size, a new config file was created.
  * Make sure the path for the `CausalList` points to the correct folder and the correct population is specified. (I already checked it; for this case, it is fine.)
  * The `chromosome` parameter can be `"all"` or an integer between `1` and `22`. Since we are operating in a cluster, we set it to `${chr}` and replace it with a number between `1` and `22` to generate in parallel.

After generating the genotype, we generate the phenotype, validate, and convert (all sequentially). The script automatically replaces the `${chr}` placeholder with `"all"`.

---

## Step 5: Submit jobs

After running the `init` and `fetch` commands, just run the `submit_jobs.sh` script to submit all the jobs.

**Important:**
This script uses two underlying template scripts and replaces the sample size and population identifier. Therefore, please check the paths in those scripts before running `submit_jobs.sh`.

The template scripts are:

* `sbatch_geno_template.sh`
* `sbatch_pheno_template.sh` (this one generates phenotypes, validates, and converts).

---

## Output

* Results can be found here:

  ```
  /data/outputs/test/{population}/{sample size}
  ```

* VCF files are in the same path but in a dedicated `VCF` folder.

* Evaluation results can be found here:

  ```
  /data/outputs/test/{population}/{sample size}/evaluation
  ```

This folder contains multiple `sumstats` files: one for each trait and each chromosome.

* For summary statistics of a trait across all chromosomes, see:

  ```
  test_chr-1.traitX.all.sumstat
  ```

* This path also includes all **Manhattan plots**, **QQ plots**, and intermediate results (`pca`, `glm.linear` files, etc.).


If you want to generate data on you own and do not use my scripts, you can do the following: 
Create a new directory and start clean and use the pull, init (to avoid this step, just copy my container) and fetch command to download the reference data. Make sure to put the `config.yaml` file in it. Alternatively to downloading you can also copy all the files from `/data/inputs/processed/1KG+HGDP` (you don't need the population specific folders; in those I stored only the CausalList Files (Causal SNPs) for the corresponding population). Then you can start generating genes with 
  ```bash
  singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif generate_geno <number-of-threads> data/config.yaml
  ```
  and phenotypes: 
  ```bash
  singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif generate_pheno data/config.yaml
  ```
  validate (set gwas to `true` in config.yaml to calculate summary statistics): 
  ```bash
  singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif validate data/config.yaml
  ```
  and finally convert to vcf files:
  ```bash
  singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif validate data/config.yaml
  ```
  The config.yaml file is rather self explenatory; the most important parameters are: 
  - `chromosome`: which chromosome to you want to generate, or all
  - `causal_list`: Path to the CausalList-files. They don't have any file specific ending, just a text file. The last character in the filename should be the trait it reffers to, so CausalList1 --> Trait1. If you want to use it, you need also to set `UseCausalList` to true
  - nsamples: To change the sample size, note that you can use a default population structure or a custom one. 
  - `nTrait`: number of traits you want to simulate. If you change this value, you also need to update the matrices (PropotionGeno and so on) as they are nTrait x nTrait or nPopulation x nTrait matrices. See the official README for more info. 



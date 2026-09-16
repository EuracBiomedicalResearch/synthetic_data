"""Executes the pipeline for evaluating synthetic data quality
"""

include("../utils/reference_data.jl")
include("metrics/eval_aats.jl")
include("metrics/eval_kinship_detail.jl")
include("metrics/eval_ld_corr.jl")
include("metrics/eval_ld_decay.jl")
include("metrics/eval_maf.jl")
include("metrics/eval_pca.jl")
include("metrics/eval_gwas.jl")
 

function run_kinship_evaluation(ibsfile_real, ibsfile_synt, ibsfile_cross)
    run_kinship(ibsfile_real, ibsfile_synt, ibsfile_cross)
end


function run_aats_evaluation(ibsfile_cross)
    run_aats(ibsfile_cross)
end


function run_ld_corr_evaluation(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
    run_ld_corr(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
end


function run_ld_decay_evaluation(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map, chromosome)
    run_ld_decay(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map, chromosome)
end


function run_maf_evaluation(real_maf_file, synt_maf_file)
    run_maf(real_maf_file, synt_maf_file)
end


function run_pca_evaluation(real_data_pca_prefix, synt_data_pca_prefix, pcaproj_file, real_data_prefix, synt_data_prefix, eval_dir)
    real_data_pop = string(real_data_prefix, ".sample")
    synt_data_pop = string(synt_data_prefix, ".sample")
    run_pca(real_data_pca_prefix, synt_data_pca_prefix, pcaproj_file, real_data_pop, synt_data_pop, eval_dir)
end


function run_gwas_evaluation(ntraits, plink2, covar, synt_data_prefix, eval_dir, chromosome)
    run_gwas(ntraits, plink2, covar, synt_data_prefix, eval_dir, chromosome)
end


function run_mia_evaluation(orig_prefix, gen_prefix, out_dir, balancing, pca_components, train_fraction, shuffle, random_state)
    @info "Running MIA (Python)"
    py_script = joinpath(@__DIR__, "metrics", "eval_mia.py")
    println(`python3 $py_script --orig_prefix $orig_prefix --gen_prefix $gen_prefix --out_dir $out_dir --balancing $balancing --pca_components $pca_components --train_fraction $train_fraction --shuffle $shuffle --random_state $random_state`)
    
    run(`python3 $py_script --orig_prefix $orig_prefix --gen_prefix $gen_prefix --out_dir $out_dir --balancing $balancing --pca_components $pca_components --train_fraction $train_fraction --shuffle $shuffle --random_state $random_state`)
end


"""Computations for MAF using PLINK
"""
function run_maf_tools(plink, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for MAF"
    reffile_out = @sprintf("%s.ref.maf", outdir)
    synfile_out = @sprintf("%s.syn.maf", outdir)
    run(`$plink --bfile $reffile_prefix --freq --out $reffile_out`)
    run(`$plink --bfile $synfile_prefix --freq --out $synfile_out`)
    real_maffile = @sprintf("%s.frq", reffile_out)
    syn_maffile = @sprintf("%s.frq", synfile_out)
    return real_maffile, syn_maffile
end


"""Computations for PCA using PLINK
"""

function run_pca_tools(plink2, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for PCA"
    reffile_out = @sprintf("%s.ref.pca", outdir)
    synfile_out = @sprintf("%s.syn.pca", outdir)
    ref_snps_file = @sprintf("%s.ref", outdir)
    syn_snps_file = @sprintf("%s.syn", outdir)

    run(`$plink2 --bfile $reffile_prefix --maf 0.05 --write-snplist --out $ref_snps_file`)
    run(`$plink2 --bfile $synfile_prefix --maf 0.05 --write-snplist --out $syn_snps_file`)

    ref_snps = readlines(ref_snps_file * ".snplist")
    syn_snps = readlines(syn_snps_file * ".snplist")
    common   = intersect(ref_snps, syn_snps)

    common_file_name = outdir * "common.snplist" 
    write(common_file_name, join(common, "\n"))

    run(`$plink2 --bfile $reffile_prefix --extract $common_file_name --pca approx allele-wts --freq counts --seed 42 --out $reffile_out`)
    run(`$plink2 --bfile $synfile_prefix --extract $common_file_name --pca approx allele-wts --freq counts --seed 42 --out $synfile_out`)

    # --- PCA projection via plink2 --score ---
    projfile_prefix = @sprintf("%s_proj", outdir)
    projfile_out = @sprintf("%spc.txt", projfile_prefix)

    ref_proj_prefix = @sprintf("%s_ref_scored", projfile_prefix)
    syn_proj_prefix = @sprintf("%s_syn_scored", projfile_prefix)

    # Dynamically determine which columns hold ID, A1, and the PCs in the
    # .eigenvec.allele header, instead of hardcoding column numbers (these
    # can shift depending on plink2 version/options and on how many PCs
    # were computed).
    allele_header = split(readline(reffile_out * ".eigenvec.allele"), "\t")
    
    id_col = findfirst(h -> h == "ID", allele_header)
    a1_col = findfirst(h -> h == "A1", allele_header)
    pc_start = findfirst(h -> startswith(h, "PC"), allele_header)
    pc_end = length(allele_header)
    if id_col === nothing || a1_col === nothing || pc_start === nothing
        error("Could not find expected ID/A1/PC columns in $(reffile_out).eigenvec.allele header: $allele_header")
    end
    score_col_range = "$pc_start-$pc_end"
    @info "Detected ID col $id_col, A1 col $a1_col, PC columns $score_col_range in eigenvec.allele header ($(pc_end - pc_start + 1) PCs)"

    run(`$plink2 --bfile $reffile_prefix --extract $common_file_name
        --score $reffile_out.eigenvec.allele $id_col $a1_col variance-standardize cols=-scoreavgs,+scoresums list-variants header-read
        --score-col-nums $score_col_range
        --read-freq $reffile_out.acount
        --out $ref_proj_prefix`)

    run(`$plink2 --bfile $synfile_prefix --extract $common_file_name
        --score $reffile_out.eigenvec.allele $id_col $a1_col variance-standardize cols=-scoreavgs,+scoresums list-variants header-read
        --score-col-nums $score_col_range
        --read-freq $reffile_out.acount
        --out $syn_proj_prefix`)

    ref_scores = CSV.File(ref_proj_prefix * ".sscore", normalizenames=true) |> DataFrame
    syn_scores = CSV.File(syn_proj_prefix * ".sscore", normalizenames=true) |> DataFrame

     # Strip leading "#" or "_" from header 
    rename!(ref_scores, [n => Symbol(replace(string(n), r"^[#_]" => "")) for n in names(ref_scores)])
    rename!(syn_scores, [n => Symbol(replace(string(n), r"^[#_]" => "")) for n in names(syn_scores)])

    ref_scores.AFF = fill(1, nrow(ref_scores))
    syn_scores.AFF = fill(2, nrow(syn_scores))

    # Normalize projected PCs

    n_score_vars_ref = countlines(ref_proj_prefix * ".sscore.vars")
    n_score_vars_syn = countlines(syn_proj_prefix * ".sscore.vars")
    eigenval = vec(readdlm(reffile_out * ".eigenval"))
    
    score_cols = filter(c -> occursin(r"^PC\d+_SUM$", string(c)), names(ref_scores))
    
    for (i, c) in enumerate(score_cols)
        ref_scores[!, c] ./= (n_score_vars_ref * sqrt(eigenval[i]))
        syn_scores[!, c]  ./= (n_score_vars_syn * sqrt(eigenval[i]))
    
        new_name = Symbol(replace(string(c), r"_SUM$" => ""))
        rename!(ref_scores, c => new_name)
        rename!(syn_scores, c => new_name)
    end

    combined = vcat(ref_scores, syn_scores, cols=:union)
    CSV.write(projfile_out, combined)

    return reffile_out, synfile_out, projfile_out
end


"""Computations for relatedness using KING
"""
function run_relatedness_tools(king, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for relatedness"
    reffile_bed = @sprintf("%s.bed", reffile_prefix)
    synfile_bed = @sprintf("%s.bed", synfile_prefix)
    reffile_out = @sprintf("%s.ref.king", outdir)
    synfile_out = @sprintf("%s.syn.king", outdir)
    crossfile_out = @sprintf("%s.cross.king", outdir)
    run(`$king -b $reffile_bed --ibs --prefix $reffile_out`)
    run(`$king -b $synfile_bed --ibs --prefix $synfile_out`)
    run(`$king -b $reffile_bed,$synfile_bed --ibs --prefix $crossfile_out`)
    real_ibsfile = @sprintf("%s.ibs0", reffile_out)
    syn_ibsfile = @sprintf("%s.ibs0", synfile_out)
    cross_ibsfile = @sprintf("%s.ibs0", crossfile_out)
    return real_ibsfile, syn_ibsfile, cross_ibsfile
end


"""Runs computations using external tools such as PLINK and KING, based on the selection of evaluation metrics
"""
function run_external_tools(metrics, reffile_prefix, synfile_prefix, filepaths)
    external_files = Dict()
    if metrics["maf"]
        external_files["real_maffile"], external_files["syn_maffile"] = run_maf_tools(filepaths.plink, reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["pca"] || metrics["gwas"]
        external_files["real_pcafile"], external_files["syn_pcafile"], external_files["pcaproj_file"] = run_pca_tools(
            filepaths.plink2, 
            reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["kinship"] || metrics["aats"]
        external_files["real_ibsfile"], external_files["syn_ibsfile"], external_files["cross_ibsfile"] = run_relatedness_tools(filepaths.king, reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    return external_files
end


"""Executes evaluation for the metrics specified in the configuration file
"""
function run_pipeline(options, chromosome, superpopulation, metrics, filepaths, genomic_metadata, reffile_prefix, nsamples_ref, synfile_prefix, external_files)
    if metrics["aats"]
        run_aats_evaluation(external_files["cross_ibsfile"])
    end
    if metrics["kinship"]
        run_kinship_evaluation(external_files["real_ibsfile"], external_files["syn_ibsfile"], external_files["cross_ibsfile"])
    end
    if metrics["ld_corr"]
        run_ld_corr_evaluation(reffile_prefix, synfile_prefix, filepaths.evaluation_output, filepaths.plink)
    end
    if metrics["ld_decay"]
        bp_to_cm_map = create_bp_cm_ref(filepaths.genetic_distfile)
        run_ld_decay_evaluation(synfile_prefix, reffile_prefix, filepaths.plink, filepaths.mapthin, filepaths.evaluation_output, bp_to_cm_map, chromosome)
    end
    if metrics["maf"]
        run_maf_evaluation(external_files["real_maffile"],  external_files["syn_maffile"])
    end
    if metrics["pca"]
        run_pca_evaluation(
            external_files["real_pcafile"], external_files["syn_pcafile"], external_files["pcaproj_file"],
            reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["mia"]
        mia_cfg = haskey(options["evaluation"], "mia") ? options["evaluation"]["mia"] : Dict{String,Any}()
        balancing = haskey(mia_cfg, "balancing") ? mia_cfg["balancing"] : "undersample"
        pca_components = haskey(mia_cfg, "pca_components") ? mia_cfg["pca_components"] : 20
        train_fraction = haskey(mia_cfg, "train_fraction") ? mia_cfg["train_fraction"] : 0.7
        shuffle = haskey(mia_cfg, "shuffle") ? mia_cfg["shuffle"] : true
        out_dir = filepaths.evaluation_output
        random_state = options["global_parameters"]["random_seed"]
        run_mia_evaluation(reffile_prefix, synfile_prefix, out_dir, balancing, pca_components, train_fraction, shuffle, random_state)
    end
    if metrics["gwas"]
        run_gwas_evaluation(options["phenotype_data"]["nTrait"], filepaths.plink2, @sprintf("%s.eigenvec", external_files["syn_pcafile"]), synfile_prefix, filepaths.evaluation_output, chromosome)
    end
    

end


"""Entry point to running the evaluation pipeline for genotype data

Note that the evaluation pipeline assumes that the synthetic data you want 
to evaluate has already been geneerated, using the setup specified in the
configuration file. It is therefore recommended to run the pipeline with the 
--genotype and --evaluation flags together, so that the program generates 
the data and then immediately evaluates it using the correct settings.
"""
function run_evaluation(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)

    metrics = options["evaluation"]["metrics"]
    gwas = metrics["gwas"] 
    if chromosome == "all"      
        covar_paths = String[]
        for chromosome_i in 1:22
            # Prepare per-chromosome filepaths and metadata
            filepaths = parse_filepaths(options, chromosome_i, superpopulation)
            genomic_metadata = parse_genomic_metadata(options, superpopulation, filepaths)

            reffile_prefix, nsamples_ref = create_reference_dataset(filepaths.vcf_input_processed, filepaths.popfile_processed, genomic_metadata.population_weights, filepaths.plink, filepaths.reference_dir, chromosome_i)
            synfile_prefix =  filepaths.synthetic_data_prefix

            # Ensure PCA is computed per chromosome if GWAS is requested
            local_metrics = copy(metrics)
            local_metrics["gwas"] = false
            if gwas
                local_metrics["pca"] = true
            end

            external_files = run_external_tools(local_metrics, reffile_prefix, synfile_prefix, filepaths)
            run_pipeline(options, chromosome_i, superpopulation, local_metrics, filepaths, genomic_metadata, reffile_prefix, nsamples_ref, synfile_prefix, external_files)

            # Collect per-chromosome covariate path (PCA eigenvectors)
            if gwas
                push!(covar_paths, @sprintf("%s.eigenvec", external_files["syn_pcafile"]))
            end
        end
        # If GWAS was originally requested, run it across all chromosomes now
        if gwas
            filepaths_all = parse_filepaths(options, 1, superpopulation)
            synfile_prefix_all = filepaths_all.synthetic_data_prefix
            run_gwas_evaluation(options["phenotype_data"]["nTrait"], filepaths_all.plink2, covar_paths, synfile_prefix_all, filepaths_all.evaluation_output, "all")
        end
    else
        filepaths = parse_filepaths(options, chromosome, superpopulation)
        genomic_metadata = parse_genomic_metadata(options, superpopulation, filepaths)

        reffile_prefix, nsamples_ref = create_reference_dataset(filepaths.vcf_input_processed, filepaths.popfile_processed, genomic_metadata.population_weights, filepaths.plink, filepaths.reference_dir, chromosome)
        synfile_prefix =  filepaths.synthetic_data_prefix

        external_files = run_external_tools(metrics, reffile_prefix, synfile_prefix, filepaths)
        run_pipeline(options, chromosome, superpopulation, metrics, filepaths, genomic_metadata, reffile_prefix, nsamples_ref, synfile_prefix, external_files)
    end
end
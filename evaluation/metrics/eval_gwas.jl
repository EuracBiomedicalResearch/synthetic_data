"""
    Run a fast GWAS and generate manhattan and qqplot
"""

using CSV
using DataFrames
using MendelPlots
using Printf


function run_gwas_tools(plink2, syngeno_prefix, synpheno_prefix, trait_idx, covar, outdir)
	@info  "Create plink phenotype file"
	# for some reason, the file name is wrong, it carries an additional "-<number>", therefore we remove it and the filename is valid again
	# this only applies to the pheno filename, the geno filename needs this additional number 	
	synpheno_prefix_removed = replace(synpheno_prefix, r"-\d+" => "")
	fam = CSV.File(syngeno_prefix * ".fam", normalizenames=true, header = 0) |> DataFrame
	pheno = CSV.File(synpheno_prefix_removed * ".pheno" * string(trait_idx), normalizenames=true) |> DataFrame
	# composed of iid and fid 
	#fam[!,"Sample"] = string.(fam[:,1], "_",fam[:,2])
	# just take the iid to join fam and pheno
	fam[!,"Sample"] = string.(fam[:,1])
	
	PhenoFam = outerjoin(fam, pheno, on = :Sample)
	PhenoFam = PhenoFam[:, ["Column1","Column2","Column3", "Column4", "Column5", "Phenotype_liability_"]]

	syngeno_prefix_phe_trait_idx = @sprintf("%s.phe%i", synpheno_prefix, trait_idx)
	CSV.write(syngeno_prefix_phe_trait_idx, DataFrame(PhenoFam), delim = "\t", header = false)
	
	@info  "GWAS using plink 2"
	syngeno_prefix_bed = @sprintf("%s.bed", syngeno_prefix)
	syngeno_prefix_bim = @sprintf("%s.bim", syngeno_prefix)
	run(`$plink2 --bed $syngeno_prefix_bed --bim $syngeno_prefix_bim --fam $syngeno_prefix_phe_trait_idx --glm hide-covar --covar $covar --covar-variance-standardize --ci 0.95 --out $outdir`)	
	
	@info  "Create summary statistics"
	GWASout = CSV.File( outdir * ".PHENO1.glm.linear", normalizenames=true) |> DataFrame
	#CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	L95	U95	T_STAT	P	ERRCODE
	# added OMMITED and FRQ, which were missing in the original code, but we don't need them, so we drop them too
	rename!(GWASout,
		[:CHR, :BP, :SNP, :A2, :A1, :A1_dup, :A1_again, :OMITTED, :FRQ, :TEST, :NMISS,
		:BETA, :SE, :L95, :U95, :STAT, :P, :ERRCODE]
	)
	select!(GWASout, Not(:A1_dup, :A1_again, :OMITTED, :FRQ)) 

	CSV.write(outdir * ".sumstat", DataFrame(GWASout), delim = "\t")
end

# Converts to Float64; if the number underflows (0.0), sets it to the minimum positive Float64 (1e-323)
parse_safe_p(val) = begin
    bf = parse(BigFloat, string(val))
    f = Float64(bf)
    return f == 0.0 ? nextfloat(0.0) : f
end

function plot_gwas(gwas_out_prefix, outdir)
    df = CSV.File(gwas_out_prefix * ".sumstat", normalizenames=true) |> DataFrame
    
    # 1. Filter out NA and missing values
	df = df[(df.P .!= "NA") .& .!ismissing.(df.P), :]
	
    # 2. Safely parse the P-value column
    df.P = parse_safe_p.(df.P)

    qq_out = @sprintf("%s.qq.png", outdir)
    man_out = @sprintf("%s.man.png", outdir)
    
    @info "Plot qq plot"
    qq(df.P, dpi = 400, fontsize = 10pt, ystep = 5, outfile = qq_out)
    
    @info "Plot manhattan plot"
    manhattan(df.P, df.CHR, df.BP, dpi = 400, fontsize = 10pt, ystep = 5, outfile = man_out)
end


"""Run GWAS for each phenotypic trait and plot results
"""
function run_gwas(ntraits, plink2, covar, synt_data_prefix, eval_dir, all_chromosomes)
    for trait_idx in 1:ntraits
        if all_chromosomes == "all"
            # Aggregate GWAS results across all chromosomes for a single plot
            combined_df = DataFrame()
            # Ensure we don't double-append a chromosome suffix like "-1-1"
            base_prefix = replace(synt_data_prefix, r"-\d+$" => "")
            # Build the base eval_dir for the aggregated output, replacing the
            # trailing "-<number>" chromosome suffix with "-all" instead of
            # stripping it, so the merged file is clearly labeled
            base_eval_dir = replace(eval_dir, r"-\d+$" => "-all")
            for chr in 1:22
                # Genotype prefix requires the "-<number>" suffix; phenotype prefix will have it removed internally
                syngeno_prefix_chr = @sprintf("%s-%d", base_prefix, chr)
                synpheno_prefix_chr = syngeno_prefix_chr
                # Ensure output path reflects the correct chromosome and avoid extra ".chr<j>" suffix
                eval_dir_chr = replace(eval_dir, r"-\d+$" => @sprintf("-%d", chr))
                outdir_chr = @sprintf("%s.trait%i", eval_dir_chr, trait_idx)
                covar_chr = isa(covar, Vector{String}) ? covar[chr] : covar
                run_gwas_tools(plink2, syngeno_prefix_chr, synpheno_prefix_chr, trait_idx, covar_chr, outdir_chr)
                # Read per-chromosome summary stats and append
                chr_sum = CSV.File(outdir_chr * ".sumstat", normalizenames=true) |> DataFrame
                combined_df = isempty(combined_df) ? chr_sum : vcat(combined_df, chr_sum; cols = :union)
            end
            # Write combined summary stats and plot once across all chromosomes
            # Use base_eval_dir (with "-all" suffix) instead of the raw eval_dir,
            # and drop the extra ".all" from the filename since "all" is now
            # already part of the prefix
            outdir_all = @sprintf("%s.trait%i", base_eval_dir, trait_idx)
            CSV.write(outdir_all * ".sumstat", DataFrame(combined_df), delim = "\t")
            plot_gwas(outdir_all, outdir_all)
        else
            # Single chromosome or pre-merged input
            outdir = @sprintf("%s.trait%i", eval_dir, trait_idx)
            run_gwas_tools(plink2, synt_data_prefix, synt_data_prefix, trait_idx, covar, outdir)
            plot_gwas(outdir, outdir)
        end
    end
end
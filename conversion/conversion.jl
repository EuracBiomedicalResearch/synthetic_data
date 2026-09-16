function run_vcf_conversion(options)
    # path of plink output and vcf files needed 
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    
    if chromosome == "all"
        for chromosome_i in 1:22
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            metadata = parse_genomic_metadata(options, superpopulation, fp)
            convert_plink_to_vcf(fp, metadata.plink)

        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        metadata = parse_genomic_metadata(options, superpopulation, fp)

        convert_plink_to_vcf(fp, metadata.plink)
    end
end

function convert_plink_to_vcf(fp, plink)
    input_filename = fp.synthetic_data_prefix
    output_filename = fp.vcf_output_prefix
    run(`$plink --bfile $input_filename --recode vcf --out $output_filename`)
end

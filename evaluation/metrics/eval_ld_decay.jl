using CategoricalArrays, Distances


"""
Function that implements the LD decay metric for ABC
"""
function LD_decay(plink_file, plink, mapthin, out_prefix, bp_to_cm_map)
    # thin snp list to speed up calculations
    # keep 20 snps every 10^6 bases and create a new bim
    bim_file = @sprintf("%s.bim", plink_file)
    snp_thin = @sprintf("%s_snp_thin", out_prefix)
    run(`$mapthin -n -b 20 $bim_file $snp_thin`)

    # compute r2 values between all pairs
    run(`$plink --bfile $plink_file --extract $snp_thin --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out $snp_thin`)

    # extract the distance between each pair and the r2 value
    ld_df = DataFrame(CSV.File(@sprintf("%s.ld", snp_thin), delim=' ', ignorerepeated=true))

    # convert bp to cm
    ld_df.CM_A = [bp_to_cm_map[x] for x in ld_df.BP_A]
    ld_df.CM_B = [bp_to_cm_map[x] for x in ld_df.BP_B]
    ld_df.dist = ld_df.CM_B - ld_df.CM_A
    
    # categorise distances into intervals of a fixed length and compute mean r2 within each interval
    interval = 0.1
    max_dist = maximum(ld_df.dist)
    breaks = collect(0:interval:max_dist)
    
    # Manually assign bins and calculate means
    bin_means = Float64[]
    bin_centers = Float64[]
    
    for i in 1:(length(breaks)-1)
        mask = (ld_df.dist .>= breaks[i]) .& (ld_df.dist .< breaks[i+1])
        if sum(mask) > 0
            push!(bin_means, mean(ld_df.R2[mask]))
            push!(bin_centers, (breaks[i] + breaks[i+1]) / 2)
        end
    end

    return [bin_centers[1:end] bin_means[1:end]]

end

function run_ld_decay(synfile, realfile, plink, mapthin, out_prefix, bp_to_cm_map, chromosome)
    ld_decay_real = LD_decay(realfile, plink, mapthin, out_prefix, bp_to_cm_map)
    ld_decay_syn = LD_decay(synfile, plink, mapthin, out_prefix, bp_to_cm_map)
    # store data points
    df = DataFrame(real_x=ld_decay_real[:,1], real_y=ld_decay_real[:,2], syn_x=ld_decay_syn[:,1], syn_y=ld_decay_syn[:,2])
    outfile = joinpath(dirname(out_prefix), @sprintf("results-chr-%s-ld-decay.csv", chromosome))
    CSV.write(outfile, df)
    @info "LD decay data saved at $outfile"
    # make plot
    outfile = joinpath(dirname(out_prefix), @sprintf("results-chr-%s-ld-decay.png", chromosome))
    fig = Plots.plot(size=(400, 400))
    plot!(fig, ld_decay_real[:,1], ld_decay_real[:,2], label="real", xaxis="Genetic distance (cm)", yaxis="LD estimate (r2)", title=@sprintf("LD decay for chromosome %s", chromosome))
    plot!(fig, ld_decay_syn[:,1], ld_decay_syn[:,2], label="synthetic")
    savefig(fig, outfile)
    @info "LD decay plot saved at $outfile"
    ld_distance = evaluate(Euclidean(), ld_decay_real[:,2], ld_decay_syn[:,2])
    @info @sprintf("LD decay distance between real and synthetic data is %f", ld_distance)
end

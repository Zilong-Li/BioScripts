#!/usr/bin/env julia

using ArgParse
using Printf
using StatsPlots
using CodecZlib
using DelimitedFiles

# disable calling display window
ENV["GKSwstype"] = "100"

function qqfly(out::String, title::String, cutoff::Float64=1e-4, bins::Int=1000)
    allbins = [Float64[-1, 0] for i in 1:bins]
    sigp = Float64[]
    N = 0
    c = -log10(cutoff)
    while !eof(stdin)
        N += 1
        if N % 1000000 == 0
            @printf "reaching the %.6e line.\n" N
        end
        
        p = -log10(parse(Float64, readline(stdin)))
        if p > c
            sigp = [sigp; p]
        else
            i = bins - Int(floor(p / c * bins))
            i == 0 ? i = 1 : i
            allbins[i][1] = max(allbins[i][1], p)
            allbins[i][2] += 1
        end
    end

    size = length(sigp)
    obs = sort(sigp, rev=true)
    exp = [-log10((i + 1 - 0.5) / N) for i in 1:size]
    for b in allbins
        if b[1] != -1
            obs = [obs; b[1]]
            exp = [exp; -log10((b[2] + size - 0.5) / N)]
            size += b[2]
        end
    end
    # plotting
    gr(size=(1000, 1000))
    qqplot(exp, obs, xlabel="Expected -log10(P)", ylabel="Observed -log10(P)", markerstrokecolor=:deepskyblue3, linecolor=:red, markercolor=:deepskyblue3, linewidth=3, title=title)
    figname = out * ".png"
    savefig(figname)
    csvname = out * ".csv.gz"
    open(GzipCompressorStream, csvname, "w") do stream
        writedlm(stream, [obs exp], ',')
    end

    return nothing
end

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--out"
            help = "prefix of output"
            default = "qqplot-jl"
        "--title"
            help = "title of the plot"
        "--cutoff"
            help = "pval bigger than this cutoff will be grouped into bins.[1e-4]"
            arg_type = Float64
            default = 1e-4
        "--bins"
            help = "the number of bins to use.[1000]"
            arg_type = Int
            default = 1000
    end

    args = parse_args(s)
    println("Parsed args:")
    for (arg, val) in args
        println("$arg  =>  $val")
    end

    qqfly(args["out"], args["title"], args["cutoff"], args["bins"])

end

main()


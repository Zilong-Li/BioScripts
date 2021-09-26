#!/usr/bin/env julia

using ArgParse
using StatsPlots
using CodecZlib
using DelimitedFiles
using Dates

# disable calling display window
ENV["GKSwstype"] = "100"

function qqfly(out::String, title::String="", cutoff::Float64=1e-4, bins::Int=1000; cmd=nothing)
    println(Dates.format(now(), "YY-mm-dd:HH:MM:SS  "), "program started.")
    allbins = [Float64[-1, 0] for i in 1:bins]
    sigp = Float64[]
    N = 0
    c = -log10(cutoff)
    if cmd == nothing
        io = stdin
    else
        io = open(cmd)
    end
    while !eof(io)
        N += 1
        if N % 100000000 == 0
            println(Dates.format(now(), "YY-mm-dd:HH:MM:SS  "), "reached the $N line.")
        end

        p = -log10(parse(Float64, readline(io)))
        if p > c
            sigp = [sigp; p]
        else
            i = bins - Int(floor(p / c * bins))
            i == 0 ? i = 1 : i
            allbins[i][1] = max(allbins[i][1], p)
            allbins[i][2] += 1
        end
    end
    println(Dates.format(now(), "YY-mm-dd:HH:MM:SS  "), "total number of lines is $N.")

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
    println(Dates.format(now(), "YY-mm-dd:HH:MM:SS  "), "program finished.")

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
            default = ""
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

# running on terminal with one-liner
# zcat pval.gz | awk 'NR>1{print $10}'| qqplot.jl --out test --title "QQ plot on the fly"
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
# the recommanded is running in julia REPL as follows
# cmd = pipeline(`zcat pval.gz`, `awk 'NR>1{print $10}'`)
# qqfly("test.out", cmd=cmd)

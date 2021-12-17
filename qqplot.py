#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from time import localtime, strftime


def _plot(x, y, ax=None, color=None, ablinecolor="r", alpha=0.8, **kwargs):
    """
    Parameters
    ----------
    x, y : array-like, x is the expected, y is the observed

    ax : matplotlib axis, optional

    color : the dots color in the plot, optional

    kwargs : key, value pairings passed to ``plt.scatter()``

    Returns
    -------
    ax : matplotlib Axes object with the plot.
    """
    # Draw the plot and return the Axes
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")

    # Get the color from the current color cycle
    if color is None:
        (line,) = ax.plot(0, x.mean())
        color = line.get_color()
        line.remove()

    ax.scatter(x, y, c=color, alpha=alpha, edgecolors="none", **kwargs)
    ax.set_xlim(xmin=x.min(), xmax=1.05 * x.max())
    ax.set_ylim(ymin=y.min())

    if ablinecolor:
        # plot the y=x line
        ax.plot(
            [x.min(), ax.get_xlim()[1]],
            [x.min(), ax.get_xlim()[1]],
            color=ablinecolor,
            linestyle="-",
        )

    return ax


def _group_bins(cutoff, bins):

    if cutoff < 0 or cutoff > 1:
        raise ValueError("`cutoff` must be a float in (0, 1).")
    # bin: [pval, size]
    allbins = [[-1, 0] for _ in range(bins)]
    # keep all the low pvals, which may be significant
    sigp = []
    N = 0
    c = -np.log10(cutoff)
    sys.stdout.write(
        strftime("%a, %d %b %Y %H:%M:%S", localtime())
        + f" : start to process the input!\n"
    )
    try:
        for line in sys.stdin:
            N += 1
            if N % 100000000 == 0:
                sys.stdout.write(
                    strftime("%a, %d %b %Y %H:%M:%S", localtime())
                    + f" : reaching the {N}-th line.\n"
                )
            p = -np.log10(float(line))
            if p > c:
                sigp.append(p)
            else:
                # find the bin where the P should go
                idx = bins - int(np.floor(p / c * bins)) - 1
                idx = 0 if idx == -1 else idx
                allbins[idx][0] = max(allbins[idx][0], p)
                allbins[idx][1] += 1
    except KeyboardInterrupt:
        pass

    sys.stdout.write(
        strftime("%a, %d %b %Y %H:%M:%S", localtime())
        + f" : reached the end! the total number of lines is {N}!\n"
    )
    size = len(sigp)
    obs = sorted(sigp, reverse=True)
    exp = [-np.log10((i + 1 - 0.5) / N) for i in np.arange(size)]

    for bi in allbins:
        if bi[0] != -1:
            obs.append(bi[0])
            exp.append(-np.log10((bi[1] + size - 0.5) / N))
            size += bi[1]
        else:
            pass

    return np.array(obs), np.array(exp)


def qq(
    x=None,
    y=None,
    figname=None,
    cutoff=1e-4,
    bins=1000,
    ax=None,
    title=None,
    color=None,
    alpha=0.8,
    ablinecolor="r",
    dpi=300,
    xlabel=None,
    ylabel=None,
    **kwargs,
):

    if x is None:
        # read data from stdin
        o, e = _group_bins(cutoff, bins)
    else:
        if y is None:
            n = len(x)
            a = 0.5
            o = -np.log10(sorted(x))
            e = -np.log10((np.arange(n, dtype=float) + 1 - a) / (n + 1 - 2 * a))
        else:
            o = x
            e = y
    ax = _plot(e, o, ax=ax, color=color, ablinecolor=ablinecolor, alpha=alpha, **kwargs)
    if xlabel is None:
        xlabel = r"Expected $-log_{10}{(P)}$"
    if ylabel is None:
        ylabel = r"Observed $-log_{10}{(P)}$"
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(figname, bbox_inches="tight", dpi=dpi)
    plt.clf()
    plt.close()

    return ax


def main():

    parser = argparse.ArgumentParser(
        description="QQ plot on the fly! only read data from the pipe!"
    )
    parser.add_argument(
        "-cutoff",
        metavar="FLOAT",
        nargs="?",
        default=1e-4,
        type=float,
        help="pval bigger than this cutoff will be grouped into bins. [%(default)s]",
    )
    parser.add_argument(
        "-bins",
        metavar="INT",
        nargs="?",
        default=1000,
        type=int,
        help="the number of bins to use. [%(default)s]",
    )
    parser.add_argument(
        "-csv",
        default=False,
        action="store_true",
        dest="csv",
        help="read csv (two columns) from stdin. [%(default)s]",
    )
    parser.add_argument(
        "-no",
        default=False,
        action="store_true",
        dest="no",
        help="no grouping. [%(default)s]",
    )
    parser.add_argument("-out", metavar="FILE", help="prefix of output files")
    parser.add_argument("-title", metavar="STRING", help="title of the plot")
    args = parser.parse_args()

    assert args.out is not None, "please specify the file of output."

    f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    if args.csv:

        import pandas as pd

        d = pd.read_csv(sys.stdin, header=0, names=["o", "e"])
        qq(
            x=d["o"].values,
            y=d["e"].values,
            figname=args.out + ".png",
            title=args.title,
            ax=ax,
        )
    elif args.no:
        o = list(map(float, filter(None, (row.rstrip() for row in sys.stdin))))
        qq(
            x = o,
            figname=args.out + ".png",
            cutoff=args.cutoff,
            bins=args.bins,
            title=args.title,
            ax=ax,
        )
    else:
        qq(
            figname=args.out + ".png",
            cutoff=args.cutoff,
            bins=args.bins,
            title=args.title,
            ax=ax,
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script estimates the sex of individuals based on the ratio of the counts
of reads mapped on Y and X chromosomes. The script uses the alignment
summary statistics reported by samtools idxstats.

The files of alignment summary statistics are assumed to have the extension
of ".txt". If the files have different extensions, please specify that by
using --ext option. If files have no extension, specify that by `--ext ""`.
"""

import os
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

__author__ = "ITO Tsuyoshi"
__version__ = "0.1.0"
__email__ = "ito.tsuyoshi.3a@kyoto-u.ac.jp"
__date__ = "2020-03-24"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    parser.add_argument("-i",
                        "--input",
                        action="store",
                        dest="idxstats_dir",
                        required=True,
                        help="Path to the directory containing idxstats files")
    parser.add_argument("-y",
                        "--ychrom",
                        action="store",
                        dest="ychrom",
                        required=True,
                        help="Name of Y chromosome in reference sequence")
    parser.add_argument("-x",
                        "--xchrom",
                        action="store",
                        dest="xchrom",
                        required=True,
                        help="Name of X chromosome in reference sequence")
    parser.add_argument("-t",
                        "--threshold",
                        action="store",
                        dest="threshold",
                        default=None,
                        type=float,
                        help="Threshold for sexing")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="out_dir",
                        default=".",
                        help="Name of output directory. It is current "
                             "directory by default.")
    parser.add_argument("-p",
                        "--plot",
                        action="store_true",
                        dest="plot",
                        default=False,
                        help="Output scatter plot and histogram, disabled by "
                             "default")
    parser.add_argument("--min_y",
                        action="store",
                        dest="min_y",
                        default=0,
                        type=int,
                        help="Minimum number of Y chromosome reads required "
                             "for sexing")
    parser.add_argument("--min_x",
                        action="store",
                        dest="min_x",
                        default=1,
                        type=int,
                        help="Minimum number of X chromosome reads required "
                             "for sexing")
    parser.add_argument("--two_thresholds",
                        action="store_true",
                        dest="two_thresholds",
                        default=False,
                        help="Two thresholds will be used for sexing. If "
                             "you use this option, you should specify upper "
                             "and lower thresholds.")
    parser.add_argument("--upper_threshold",
                        action="store",
                        dest="upper_threshold",
                        type=float,
                        help="Upper threshold for sexing")
    parser.add_argument("--lower_threshold",
                        action="store",
                        dest="lower_threshold",
                        type=float,
                        help="Lower threshold for sexing")
    parser.add_argument("--ext",
                        action="store",
                        dest="ext",
                        default=".txt",
                        help="Extension")
    parser.add_argument("--nbins",
                        action="store",
                        dest="n_bins",
                        default=100,
                        type=int,
                        help="Number of bins for histogram")
    parser.add_argument("--dpi",
                        action="store",
                        dest="dpi",
                        default=300,
                        type=int,
                        help="DPI for plots")
    parser.add_argument("-v",
                        "--version",
                        action="version",
                        dest="out_file",
                        version=__version__)

    args = parser.parse_args()

    if args.threshold is None:
        print("Threshold was automatically set to be the one fourth of the "
              "ratio of sequence length of Y and X chromosomes. Use -t "
              "option if you want to manually specify threshold.\n")

    col_names = ["id", "sex", "nreads_ychrom", "nreads_xchrom", "ratio",
                 "threshold"]
    data = []

    for file_name in os.listdir(args.idxstats_dir):

        # This avoids reading files with different extension or hidden files
        if not file_name.startswith('.') and \
                os.path.splitext(file_name)[1] == args.ext:

            # This removes extension from file name
            sample_name = os.path.splitext(file_name)[0]

            # This counts the umber of reads mapped on Y and X chromosomes
            nreads_ychrom = None
            nreads_xchrom = None
            with open(args.idxstats_dir + "/" + file_name) as lines:
                for line in lines:
                    if args.ychrom in line:
                        nreads_ychrom = int(line.split('\t')[2])
                        seqlen_ychrom = int(line.split('\t')[1])
                    elif args.xchrom in line:
                        nreads_xchrom = int(line.split('\t')[2])
                        seqlen_xchrom = int(line.split('\t')[1])

            # If the specified name of Y and X chromosome cannot be found in
            # a target file, warning will be raised and sex will be
            # noted as undetermined.
            if nreads_ychrom is None or nreads_xchrom is None:
                print("Warning: Cant't find the specified name of Y or X "
                      "chromosome in " + file_name)
                sex = "undetermined"
                ratio = None
                threshold = None

            else:
                # This is for the option of minimum number of Y and X
                # chromosome reads required for sexing
                if nreads_ychrom < args.min_y or nreads_xchrom < args.min_x:
                    sex = "undetermined"

                else:
                    ratio = float(nreads_ychrom) / float(nreads_xchrom)

                    # This is for the option of two thresholds. Upper and
                    # lower thresholds must be specified.
                    if args.two_thresholds:
                        if args.lower_threshold is None or \
                                args.upper_threshold is None:
                            raise Exception("You must specify "
                                            "--lower_threshold "
                                            "and --upper_threshold when using "
                                            "--two_thresholds option")
                        else:
                            if args.upper_threshold < ratio:
                                sex = "male"
                            elif args.lower_threshold > ratio:
                                sex = "female"
                            else:
                                sex = "undetermined"

                    else:
                        # Threshold was automatically set if it was not
                        # specified by an user
                        if args.threshold is None:
                            threshold = float(seqlen_ychrom) * 0.25 / float(
                                seqlen_xchrom)
                        else:
                            threshold = args.threshold
                        if threshold < ratio:
                            sex = "male"
                        else:
                            sex = "female"

            data += [[sample_name, sex, nreads_ychrom, nreads_xchrom, ratio,
                      threshold]]

    df_sexing = pd.DataFrame(data, columns=col_names)

    # This checks whether thresholds are different among samples,
    # which indicates that sequence length is different among samples
    if len(df_sexing["threshold"].dropna().unique()) == 1:
        df_sexing.drop("threshold", axis=1).to_csv(args.out_dir + "/" +
                                                   "sexing.txt", index=False,
                                                   sep='\t')
        print("threshold = " + str(threshold) + "\n")
    else:
        df_sexing.to_csv(args.out_dir + "/" + "sexing.txt", index=False,
                         sep='\t')
        print("Threshold was differently set among samples")

    print("Results:")
    print(df_sexing["sex"].value_counts().to_frame())

    # This is for the option of plotting figures+
    if args.plot:
        try:
            # Scatter plot of the number of reads mapped on Y and X chromosomes
            sns.scatterplot(x="nreads_xchrom", y="nreads_ychrom", hue="sex",
                            data=df_sexing)
            plt.savefig(args.out_dir + "/scatter_plot.png", dpi=args.dpi)
            plt.close()

            # Histogram of the ratio of the counts of reads mapped on Y and X
            # chromosomes
            sns.distplot(df_sexing["ratio"].dropna(), kde=False,
                         bins=args.n_bins)
            if args.two_thresholds:
                plt.axvline(args.upper_threshold, label="upper_threshold",
                            linewidth=0.5, color=sns.color_palette()[1])
                plt.axvline(args.lower_threshold, label="lower_threshold",
                            linewidth=0.5, color=sns.color_palette()[2])
                plt.legend()
                plt.xlabel("Y/X ratio")
            else:
                # If threshold is different among sample, the threshold line
                # will not be shown
                if len(df_sexing["threshold"].dropna().unique()) == 1:
                    plt.axvline(threshold, linewidth=0.5, label="threshold")
                    plt.legend()
                plt.xlabel("Y/X ratio")
            plt.savefig(args.out_dir + "/hist_plot.png", dpi=args.dpi)
            plt.close()
        except Exception as inst:
            print("Error in plotting: {}".format(inst))


if __name__ == "__main__":
    main()

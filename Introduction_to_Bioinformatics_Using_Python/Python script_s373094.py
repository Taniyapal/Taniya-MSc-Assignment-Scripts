#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 06:56:43 2021

@author: taniyapal
"""
# importing libraries
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# defining the function modules()
def modules(lines, start_text, end_text):
    """For extracting the information from each module in the FASTQC text file"""

    start_count = None
    end_count = None

    for count, line in enumerate(lines):

        if line.startswith(start_text):
            filter_information = line.split("\t")[
                1
            ]  # storing the pass fail or warn information about the file

            start_count = count

        if line == end_text + "\n":
            if start_count is not None:  # so that "start_text" is not an empty string
                end_count = count

                break

    join_line = lines[
        start_count : end_count + 1
    ]  # extracting the string from the list "lines" at the index "[start_count:end_count+1]"
    text = "".join(join_line)
    return (text, start_count, end_count, filter_information)


def main():

    parser = argparse.ArgumentParser(
        description="""This script takes a Fastqc text file and parses the modules from it giving back its corresponding plot."""
    )
    parser.add_argument("--file", help="inputs the Fastqc text file")
    parser.add_argument(
        "--output", help="stores the output of the script in that particular file"
    )
    parser.add_argument(
        "-c",
        "--Basic_Statistics",
        action="store_true",
        help="displays the Basic Statistics module",
    )
    parser.add_argument(
        "-p",
        "--Per_base_sequence-quality",
        action="store_true",
        help="displays the Per Base Sequence quality module in the fast qc text file",
    )
    parser.add_argument(
        "-t",
        "--Per_tile_sequence-quality",
        action="store_true",
        help="displays the Per tile sequence quality module in the fast qc text file",
    )
    parser.add_argument(
        "-s",
        "--Per_sequence_quality-scores",
        action="store_true",
        help="displays the Per sequece quality score module in the fast qc text file",
    )
    parser.add_argument(
        "-b",
        "--Per_base_sequence_content",
        action="store_true",
        help="displays the Per base sequence content module in the fast qc text file",
    )
    parser.add_argument(
        "-g",
        "--Per_sequence_GC_content",
        action="store_true",
        help="displays the Per sequence GC content module in the fast qc text file",
    )
    parser.add_argument(
        "-n",
        "--Per_base_N_content",
        action="store_true",
        help="displays the Per base N content module in the fast qc text file",
    )
    parser.add_argument(
        "-l",
        "--Sequence_length_distribution",
        action="store_true",
        help="displays the sequence length distribution module in the fast qc text file",
    )
    parser.add_argument(
        "-d",
        "--Sequence_duplication_levels",
        action="store_true",
        help="displays the sequence duplication level module in the fast qc text file",
    )
    parser.add_argument(
        "-o",
        "--Overrepresented_Sequences",
        action="store_true",
        help="displays the overrepresented sequences module in the fast qc text file",
    )
    parser.add_argument(
        "-k",
        "--Kmer_content",
        action="store_true",
        help="displays the kmer content module in the fast qc text file",
    )
    parser.add_argument(
        "-r",
        "--Adapter_content",
        action="store_true",
        help="displays the Adapter Content module in the fast qc text file",
    )
    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        help="Prints the information under all modules",
    )
    args = parser.parse_args()
    input_file = args.file
    input_file_open = open(input_file, "r")  # open the input file from the command line
    lines = (
        input_file_open.readlines()
    )  # reading the lines from the file stored in the input_file variable declared in the function main()
    input_file_open.close()
    output_file = args.output
    open(
        output_file, "w"
    ).close()  # so that the output is overwritten everytime the arguments are parsed

    if args.all:
        args.Basic_Statistics = True
        args.Per_base_sequence_quality = True
        args.Per_tile_sequence_quality = True
        args.Per_sequence_quality_scores = True
        args.Per_base_sequence_content = True
        args.Per_sequence_GC_content = True
        args.Per_base_N_content = True
        args.Sequence_length_distribution = True
        args.Sequence_duplication_levels = True
        args.Overrepresented_Sequences = True
        args.Adapter_content = True
        args.Kmer_content = True

    if args.Basic_Statistics:
        module_0, start_count, end_count, filter_information = modules(
            lines, ">>Basic Statistics", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_0 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_0)
        print(filter_info)
        table_0_array = table_0.to_numpy()

        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(
            filter_info + "\n"
        )  # appending the filter information to the output file
        np.savetxt(
            output_file_pointer, table_0_array, fmt="%s"
        )  # saving the output into a file in the format of string

    if args.Per_base_sequence_quality:
        module_1, start_count, end_count, filter_information = modules(
            lines, ">>Per base sequence quality", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_1 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_1)
        print(filter_info)
        row = table_1.loc[
            :, "#Base"
        ]  # extracting the particular column from the dataframe "table_1"
        column = table_1.loc[:, "Mean"]
        plot_1 = sns.barplot(x=row, y=column, data=table_1)
        plot_1.set_title("Per Base Sequence Quality plot")
        plot_1.set_xlabel("Base")
        plot_1.set_ylabel("Mean Quality Score")
        plt.savefig("Per_base_sequence_quality_plot.png")
        plt.close()

        table_1_array = table_1.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_1_array, fmt="%s")

    if args.Per_tile_sequence_quality:
        module_2, start_count, end_count, filter_information = modules(
            lines, ">>Per tile sequence quality", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_2 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_2)
        print(filter_info)
        plot_2 = sns.heatmap(table_2, annot=True)

        plot_2.set_title("Per tile sequence quality plot")
        plt.savefig("Per_tile_sequence_quality_plot.png")
        plt.close()

        table_2_array = table_2.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_2_array, fmt="%s")

    if args.Per_sequence_quality_scores:
        module_3, start_count, end_count, filter_information = modules(
            lines, ">>Per sequence quality scores", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_3 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_3)
        print(filter_info)
        quality = table_3.loc[:, "#Quality"]
        count = table_3.loc[:, "Count"]
        plot_3 = sns.lineplot(x=quality, y=count)
        plot_3.set_xlabel("Base Count", fontsize=12)
        plot_3.set_ylabel("Per sequence quality", fontsize=12)
        plot_3.set_title("Quality score distribution over all sequences")
        plt.savefig("Per_sequence_quality_scores_plot.png")
        plt.close()
        table_3_array = table_3.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_3_array, fmt="%s")

    if args.Per_base_sequence_content:
        module_4, start_count, end_count, filter_information = modules(
            lines, ">>Per base sequence content", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_4 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_4)
        print(filter_info)
        base = table_4.loc[:, "#Base"]
        percent_G = table_4.loc[:, "G"]
        percent_A = table_4.loc[:, "A"]
        percent_T = table_4.loc[:, "T"]
        percent_C = table_4.loc[:, "C"]
        plot_4 = sns.lineplot(x=base, y=percent_G, legend="full", label="G")
        plot_4 = sns.lineplot(x=base, y=percent_A, legend="full", label="A")
        plot_4 = sns.lineplot(x=base, y=percent_T, legend="full", label="T")
        plot_4 = sns.lineplot(x=base, y=percent_C, legend="full", label="C")
        plot_4.set_xlabel("Base Count", fontsize=12)
        plot_4.set_ylabel("Sequence content", fontsize=12)
        plot_4.set_title("Per base sequence content")
        plt.savefig("Per_base_sequence_content_plot.png")
        plt.close()
        table_4_array = table_4.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_4_array, fmt="%s")

    if args.Per_sequence_GC_content:
        module_5, start_count, end_count, filter_information = modules(
            lines, ">>Per sequence GC content", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_5 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_5)
        print(filter_info)
        sample_x = range(35, 75)
        gc_content = table_5.loc[:, "#GC Content"]
        count = table_5.loc[:, "Count"]
        plot_5 = sns.lineplot(
            x=gc_content, y=count, color="red", legend="full", label="GC Count per read"
        )
        plot_5_1 = plot_5.twinx()
        plot_5_1 = sns.distplot(sample_x, color="blue", hist=False)
        plot_5_1.set(yticklabels=[])
        plot_5_1.set(ylabel=None)
        plot_5_1.tick_params(right=False)
        plot_5.set_xlabel("GC Content in each sequence", fontsize=12)
        plot_5.set_ylabel("Number of Sequence", fontsize=12)
        plot_5.set_title("Per Sequence GC Content")
        plt.savefig("Per_Sequence_GC_content_plot.png")
        plt.close()
        table_5_array = table_5.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_5_array, fmt="%s")

    if args.Per_base_N_content:
        module_6, start_count, end_count, filter_information = modules(
            lines, ">>Per base N content", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_6 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        plot_6 = sns.lineplot(x="#Base", y="N-Count", data=table_6)
        print(table_6)
        print(filter_info)
        plot_6.set_xlabel("Base Count")
        plot_6.set_ylabel("%N")
        plot_6.set_title("Per base N content")
        plt.savefig("Per_base_N_content_plot.png")
        plt.close()
        table_6_array = table_6.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_6_array, fmt="%s")

    if args.Sequence_length_distribution:
        module_7, start_count, end_count, filter_information = modules(
            lines, ">>Sequence Length Distribution", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_7 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_7)
        print(filter_info)
        len_data = table_7.loc[:, "#Length"]
        length = [100, len_data, 50]
        plot_7 = sns.distplot(length, hist=False)
        plot_7.set_xlabel("Length of Sequence")
        plot_7.set_title("Sequence Length Distribution plot")
        plot_7.set(yticklabels=[])
        plot_7.set(ylabel=None)
        plot_7.tick_params(left=False)

        plt.savefig("Sequence_Length_Distribution_plot.png")
        plt.close()
        table_7_array = table_7.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_7_array, fmt="%s")

    if args.Sequence_duplication_levels:
        module_8, start_count, end_count, filter_information = modules(
            lines, ">>Sequence Duplication Levels", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_8 = pd.read_csv(
            input_file,
            skiprows=start_count + 2,
            nrows=end_count - (start_count + 3),
            delimiter="\t",
        )
        print(table_8)
        print(filter_info)
        plot_8 = sns.lineplot(
            x="#Duplication Level",
            y="Percentage of deduplicated",
            data=table_8,
            color="red",
            legend="full",
            label="% Deduplicated Sequences",
        )
        plot_8_1 = plot_8.twinx()
        plot_8_1 = sns.lineplot(
            x="#Duplication Level",
            y="Percentage of total",
            data=table_8,
            color="blue",
            legend="auto",
            label="% Total Sequences",
        )
        plot_8_1.set(yticklabels=[])
        plot_8_1.set(ylabel=None)
        plot_8_1.tick_params(right=False)
        plot_8.set_xlabel("Sequence Duplication Level")
        plot_8.set(ylabel=None)
        plot_8.set_title("Sequence Duplication Level plot")
        plt.savefig("Sequence_Duplication_Level_plot.png")
        plt.close()
        table_8_array = table_8.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_8_array, fmt="%s")

    if args.Overrepresented_Sequences:
        module_9, start_count, end_count, filter_information = modules(
            lines, ">>Overrepresented sequences", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_9 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_9)
        print(filter_info)
        table_9_array = table_9.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_9_array, fmt="%s")

    if args.Adapter_content:
        module_10, start_count, end_count, filter_information = modules(
            lines, ">>Adapter Content", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_10 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_10)
        print(filter_info)
        Illumina_univ = table_10.loc[:, "Illumina Universal Adapter"]
        Illumina_small = table_10.loc[:, "Illumina Small RNA Adapter"]
        Nextera = table_10.loc[:, "Nextera Transposase Sequence"]
        solid = table_10.loc[:, "SOLID Small RNA Adapter"]
        row = table_10.loc[:, "#Position"]
        plot_10 = sns.lineplot(
            x=row,
            y=Illumina_univ,
            data=table_10,
            legend="full",
            label="Illumina Universal Adapter",
        )
        plot_10 = sns.lineplot(
            x=row,
            y=Illumina_small,
            data=table_10,
            legend="full",
            label="Illumina Small RNA Adapter",
        )
        plot_10 = sns.lineplot(
            x=row,
            y=Nextera,
            data=table_10,
            legend="full",
            label="Nextera Transposase Sequence",
        )
        plot_10 = sns.lineplot(
            x=row,
            y=solid,
            data=table_10,
            legend="full",
            label="SOLID Small RNA Adapter",
        )
        plot_10.set_xlabel("Position")
        plot_10.set_ylabel("Adapter")
        plot_10.set_title("Adapter content plot")

        plt.savefig("Adapter_Content_plot.png")
        plt.close()
        table_10_array = table_10.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_10_array, fmt="%s")

    if args.Kmer_content:
        module_11, start_count, end_count, filter_information = modules(
            lines, ">>Kmer Content", ">>END_MODULE"
        )
        filter_info = "This module has the filter " + filter_information
        table_11 = pd.read_csv(
            input_file,
            skiprows=start_count + 1,
            nrows=end_count - (start_count + 2),
            delimiter="\t",
        )
        print(table_11)
        print(filter_info)

        row = table_11.loc[:, "#Sequence"]
        column = table_11.loc[:, "Obs/Exp Max"]
        plot_11 = sns.barplot(x=row, y=column, data=table_11)
        plot_11.set_title("Kmer Content Plot")
        plot_11.set_xlabel("Sequence Number")
        plot_11.set_ylabel("Expected/Observed Frequency")
        plt.savefig("Kmer_Content_plot.png")
        plt.close()

        table_11_array = table_11.to_numpy()
        output_file_pointer = open(output_file, "a")
        output_file_pointer.write(filter_info + "\n")
        np.savetxt(output_file_pointer, table_11_array, fmt="%s")

    input_file_open.close()  # closing the input file given in the command file
    output_file_pointer.close()  # closing the ouput file generated


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot genome depth with nucleotide polymorphisms"""

import argparse
import os
import sys
import numpy as np
import pandas as pd


from matplotlib import pyplot as plt

def get_args():
    """Initiate arguments"""
    parser = argparse.ArgumentParser(
            description='Plot genome depth with nucleotide polymorphisms',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
            'vcf_file',
            help='VCF file with AD format',
            default=sys.stdin,
            type=str,
            )

    def float01(x):
        try:
            x = float(x)
        except ValueError as error:
            raise argparse.ArgumentTypeError('%r must be float' % (x,)) from error

        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError('%r is out of range, must be in [0.0, 1.0]' % (x,))
        return x

    parser.add_argument(
            '-M',
            '--major',
            help='Threshold to include variant at position',
            type=float01,
            metavar='FLOAT',
            default=0.1,
            )

    parser.add_argument(
            '-m',
            '--minor',
            help='Threshold to include base variant',
            type=float01,
            metavar='FLOAT',
            default=0.1,
            )

    parser.add_argument(
            '-e',
            '--extend',
            help='Extend variant bar for INT position to both sides',
            type=int,
            metavar='INT',
            default=4,
            )

    parser.add_argument(
            '-f',
            '--figure',
            help='Output figure',
            type=str,
            metavar='STR',
            )

    parser.add_argument(
            '-F',
            '--force',
            help='Force overwite output',
            action='store_true',
            )

    parser.add_argument(
            '-l',
            '--log',
            help='Plot depth in log scale',
            action='store_true'
            )

    parser.add_argument(
            '-D',
            '--dpi',
            help='Image DPI',
            type=int,
            metavar='INT',
            default=300
            )

    parser.add_argument(
            '-s',
            '--size',
            help='Size of output',
            metavar='FLOAT,FLOAT',
            default='16,9',
            )

    parser.add_argument(
            '-g',
            '--gap',
            help='Gap position and size',
            type=str,
            metavar='INT,INT[,INT,INT[,INT,INT,...]]',
            )

    parser.add_argument(
            '-x',
            '--x-tick',
            help='Tick interval',
            type=int,
            default=500,
            )

    parser.add_argument(
            '-t',
            '--table-relative',
            help="Save relative variant table to STR. Use '' for STDOUT",
            type=str,
            metavar='STR',
            )

    parser.add_argument(
            '-T',
            '--table-absolute',
            help="Save aboslute variant table to STR. Use '' for STDOUT",
            type=str,
            metavar='STR',
            )

    parser.add_argument(
            '-d',
            '--depth',
            help='Depth of position to report variant',
            type=int,
            metavar='INT',
            default=20,
            )

    args = parser.parse_args()

    # Check input
    if not os.path.isfile(args.vcf_file):
        parser.error('Input not exist')

    # Check output
    if args.figure:
        if os.path.isfile(args.figure) and not args.force:
            parser.error(f'{args.figure} exist. Use flag -F to overwrite.')

    if args.table_relative:
        if os.path.isfile(args.table_relative) and not args.force:
            parser.error(f'{args.table_relative} exist. Use flag -F to overwrite.')

    if args.table_relative == '':
        args.table_relative = sys.stdout

    if args.table_absolute:
        if os.path.isfile(args.table_absolute) and not args.force:
            parser.error(f'{args.table_absolute} exist. Use flag -F to overwrite.')

    if args.table_absolute == '':
        args.table_absolute = sys.stdout

    # Check size
    try:
        fig_size = args.size
        fig_size = fig_size.split(',')
        if len(fig_size) != 2:
            parser.error('Size must contain 2 integers and seperate by ","')
        fig_size = [float(x) for x in fig_size]
        args.size = fig_size
    except ValueError:
        parser.error('Size must be integer')

    # Check gap
    if args.gap:
        try:
            gap = args.gap
            gap = gap.split(',')
            gap = [float(x) for x in gap]
            gap_decimal = [x%1 for x in gap]
            if sum(gap_decimal):
                parser.error('Gap must be integer')
            if len(gap)%2 != 0:
                parser.error('Gap size missing')
            gap = [int(x) for x in gap]
            gap = [[gap[x], gap[x+1]] for x in range(0, len(gap), 2)]
            args.gap = gap
        except ValueError:
            parser.error('Gap must be integer')

    # Check depth
    if args.depth < 0:
        parser.error('Depth argument must not be negative')

    return args


def parse_vcf(vcf_file):
    """
    Read VCF file and return vcf dataframe
    Extract pos, ref, alt, and DP4
    """
    vcf_lines = [line.rstrip().split() for line in vcf_file if line[0] != '#' and "INDEL" not in line]
    pos = [int(line[1]) for line in vcf_lines]
    refs = [line[3] for line in vcf_lines]
    alts = [line[4].split(',') for line in vcf_lines] # nested list
    info = [dict(
            [j if len(j)==2 else j+[""] for j in map(lambda i: str.split(i, "="), line[7].split(";"))]
            ) for line in vcf_lines
            ]
    af = [line["AF"] for line in info]
    vcf_df = pd.DataFrame({
            'pos': pos,
            'ref': refs,
            'alt': alts,
            'af': af,
            })
    vcf_df["af"] = pd.to_numeric(vcf_df["af"])
    return vcf_df


def add_gaps(depth_df, gaps):
    """
    Add regions with depth = 0 from a list of gaps and return updated DataFrame
    Parameters
    ----------
    depth_df    DataFrame with ["id", "pos", "depth", "ref", "alt", "af"]
    gaps        List of [position, gap_length]. Ex. [[1, 30], [1024, 600]]
    """
    for gap_pos, gap_len, in gaps:
        depth_df["pos"] = depth_df.apply(
                lambda i: i["pos"] + gap_len if i["pos"] >= gap_pos else i["pos"],
                axis=1
                )

    gap_insert = pd.DataFrame(
            {
                    "id": [merged_depth_df["id"][0]]*gap_len,
                    "pos": [gap_pos+i for i in range(gap_len)],
                    "depth":[np.nan]*gap_len
                    }
            )
    merged_depth_df = pd.concat([merged_depth_df, gap_insert], ignore_index=True)
    merged_depth_df = merged_depth_df.sort_values(by="pos", ignore_index=True)
    merged_depth_df.reset_index()
    return depth_df


def make_variant_df(depth_df, extend=4):
    """Return two df with true variant depth and extended variant depth"""
    variant_list = []
    for index, row in depth_df.iterrows():
        row = list(row)
        row = row[:-1]
        if row[-1]:
            row[-1] = 100
            variant_list.append(row)
        else:
            base_sum = sum(row[1:5])
            entry = [row[0]]
            for base in row[1:5]:
                relative_depth = base/base_sum*100 if base_sum > 0 else 0
                entry.append(relative_depth)
            entry.append(0)
            variant_list.append(entry)

    extended_variant_list = list(variant_list)
    if extend:
        takeover_pos = []
        for i, entry in enumerate(extended_variant_list):
            if sum(entry[1:5]):
                takeover_pos.append(i)
        takeover_entry = [extended_variant_list[x] for x in takeover_pos]
        for pos, entry in zip(takeover_pos, takeover_entry):
            for i in range(-extend,extend+1):
                extended_variant_list[pos+i] = entry

    variant_df = pd.DataFrame(
            variant_list,
            columns=['pos', 'A', 'T', 'C', 'G', 'noVar'],
            )

    extended_variant_df = pd.DataFrame(
            extended_variant_list,
            columns=['pos', 'A', 'T', 'C', 'G', 'noVar'],
            )

    return variant_df, extended_variant_df


def main():
    """Plot please"""
    args = get_args()
    with open(args.vcf_file) as vcf_file:
        vcf_df = parse_vcf(vcf_file)
    deep_vcf_df = vcf_df[vcf_df["af"] >= args.threshold]
    depth_df = pd.read_csv(
            args.depth_file,
            sep="\t",
            names=["id", "pos", "depth"]
            )
    depth_df = depth_df.merge(deep_vcf_df, how="left")
    if args.gaps:
        depth_df = add_gaps(depth_df, args.gaps)
    if args.extend:
        depth_df = add_extension(depth_df, args.extend)
    variant_df, extended_variant_df = make_variant_df(depth_df, args.extend)
    if args.figure:
        colors = ['#5772B2', '#3A9276', '#F0430F', '#B615D6', '#DEE0E3']
        # fig = plt.figure(figsize=args.size)
        plt.rcParams["figure.dpi"] = args.dpi
        ax = extended_variant_df.plot.bar(
                x='pos',
                stacked=True,
                color=colors,
                figsize=args.size, #(10.5,0.75),
                rot=30,
                width=1,
                ylim=[0,100],
                yticks=[0, 50, 100],
                xticks=list(range(0, len(extended_variant_df), args.x_tick)),
                # legend=False,
                )
        ax2 = ax.twinx()
        depth_df.plot.line(
                x='pos',
                y='sum_depth',
                lw=0.5,
                secondary_y=True,
                # ylim=[1, max(list(depth_df['sum_depth']).remove(np.inf))],
                logy=args.log,
                color=['black'],
                legend=False,
                # xticks=list(range(-1, len(depth_df), 500)) + [0],
                ax=ax2,
                )
        ax.ticklabel_format(axis='x', style='plain')
        ax.minorticks_off()
        ax.figure.savefig(args.figure, bbox_inches='tight')

    if args.table_relative:
        relative_variant_df = variant_df[
                (variant_df['A'] > 0) | # or
                (variant_df['T'] > 0) | # or
                (variant_df['C'] > 0) | # or
                (variant_df['G'] > 0)
                ]
        relative_variant_df.to_csv(
                args.table_relative,
                columns=['pos', 'A', 'T', 'C', 'G'],
                index=False,
                )

    if args.table_absolute:
        absolute_variant_df = depth_df[
                (depth_df['A'] > 0) |
                (depth_df['T'] > 0) |
                (depth_df['C'] > 0) |
                (depth_df['G'] > 0)
                ]
        absolute_variant_df.to_csv(
                args.table_absolute,
                columns=['pos', 'A', 'T', 'C', 'G'],
                index=False,
                )

if __name__ == "__main__":
    main()

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
    Extract only pos, ref, alt, and AD in info column
    """
    vcf_lines = [line.rstrip().split() for line in vcf_file if line[0] != '#']
    pos = [int(line[1]) for line in vcf_lines]
    refs = [line[3] for line in vcf_lines]
    alts = [line[4].split(',') for line in vcf_lines] # nested list`
    ad_pos = vcf_lines[0][8].split(':')
    ad_pos = ad_pos.index('AD')
    ads = [
            [int(x) for x in line[9].split(':')[ad_pos].split(',')[:-1]]
            for line in vcf_lines] # nested list
    vcf_df = pd.DataFrame({
            'pos': pos,
            'ref': refs,
            'alt': alts,
            'ads': ads,
            })
    return vcf_df


def find_variants(freq_list, major_threshold=0.2, depth=20):
    """Return Boolean value if proportion of variant >= threshold"""
    total = sum(freq_list)
    if total < depth:
        return False
    variant_sum = sum(freq_list[1:])
    return variant_sum/total > major_threshold # make parameter


def make_depth_df(vcf_df, minor_threshold=0.1, gaps=None):
    """Return df with reference depth and alternative depth"""
    depth_list = [] # [A, T, C, G, noVar] for pos x
    base_dict = {
            'A': 1,
            'T': 2,
            'C': 3,
            'G': 4,
            'noVar': 5,
            }
    for index, row in vcf_df.iterrows():
        if row['has_variant']:
            entry = [int(row['pos']), 0, 0, 0, 0, 0]
            base = row['ref']
            base_index = base_dict[base]
            depth = row['ads'][0]
            entry[base_index] = depth
            for i in range(len(row['ads'][1:])):
                depth = row['ads'][i+1]
                if depth/sum(row['ads']) > minor_threshold: # make parameter
                    base = row['alt'][i]
                    base_index = base_dict[base]
                    entry[base_index] = depth
        else:
            entry = [int(row['pos']), 0, 0, 0, 0, row['ads'][0]]
        entry.append(sum(row['ads']))
        depth_list.append(entry)

    if gaps:
        gaps.sort(key=lambda x: x[0], reverse=True)
        for pos, length in gaps:
            for entry in depth_list:
                if entry[0] >= pos:
                    entry[0] += length
            gap_entries = [[x, 0, 0, 0, 0, 0, np.nan] for x in range(pos, pos+length)]
            depth_list = depth_list[:pos-1] + gap_entries + depth_list[pos-1:]

    depth_df = pd.DataFrame(
           depth_list,
           columns=['pos', 'A', 'T', 'C', 'G', 'noVar', 'sum_depth'],
           )

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
    vcf_df['has_variant'] = vcf_df.apply(
            lambda x: find_variants(x['ads'], args.major, args.depth), axis=1
            )
    depth_df = make_depth_df(vcf_df, minor_threshold=args.minor, gaps=args.gap)
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

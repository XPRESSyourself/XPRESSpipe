"""
XPRESSpipe
An alignment and analysis pipeline for RNAseq data
alias: xpresspipe

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""IMPORT DEPENDENCIES"""
import os
import sys
import pandas as pd
import numpy as np
from math import ceil
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import add_directory

input_count = 'coverage'
window = 20

"""Compile images from a list of metrics matrices"""
def compile_matrix_metrics(
        path,
        file_list,
        lab_x,
        lab_y,
        plot_type,
        experiment,
        plot_output,
        dpi=600):

    # Keep axes happy to avoid 'IndexError: too many indices for array' error
    # Auto formats figure summary size based on number of plots
    if (len(file_list) / 2) < 2:
        plot_rows = 2
        fig_size = (15, 16)
    else:
        plot_rows = ceil(len(file_list) / 2)
        fig_size = (15, (8 * (int(len(file_list) / 2))))

    # Set up figure space
    fig, axes = plt.subplots(
        nrows = plot_rows,
        ncols = 2,
        figsize = fig_size,
        subplot_kw = {'facecolor':'none'})
    plt.subplots_adjust(
        bottom = 0.3)

    # Initialize file and axis counters for formatting summary figure
    file_number = 0
    ax_y = 0

    for file in file_list:
        x = 0
        df = pd.read_csv(
            str(path) + str(file),
            sep = '\t',
            index_col = 0) # Initialize dataframe for relevant data
        df = df.dropna(
            axis = 0,
            subset = [str(lab_x), str(lab_y)]) # Remove rows where pertinent information is missing
        df = df.sort_index()

        # Prepare subplots
        if (file_number % 2) == 0:
            ax_x = 0
        else:
            ax_x = 1

        if file_number != 0:
            if (file_number % 2) == 0:
                ax_y += 1

        # Plot figure
        if 'read_distribution' in plot_type:
            df.plot.bar(
                x = lab_x,
                y = lab_y,
                title = file.rsplit('.',1)[0],
                ax = axes[ax_y, ax_x],
                grid = None,
                width = 0.8)
        else:
            df.plot.line(
                x = lab_x,
                y = lab_y,
                title = file.rsplit('.',1)[0],
                ax = axes[ax_y, ax_x],
                grid = None)

        axes[ax_y, ax_x].axhline(
            0,
            xmin = 0,
            ls = '-',
            color = 'black')

        file_number += 1

    # Save catenated figure
    plot_title = str(experiment) + '_' + str(plot_type) # Make figure title to save as from experiment name and plot type
    fig.savefig(
        str(plot_output) + plot_title + '_summary.pdf',
        dpi = dpi,
        bbox_inches = 'tight')


"""
"""
def compile_periodicity_metrics(
        path,
        file_list,
        plot_type,
        experiment,
        plot_output,
        dpi=600):

    # Keep axes happy to avoid 'IndexError: too many indices for array' error
    # Auto formats figure summary size based on number of plots
    if (len(file_list) / 2) < 2:
        plot_rows = 2
        fig_size = (15, 12)
    else:
        plot_rows = len(file_list)
        fig_size = (15, (6 * (int(len(file_list)))))

    # Set up figure space
    fig, axes = plt.subplots(
        nrows = plot_rows,
        ncols = 2,
        figsize = fig_size,
        sharey = True,
        subplot_kw = {'facecolor':'none'})

    tix_5prime = ['']
    for t in range(-24,76,3):
        if t == 0:
            tix_5prime.append('START')
            tix_5prime.append('')
            tix_5prime.append('')
        elif t == 75:
            tix_5prime.append(str(t))
        else:
            tix_5prime.append(str(t))
            tix_5prime.append('')
            tix_5prime.append('')

    tix_3prime = []
    for t in range(-75,26,3):
        if t == 0:
            tix_3prime.append('STOP')
            tix_3prime.append('')
            tix_3prime.append('')
        elif t == 24:
            tix_3prime.append(str(t))
            tix_3prime.append('')
        else:
            tix_3prime.append(str(t))
            tix_3prime.append('')
            tix_3prime.append('')

    # Initialize file and axis counters for formatting summary figure
    ax_y = 0

    for file in file_list:

        # Get data
        df = pd.read_csv(
            str(path) + str(file),
            sep = '\t',
            index_col = 0) # Initialize dataframe for relevant data
        df = df.dropna(axis = 0)
        df['distance'] = df['distance'].astype(int)

        # Get start data and fill missing points
        df_start = df[df['reg'].str.contains('start')]
        counter = max(df_start.index.tolist())
        indexer = df_start['distance'].tolist()
        for i in range(-25, 76):
            if i not in indexer:
                counter += 1
                df_start.loc[counter] = [i,0.2,'start']
        df_start = df_start.sort_values('distance')
        df_start = df_start.reset_index(drop=True)

        # Get stop data and fill missing points
        df_stop = df[df['reg'].str.contains('stop')]
        counter = max(df_stop.index.tolist())
        indexer = df_stop['distance'].tolist()
        for i in range(-75, 26):
            if i not in indexer:
                counter += 1
                df_stop.loc[counter] = [i,0.2,'stop']
        df_stop = df_stop.sort_values('distance')
        df_stop = df_stop.reset_index(drop=True)

        # Plot 5prime figure
        df_start.plot.bar(
            x = 'distance',
            y = 'reads',
            ax = axes[ax_y, 0],
            width = 0.9)
        axes[ax_y, 0].set_xticklabels(tix_5prime)
        axes[ax_y, 0].set_xlabel('')

        title = file.rsplit('.',1)[0].replace('_metrics','')
        axes[ax_y, 0].set_title(title)

        # Plot 3prime figure
        df_stop.plot.bar(
            x = 'distance',
            y = 'reads',
            ax = axes[ax_y, 1],
            width = 0.9)
        axes[ax_y, 1].set_xticklabels(tix_3prime)
        axes[ax_y, 1].set_xlabel('')

        # Next file/plot line counter
        ax_y += 1

    # Save catenated figure
    plot_title = str(experiment) + '_' + str(plot_type) # Make figure title to save as from experiment name and plot type
    fig.savefig(
        str(plot_output) + plot_title + '_summary.pdf',
        dpi = dpi,
        bbox_inches = 'tight')

""""""
def compile_complexity_metrics(
        path,
        file_list,
        column_x,
        column_y,
        plot_type,
        experiment,
        plot_output,
        dpi=600):

    # Keep axes happy to avoid 'IndexError: too many indices for array' error
    # Auto formats figure summary size based on number of plots
    if (len(file_list) / 2) < 2:
        plot_rows = 2
        fig_size = (15, 16)
    else:
        plot_rows = ceil(len(file_list) / 2)
        fig_size = (15, (8 * (int(len(file_list) / 2))))

    # Set up figure space
    fig, axes = plt.subplots(
        nrows = plot_rows,
        ncols = 2,
        figsize = fig_size,
        subplot_kw = {'facecolor':'none'})
    plt.subplots_adjust(
        bottom = 0.3)

    # Initialize file and axis counters for formatting summary figure
    file_number = 0
    ax_y = 0

    for file in file_list:
        x = 0
        df = pd.read_csv(
            str(path + file),
            sep = '\t') # Initialize dataframe for relevant data
        df = df[[column_x, column_y]]
        df = df.dropna(axis = 0) # Remove rows where pertinent information is missing
        df[column_x] = np.log10(df[column_x] + 1)
        df[column_y] = df[column_y] * 100

        # Prepare subplots
        if (file_number % 2) == 0:
            ax_x = 0
        else:
            ax_x = 1

        if file_number != 0:
            if (file_number % 2) == 0:
                ax_y += 1

        # Calculate the point density (adapted from: Joe Kington, https://stackoverflow.com/a/20107592/9571488)
        x = df[str(column_x)].values
        y = df[str(column_y)].values

        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)

        idx = z.argsort() # Sort points by density
        x, y, z = x[idx], y[idx], z[idx]

        axes[ax_y, ax_x].axhline(
            -2,
            xmin = 0.048,
            ls = '-',
            color = 'black')
        axes[ax_y, ax_x].axvline(
            0,
            ymin = 0.060,
            ls = '-',
            color = 'black')
        axes[ax_y, ax_x].scatter(
            x,
            y,
            c = z,
            s = 3,
            alpha = 0.7)

        axes[ax_y, ax_x].set_xlabel('log$_1$$_0$expression (reads/kb)')
        axes[ax_y, ax_x].set_ylabel('% duplicate reads')
        axes[ax_y, ax_x].set_title(str(file.rsplit('.',1)[0]))
        axes[ax_y, ax_x].grid(False)

        file_number += 1
        del df

    # Save catenated figure
    plot_title = str(experiment) + '_' + str(plot_type) # Make figure title to save as from experiment name and plot type
    fig.savefig(
        str(plot_output) + plot_title + '_summary.pdf',
        dpi = dpi,
        bbox_inches = 'tight')

def compile_coverage(
    path,
    file_list,
    gene_name,
    feature_regions,
    sample_names,
    plot_type,
    plot_output,
    plot_color='red',
    dpi=600):

    # Get feature regions data
    exon_list = [1]
    sum = 0
    for length in feature_regions[feature_regions['feature'] == 'exon']['length'].tolist():
        sum = sum + length
        exon_list.append(sum)

    if 'five_prime_utr' in feature_regions['feature'].tolist():
        cds_start = feature_regions[feature_regions['feature'] == 'five_prime_utr']['length'].sum()
    else:
        cds_start = exon_list[0]

    if 'three_prime_utr' in feature_regions['feature'].tolist():
        cds_stop = cds_start + feature_regions[feature_regions['feature'] == 'CDS']['length'].sum()
    else:
        cds_stop = exon_list[-1] # Mark CDS until end of transcript

    # Keep axes happy to avoid 'IndexError: too many indices for array' error
    # Auto formats figure summary size based on number of plots

    # Set up figure space
    row_number = len(file_list)
    if row_number < 2:
        row_number = 2

    fig, axes = plt.subplots(
        nrows = row_number,
        ncols = 1,
        figsize = (16, 2 * row_number),
        sharey=True,
        sharex=True,
        subplot_kw = {'facecolor':'none'})
    plt.subplots_adjust(
        bottom = 0.0)

    # Initialize file and axis counters for formatting summary figure
    ax_y = 0

    for file in file_list:
        df = pd.read_csv(
            str(path) + str(file),
            sep = '\t')

        # Smoothen coverage data
        df['coverage'] = df[input_count].rolling(window = window, min_periods=1).mean()
        df['coverage'] = df['coverage'].fillna(0)
        df['feature'] = ''

        exon_number = 1
        for index, row in df.iterrows():
            # If first index, call it exon 1
            if index == 0:
                axes[ax_y].axvline(index, ls='-', linewidth=5, color='#585555', ymin=0, ymax=0.5)
                df.at[index, 'feature'] = 'Exon ' + str(exon_number)
                exon_number += 1

            # Label all other exons, but do not label end of last exon (but it will label at the beginning still)
            if index in exon_list[1:-1]:
                axes[ax_y].axvline(index, ls='-', linewidth=2, color='#585555', ymin=0, ymax=0.5)

                # If exon is too close to CDS label, skip
                if abs(index - cds_start) < 50 or abs(index - cds_stop) < 50:
                    pass
                else:
                    df.at[index, 'feature'] = 'Exon ' + str(exon_number)

                exon_number += 1

            # If last nt, add closing line
            if index == df.index[-1] - 1:
                axes[ax_y].axvline(index, ls='-', linewidth=5, color='#585555', ymin=0, ymax=0.5)

            # Label CDS start -- in this order, will overwrite previous exon label
            if index == cds_start:
                axes[ax_y].axvline(index, ls='-', linewidth=3, color='black', ymin=0, ymax=0.75)
                df.at[index, 'feature'] = 'CDS Start'

            # Label CDS stop -- in this order, will overwrite previous exon label
            if index == cds_stop:
                axes[ax_y].axvline(index, ls='-', linewidth=3, color='black', ymin=0, ymax=0.75)
                df.at[index, 'feature'] = 'CDS Stop'

        # Create x axis line
        axes[ax_y].axhline(0, ls='-', linewidth=3, color='black', xmin=0, xmax=1)

        # Plot coverage data and exon / CDS labels
        df.plot.bar(
            x = 'feature',
            y = 'coverage',
            ax = axes[ax_y],
            grid = None,
            width=1.0,
            fontsize=12,
            color = plot_color,
            linewidth=0,
            edgecolor=None)

        # Set sample labels
        axes[ax_y].set_xlabel('')
        if sample_names != None:
            axes[ax_y].set_ylabel(str(sample_names[ax_y]), fontsize=12)
        else:
            axes[ax_y].set_ylabel(str(file)[:-12], fontsize=12)

        # Only print exon / CDS labels on last sample for the overall plot
        if ax_y == row_number - 1:
            pass
        else:
            axes[ax_y].tick_params(labelbottom=False)

        ax_y += 1

    # Save catenated figure
    fig.suptitle(gene_name)
    fig.savefig(
        str(plot_output) + str(plot_type) + '_summary.pdf',
        dpi = dpi,
        bbox_inches = 'tight')

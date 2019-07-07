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
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import add_directory

input_count = 'raw_count'
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

    if plot_type == 'periodicity':
        tix = ['']
        for t in range(0,101,3):
            tix.append(str(t))
            tix.append('')
            tix.append('')

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
        if plot_type == 'periodicity':
            df.plot.bar(
                x = lab_x,
                y = lab_y,
                title = file[:-4],
                ax = axes[ax_y, ax_x],
                grid = None)
            axes[ax_y, ax_x].set_xticklabels(tix)
            axes[ax_y, ax_x].tick_params(which='major', length=7, width=2, direction='out')
        elif plot_type == 'read_distribution':
            df.plot.bar(
                x = lab_x,
                y = lab_y,
                title = file[:-4],
                ax = axes[ax_y, ax_x],
                grid = None,
                width = 0.8)
        else:
            df.plot.line(
                x = lab_x,
                y = lab_y,
                title = file[:-4],
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

""""""
def compile_complexity_metrics(
        path,
        file_list,
        column_x,
        column_y,
        strand,
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
        axes[ax_y, ax_x].set_title(str(file[:-4]))
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
    record_type,
    sample_names,
    strand,
    plot_type,
    experiment,
    plot_output,
    plot_color='red',
    dpi=600):

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

        df['coverage'] = df[input_count].rolling(window = window).mean()
        df['coverage'] = df['coverage'].fillna(0)

        if strand == '-':
            df = df.reindex(index=df.index[::-1])
            df = df.reset_index(drop=True)

        last = 0
        exon_count = 1
        start = 'start'
        for index, row in df.iterrows():
            if start == 'start':
                df.at[index, record_type] = str(record_type) + ' ' + str(exon_count)
                exon_count += 1
                last = row[0]
                start = 'stop'
            elif abs(row[0] - last) > 1:
                df.at[index, record_type] = str(record_type) + ' ' + str(exon_count)
                exon_count += 1
                last = row[0]
            else:
                df.at[index,record_type] = ''
                last = row[0]

        for index, row in df.iterrows():
            if str(record_type) in row[3]:
                axes[ax_y].axvline(index, ls='-', linewidth=2, color='black', ymin=0, ymax=1)

        axes[ax_y].axhline(0, ls='-', linewidth=5, color='black', xmin=0, xmax=1)
        axes[ax_y].axvline(0, ls='-', linewidth=5, color='black', ymin=0, ymax=1)

        df.plot.bar(
            x = record_type,
            y = 'coverage',
            ax = axes[ax_y],
            grid = None,
            width=1.0,
            fontsize=12,
            color = plot_color,
            linewidth=0,
            edgecolor=None)

        axes[ax_y].set_xlabel('')
        if sample_names != None:
            axes[ax_y].set_ylabel(str(sample_names[ax_y]), fontsize=12)
        else:
            axes[ax_y].set_ylabel(str(file)[:-12], fontsize=12)

        if ax_y == row_number - 1:
            pass
        else:
            axes[ax_y].tick_params(labelbottom=False)

        ax_y += 1

    # Save catenated figure
    fig.suptitle(gene_name)
    plot_title = str(experiment) + '_' + str(plot_type) # Make figure title to save as from experiment name and plot type
    fig.savefig(
        str(plot_output) + plot_title + '_summary.pdf',
        dpi = dpi,
        bbox_inches = 'tight')

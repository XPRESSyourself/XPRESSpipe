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
if str(matplotlib.get_backend()).lower() != 'agg':
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
else:
    import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import add_directory

input_count = 'coverage'
window = 20
conversion_table = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
    'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',
}
#Format legends
legend_colors = {
    'Start Codon': '#1b9e77',
    'Stop Codon': '#d95f02',
}

def compile_matrix_metrics(
        path,
        file_list,
        lab_x,
        lab_y,
        plot_type,
        experiment,
        plot_output,
        dpi=600):
    """Compile images from a list of metrics matrices"""


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

        if 'read_distribution' in plot_type:
            df = pd.read_csv(
                str(path) + str(file),
                sep = '\t') # Initialize dataframe for relevant data
            df.index = df[str(lab_x)]
            del df.index.name

        else:
            df = pd.read_csv(
                str(path) + str(file),
                sep = '\t',
                index_col = 0) # Initialize dataframe for relevant data

        df = df.dropna(
            axis = 0,
            subset = [str(lab_x), str(lab_y)]) # Remove rows where pertinent information is missing
        df = df.sort_index()

        # Fill in missing index info
        last_index = None

        for index in df.index.tolist():

            if last_index == None:
                last_index = index

            elif abs(index - last_index) != 1:

                for x in range(last_index + 1, index):

                    df.loc[x] = 0

            else:
                pass

            last_index = index

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

            if plot_type[:3].upper() == 'CDS':
                axes[ax_y, ax_x].set_xlabel("Representative CDS")

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


# Correct P-sites
def correct_psites(
        data,
        fasta,
        anno_dict,
        length_min=None,
        length_max=None,):
    """"""

    if length_min != None and length_max != None:
        data = data.loc[data['length'].between(length_min, length_max)]
    elif length_min != None:
        data = data.loc[data['length'] >= int(length_min)]
    elif length_max != None:
        data = data.loc[data['length'] <= int(length_max)]
    else:
        pass

    # Map codon sequence to P-sites
    seq_failed = []
    def map_sequence(name, position):

        # Convert nucleotide sequence (i.e. position 1)
        # to Python (array position 0)
        try:
            convert_position = position - 1

            codon_seq = fasta[name][convert_position:convert_position + 3]

            if len(codon_seq) != 3:
                pass
            else:
                return codon_seq.replace('T','U')

        except:
            try:
                seq_failed.append(name)
            except:
                seq_failed = []

    data['sequence'] = np.vectorize(map_sequence)(
        data['transcript'],
        data['psite'])

    # Correct P-site position in regard to CDS start and stop
    def correct_to_cds_start(transcript, psite):

        return (psite - anno_dict['l_utr5'][transcript] - 1)

    data['psite_corrected_5prime'] = np.vectorize(correct_to_cds_start)(
        data['transcript'],
        data['psite'])

    def correct_to_cds_stop(transcript, psite):

        return (
            correct_to_cds_start(transcript, psite)
            - (
            anno_dict['length'][transcript]
            - anno_dict['l_utr3'][transcript]
            )
        )

    data['psite_corrected_3prime'] = np.vectorize(correct_to_cds_stop)(
        data['transcript'],
        data['psite'])

    return data

def prep_codon(
        df_corrected,
        col_name='psite_corrected'):
    """"""

    df_corrected = df_corrected.loc[df_corrected[col_name] >= 0]
    df_corrected = df_corrected['sequence'].value_counts().drop('None')

    #### Need to fix this to remove those mapping to utr3

    return df_corrected[::-1]

def prep_periodicity_5prime(
        df_corrected,
        left_range=-20,
        right_range=90):
    """"""

    df_corrected = df_corrected.loc[df_corrected['sequence'] != "None"]
    df_corrected = df_corrected.loc[(
        df_corrected['psite_corrected_5prime'] >= left_range)
        & (df_corrected['psite_corrected_5prime'] <= right_range)]
    df_corrected = df_corrected['psite_corrected_5prime'].value_counts().sort_index()

    return df_corrected

def prep_periodicity_3prime(
        df_corrected,
        left_range=-40,
        right_range=20):
    """"""

    df_corrected = df_corrected.loc[df_corrected['sequence'] != "None"]
    df_corrected = df_corrected.loc[(
        df_corrected['psite_corrected_3prime'] >= left_range)
        & (df_corrected['psite_corrected_3prime'] <= right_range)]
    df_corrected = df_corrected['psite_corrected_3prime'].value_counts().sort_index()

    return df_corrected

def prep_bar_colors(
        df_codon):

        colors = []
        for value in df_codon.index.tolist(): # keys are the names of the boys
            if value == 'UAG' or value == 'UAA' or value == 'UGA':
                colors.append('#d95f02')
            elif value == 'AUG':
                colors.append('#1b9e77')
            else:
                colors.append('#c2c5cc')

        return colors

def prep_periodicity_steps(
        data,
        number_of_steps = 10):

    min_value = data.index.min()
    max_value = data.index.max()

    return np.arange(min_value, max_value+1, number_of_steps)

def set_lines(
        ax,
        data,
        prime5=True):

    min_value = data.index.min()
    max_value = data.index.max()

    ax.axvline(x=0, c='red')

    if prime5 == True:

        x_value = 0
        while x_value < max_value:
            x_value += 3
            ax.axvline(x=x_value, c='pink')

        x_value = 0
        while x_value > min_value:
            x_value -= 3
            ax.axvline(x=x_value, c='lightblue')

    else:

        x_value = 0
        while x_value < max_value:
            x_value += 3
            ax.axvline(x=x_value, c='lightblue')

        x_value = 0
        while x_value > min_value:
            x_value -= 3
            ax.axvline(x=x_value, c='pink')

    return ax

"""
"""
def compile_p_site_qc_metrics(
        path,
        file_list,
        fasta,
        anno_dict,
        plot_periodicity,
        plot_codon,
        experiment,
        periodicity_output,
        codon_output,
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
    fig_codon, axes_codon = plt.subplots(
        nrows = plot_rows,
        ncols = 1,
        figsize = fig_size,
        sharey = True,
        sharex = False,
        subplot_kw = {'facecolor':'none'})

    fig_period, axes_period = plt.subplots(
        nrows = plot_rows,
        ncols = 2,
        figsize = fig_size,
        sharey = True,
        sharex = False,
        subplot_kw = {'facecolor':'none'})

    # Initialize file and axis counters for formatting summary figure
    ax_y = 0

    for file in file_list:

        # Get data
        df = pd.read_csv(
            str(path) + str(file),
            sep = '\t',
            index_col = 0) # Initialize dataframe for relevant data

        df_corrected = correct_psites(
            df,
            fasta,
            anno_dict)

        # Plot codon
        df_codon = prep_codon(
            df_corrected,
            col_name='psite_corrected_5prime'
        )

        colors = prep_bar_colors(
            df_codon
        )
        
        axes_codon[ax_y, 0] = df_codon.plot.bar(width=0.9, color=colors)
        plot_position = 0
        for p in df_codon[ax_y, 0].patches:
            codon = df_codon.keys()[plot_position]
            ax.annotate(
                conversion_table[codon],
                (p.get_x() + 0.22, p.get_height() + 13000))
            plot_position += 1

        if ax_y == 0:
            f = lambda m,c: plt.plot(
                [],[],
                marker='o',
                color=c,
                markersize=10,
                ls="none")[0]
            handles = [f("s", list(legend_colors.values())[i]) for i in range(len(list(legend_colors.values())))]
            first_legend = plt.legend(
                handles,
                list(legend_colors.keys()),
                bbox_to_anchor=(0.02, 0.95),
                loc=2,
                prop={'size': 15},
                borderaxespad=0.)

            # Add the legend manually to the current Axes.
            ax = plt.gca().add_artist(first_legend)

        axes_codon[ax_y, 0].set_ylabel('# P-Sites', size=25)
        axes_codon[ax_y, 0].set_xlabel('Codon', size=25)

        # Plot periodicity
        df_periodicity_5prime = prep_periodicity_5prime(
            df_corrected
        )
        df_periodicity_3prime = prep_periodicity_3prime(
            df_corrected
        )

        # 5' and 3'
        axes_period[ax_y, 0] = df_periodicity_5prime.plot.line()
        axes_period[ax_y, 1] = df_periodicity_3prime.plot.line()

        steps_5prime = prep_periodicity_steps(df_periodicity_5prime)
        steps_3prime = prep_periodicity_steps(df_periodicity_3prime)

        axes_period[ax_y, 0].set(xticks=steps_5prime, xticklabels=steps_5prime)
        axes_period[ax_y, 1].set(xticks=steps_3prime, xticklabels=steps_3prime)

        axes_period[ax_y, 0] = set_lines(
            axes_period[ax_y, 0],
            df_periodicity_5prime,
            prime5=True)
        axes_period[ax_y, 1] = set_lines(
            axes_period[ax_y, 1],
            df_periodicity_3prime,
            prime5=False)

        # Next file/plot line counter
        ax_y += 1

    # Save catenated figures
    fig_codon.savefig(
        (
            str(codon_output)
            + str(experiment)
            + '_' + plot_codon
            + '_summary.pdf'),
        dpi = dpi,
        bbox_inches = 'tight')

    fig_period.savefig(
        (
            str(periodicity_output)
            + str(experiment)
            + '_' + plot_periodicity
            + '_summary.pdf'),
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

        stop_written = False
        start_written = False
        exon_number = 1
        for index, row in df.iterrows():

            # If first exon is also CDS start, call it so
            if cds_start == exon_list[0] and start_written == False:
                axes[ax_y].axvline(index, ls='-', linewidth=3, color='black', ymin=0, ymax=0.75)
                df.at[index, 'feature'] = 'CDS Start'
                exon_number += 1
                start_written = True

            # If first index, call it exon 1
            elif index == 0 and start_written == False:
                axes[ax_y].axvline(index, ls='-', linewidth=5, color='#585555', ymin=0, ymax=0.5)
                df.at[index, 'feature'] = 'Exon ' + str(exon_number)
                exon_number += 1

            else:
                pass

            # Label all other exons, but do not label end of last exon (but it will label at the beginning still)
            if index in exon_list[1:-1]:
                axes[ax_y].axvline(index, ls='-', linewidth=2, color='#585555', ymin=0, ymax=0.5)

                # If exon is too close to CDS label, skip
                if abs(index - cds_start) < 50 or abs(index - cds_stop) < 50:
                    pass
                else:
                    df.at[index, 'feature'] = 'Exon ' + str(exon_number)

                exon_number += 1

            # Label CDS start -- in this order, will overwrite previous exon label
            if index == cds_start and start_written == False:
                axes[ax_y].axvline(index, ls='-', linewidth=3, color='black', ymin=0, ymax=0.75)
                df.at[index, 'feature'] = 'CDS Start'
                start_written = True

            # Label CDS stop -- in this order, will overwrite previous exon label
            if index == cds_stop and stop_written == False:
                axes[ax_y].axvline(index, ls='-', linewidth=3, color='black', ymin=0, ymax=0.75)
                df.at[index, 'feature'] = 'CDS Stop'
                stop_written = True

        if stop_written == False:
            axes[ax_y].axvline(index, ls='-', linewidth=3, color='black', ymin=0, ymax=0.75)
            df.at[df.index[-1] - 1, 'feature'] = 'CDS Stop'
            stop_written = True

        # If last nt, add closing line
        elif index == df.index[-1] - 1:
            axes[ax_y].axvline(index, ls='-', linewidth=5, color='#585555', ymin=0, ymax=0.5)
        else:
            pass

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

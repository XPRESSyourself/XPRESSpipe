"""
XPRESSplot
A toolkit for navigating and analyzing gene expression datasets
alias: xpressplot

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

#rpm, r_fpkm, batch_normalize, convert_names, tpm
#from .xpressplot import batch_normalize, rpm, tpm, r_fpkm
#from .xpressplot import count_table
"""IMPORT DEPENDENCIES"""
import os
import sys
import csv
import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib
import matplotlib.pyplot as plt
if str(matplotlib.get_backend()).lower() != 'agg':
    plt.switch_backend('agg')
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'sans-serif']
plt.rcParams['font.family'] = 'sans-serif'

"""INITIALIZATION PARAMETERS"""
# Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpressplot_arguments = os.path.split(__file__)
gtf_type_column = 2
gtf_leftCoordinate_column = 3
gtf_rightCoordinate_column = 4
gtf_annotation_column = 8
gtf_transcript_column = 9
gtf_gene_id_column = 10
gtf_gene_name_column = 11
id_search = r'transcript_id \"(.*?)\"; '
gene_id_search = r'gene_id \"(.*?)\"; '
gene_name_search = r'gene_name \"(.*?)\"; '


"""Collate counts files"""
def count_table(
    file_list,
    gene_column=0, sample_column=1,
    sep='\t', drop_rows=5):

    # Read in first count file to get gene names
    df = pd.read_csv(str(file_list[0]), sep=sep, comment='#', header=None)
    pos_starter = [gene_column,sample_column]
    colname = df.columns[pos_starter]
    df = df[colname]

    # For the rest of the files in the file list, add the counts for that sample only
    for f in file_list[1:]:
        df_pull = pd.read_csv(str(f), sep=sep, comment='#', header=None)
        df = pd.concat([df, df_pull[df_pull.columns[sample_column]]], axis=1)
        df_pull = None

    # Final formatting clean up of table
    df_counts = df.copy()
    df = None
    df_counts = df_counts.set_index(0)
    df_counts.index.name = None

    # Remove path and file suffix from each file's name before adding as column names to table
    c = 0
    for x in file_list:
        file_list[c] = x[(x.rfind('/')+1):(x.find('.'))]
        c += 1
    df_counts.columns = file_list
    df_counts = df_counts[:-drop_rows]

    return df_counts

"""Convert row names (genes) of dataframe using GTF as reference for new name"""
def convert_names(
    data,
    gtf,
    orig_name_label='gene_id',
    orig_name_location=0,
    new_name_label='gene_name',
    new_name_location=2,
    refill=None,
    add_space=True,
    sep='\t'):

    if add_space == True:
        orig_name_label = orig_name_label + ' \"'
        new_name_label = new_name_label + ' \"'

    # Import reference GTF
    gtf = pd.read_csv(str(gtf),sep=sep,comment='#', low_memory=False, header=None)

    # Parse out old and new names from GTF
    gtf_genes = gtf.loc[gtf[2] == 'gene']
    gtf_genes['original'] = gtf[8].str.split(';').str[orig_name_location]
    gtf_genes['new'] = gtf[8].str.split(';').str[new_name_location]
    gtf_genes['original'] = gtf_genes['original'].map(lambda x: x.lstrip(str(orig_name_label)).rstrip('\"').rstrip(' '))
    gtf_genes['new'] = gtf_genes['new'].map(lambda x: x.lstrip(str(new_name_label)).rstrip('\"').rstrip(' '))
    gtf_genes = gtf_genes[['original','new']].copy()

    # Create dictionary
    if refill != None:
        gene_dict = {}
        for index, row in gtf_genes.iterrows():
            if row[1] == str(refill):
                gene_dict[row[0]] = row[0]
            else:
                gene_dict[row[0]] = row[1]
    else:
        gene_dict = pd.Series(gtf_genes['new'].values,index=gtf_genes['original']).to_dict()

    # Replace old gene names/ids with new
    data_names = data.copy()
    data_names['new'] = data_names.index.to_series().map(gene_dict).fillna(data_names.index.to_series())
    data_names = data_names.set_index('new')
    data_names.index.name = None

    return data_names

"""Create gene dictionary"""
# Provide Longest GTF to get longest canonical record
def gene_length_dictionary(
    gtf,
    feature_type='exon', # or CDS
    identifier='gene_name', # or gene_id or transcript_id
    sep='\t'):

    if not str(gtf).endswith('.gtf'):
        print(str(gtf) + ' does not appear to be a GTF file')

    if not str(feature_type).lower() in ['cds', 'exon']:
        print('Must provide CDS or exon as feature_type')

    # Process gtf data for gene_name and gene_length
    gtf = pd.read_csv(
        str(gtf),
        sep = '\t',
        header = None,
        comment = '#',
        low_memory = False)

    # Get relevant metadata
    gtf = gtf[(gtf[gtf_type_column] == feature_type)]
    gtf[gtf_transcript_column] = gtf[gtf_annotation_column].str.extract(id_search)
    gtf[gtf_gene_id_column] = gtf[gtf_annotation_column].str.extract(gene_id_search)
    gtf[gtf_gene_name_column] = gtf[gtf_annotation_column].str.extract(gene_name_search)
    gtf['length'] = abs(gtf[gtf_rightCoordinate_column] - gtf[gtf_leftCoordinate_column]) + 1

    # Return a transcript id mapping dictionary for isoform normalization
    if identifier == 'transcript_id':
        gtf = gtf[[gtf_transcript_column, 'length']]
        gtf.columns = ['transcript', 'length']
        gtf = gtf.dropna()
        length_dict = pd.Series(gtf['length'].values,index=gtf['transcript']).to_dict()
        return length_dict

    else:
        if identifier == 'gene_id':
            col = gtf_gene_id_column
        else:
            col = gtf_gene_name_column

        gtf = gtf[[col, gtf_transcript_column, 'length']]
        gtf.columns = ['gene', 'transcript', 'length']

        # Make transcript gene mapping dictionary
        ref = gtf[['transcript', 'gene']]
        reference = pd.Series(ref['gene'].values,index=ref['transcript']).to_dict()

        # Group by transcript ID and get feature length sum
        records = gtf[['transcript', 'length']]
        records = records.sort_values('transcript')
        records = records.groupby('transcript').sum()
        records = records.reset_index()

        # Map exon space to each bam record based on its transcript ID
        records['gene'] = records['transcript'].map(reference)

        # Get longest transcript based on max exon or CDS space size along (not Ensembl canonical)
        records = records.sort_values('gene')
        records = records.groupby('gene').max()
        records = records.reset_index()

        length_dict = pd.Series(records['length'].values,index=records['gene']).to_dict()
        return length_dict


"""Perform gene kilobase normalization"""
def rpk(
    data,
    length_dict):

    data_c = data.copy()

    data_rpk = data_c.div(data_c.index.map(length_dict) / 1000, axis = 0)
    data_rpk['index'] = data_rpk.index
    data_rpk = data_rpk.drop_duplicates(subset='index', keep='first')
    #data_rpk = data_rpk.reindex(list(length_dict.keys()), axis=0)
    data_rpk = data_rpk.drop('index', axis=1)
    data_rpk = data_rpk.fillna(0)

    return data_rpk

"""Perform reads per million sample normalization on RNAseq data"""
def rpm(
    data):

    data_rpm = data / \
        (data.sum() / 1e6)
    data_rpm = data_rpm.fillna(0)

    return data_rpm

"""Perform reads/fragments per kilobase million sample normalization on RNAseq data"""
def r_fpkm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    length_dict = gene_length_dictionary(
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

    data_rpm = rpm(data)
    data_rpkm = rpk(
        data_rpm,
        length_dict)

    return data_rpkm

def rpkm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    r_fpkm(
        data,
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

def fpkm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    r_fpkm(
        data,
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

"""Perform transcripts per million normalization on RNAseq data"""
def tpm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    length_dict = gene_length_dictionary(
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

    data_rpk = rpk(
        data,
        length_dict)
    data_tpm = rpm(data_rpk)

    return data_tpm

"""Normalize out batch effects from RNAseq data"""
def batch_normalize(
    input_file,
    batch_file):

    # Get output file name
    if input_file.endswith('.txt') or input_file.endswith('.tsv'):
        output_file = str(input_file[:-4]) + '_batched.tsv'
    else:
        raise Exception('Unrecognized input_file delimiter type. Files must be tab-delimited')

    # Run sva combat in R
    os.system('rscript' \
        + ' ' + str(__path__) + '/Rbatch_normalize.r' \
        + ' ' + str(input_file) \
        + ' ' + str(batch_file) \
        + ' ' + str(output_file))
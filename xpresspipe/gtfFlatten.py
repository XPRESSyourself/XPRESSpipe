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

"""IMPORT DEPENDENCIES"""
import pandas as pd

"""IMPORT INTERNAL DEPENDENCIES"""
from .gtfModify import get_chunks, run_chunks, longest_transcripts, protein_gtf

"""Make a flat reference from parsed GTF file"""
def make_flatten(gtf):

    records = []

    for index, row in gtf.iterrows():

        if row[2] == 'transcript':

            gene = gtf.at[index, 8][(gtf.at[index, 8].find('gene_id \"') + 9):].split('\";')[0]
            strand = row[6]
            coordinates = []

            # Scan for exons for the given
            n = 0
            item = gene
            while item == gene:
                n += 1

                if index + n > len(gtf.index) - 1 or gtf.at[index + n, 2] == 'transcript':
                    break
                else:
                    item = gtf.at[index + n, 8][(gtf.at[index + n, 8].find('gene_id \"') + 9):].split('\";')[0]

                    if gtf.at[index + n, 2] == 'exon': # Append coordinate paires for each exon of the transcript
                        coordinates.append(gtf.at[index + n, 3])
                        coordinates.append(gtf.at[index + n, 4])

            # Get start and end positions for transcript/gene
            # Assumes start and stop will be start and end of a protein coding transcript
            start = min(coordinates)
            end = max(coordinates)

            # Push information to record
            records.append([gene, strand, start, end, coordinates])

    # Push flattened reference into pandas dataframe
    headers = ['gene', 'strand', 'start', 'end', 'coordinates']
    reference = pd.DataFrame(records, columns=headers)

    return reference

"""Read in coordinate reference from GTF"""
def flatten_reference(
    gtf,
    threads=None):

    # Import GTF reference file
    if isinstance(gtf, pd.DataFrame) and len(gtf.columns) == 9:
        pass
    elif str(gtf).endswith('.gtf'):
        gtf = pd.read_csv(
                str(gtf),
                sep = '\t',
                header = None,
                comment = '#',
                low_memory = False)
    else:
        raise Exception('Error: A GTF-formatted file or dataframe was not provided')

    # Get chunks
    chunks = get_chunks(
                    gtf,
                    threads = threads)

    # Get longest transcripts
    chunks = run_chunks(
                longest_transcripts,
                chunks,
                target_message = 'longest transcripts')

    # Get only protein coding annotated records
    chunks = run_chunks(
                protein_gtf,
                chunks,
                target_message = 'protein coding genes')

    # Rejoin chunks into single GTF and flatten
    if len(chunks) > 0:
        gtf = pd.concat(chunks)
        gtf = gtf.reset_index(drop=True)
        gtf = make_flatten(gtf)

        return gtf
    else:
        return

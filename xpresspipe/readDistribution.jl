# """
# XPRESSpipe
# An alignment and analysis pipeline for RNAseq data
# alias: xpresspipe
#
# Copyright (C) 2019  Jordan A. Berg
# jordan <dot> berg <at> biochem <dot> utah <dot> edu
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
# """

# Import dependencies
import Pkg

# Install CSV if not found and load
try
    using CSV

catch e
    Pkg.add("CSV")

    using CSV

end

# Install DataFrames if not found and load
try
    using DataFrames

catch e
    Pkg.add("DataFrames")

    using DataFrames

end

# Install DataStructures if not found and load
try
    using DataStructures

catch e
    Pkg.add("DataStructures")

    using DataStructures

end

# Parse FASTQ line by line and accumulate stats
function readFastq(filename::String)
    counter = 0
    seq_dict = SortedDict()

    for seq in eachline(filename)
        if (counter - 1) % 4 == 0
            seq_len = length(seq)

            if haskey(seq_dict, seq_len)
                seq_dict[seq_len] += 1
            else
                seq_dict[seq_len] = 1
            end


        end
        counter += 1

    end

    return seq_dict

end

# Parse a single FASTQ file from PE pair
function parseFastq(filename::String)

    counter = 0
    seq_list = []

    for seq in eachline(filename)
        if (counter - 1) % 4 == 0
            seq_len = length(seq)

            append!( seq_list, seq_len )

        end
        counter += 1

    end

    return seq_list

end

# Parse PE FASTQ line by line and accumulate stats
function readFastqPE(filename1::String, filename2::String)

    seq_list1 = parseFastq(filename1)
    seq_list2 = parseFastq(filename2)

    try
        seq_comb = map((seq_list1, seq_list2) -> seq_list1 + seq_list2, seq_list1, seq_list2)

        seq_dict = SortedDict()

        for seq_len in seq_comb
            if haskey(seq_dict, seq_len)
                seq_dict[seq_len] += 1
            else
                seq_dict[seq_len] = 1
            end
        end

        return seq_dict

    catch e
        println("It appears your paired-end reads may be of different lengths")

    end

end

# Save metrics to output file
function saveMetrics(seq_dict::SortedDict, output::String)

    metrics = DataFrame(A = collect(keys(seq_dict)), B = collect(values(seq_dict)))

    colnames = ["read size (bp)", "count"]
    rename!(metrics, Symbol.(colnames))

    metrics |> CSV.write(output, delim='\t')

end

# """
# Run main
# file1::String
# file2::String (must be none if SE)
# output::String
# """
function main(args)
    fastq1 = args[1]
    fastq2 = args[2]
    output = args[3]

    if fastq2 == "None"
        seq_dict = readFastq(fastq1)
    else
        seq_dict = readFastqPE(fastq1, fastq2)
    end

    saveMetrics(seq_dict, output)

end

main(ARGS)

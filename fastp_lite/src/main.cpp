//
//  main.cpp
//  fastp_lite
//
//  Created by Jordan Berg on 3/12/20.
//  Copyright © 2020 Jordan Berg. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <ctime>
#include <string.h>
#include "fastqreader.h"
#include "cmdline.h"

#define FASTP_LITE_VER "0.0.1"

using namespace std;

// Timer functions
clock_t START_TIMER;

clock_t start()
{
    return START_TIMER = clock();
}

void stop(clock_t start = START_TIMER)
{
    cout
        << "Elapsed time:             "
        << (clock() - start) / (double)CLOCKS_PER_SEC
        << "s\n"
        << endl;
}

int main(int argc, char* argv[]) {

    // Parse user arguments
    if(argc == 1) {
        cerr
            << "fastp_lite: an ultra-fast 3' UMI preprocessor"
            << endl
            << "version "
            << FASTP_LITE_VER
            << endl;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr
            << "fastp_lite "
            << FASTP_LITE_VER
            << endl;
        return 0;
    }

    cmdline::parser cmd;
    cmd.add<string>("input", 'i', "Input file name", false, ""); //file
    cmd.add<string>("output", 'o', "Output file name (cannot be the same name as the input file)", false, ""); //out_file
    cmd.add<int>("umi-length", 'l', "Length of UMI, default is 0", false, 0); //_3prime_umi_len
    cmd.add<int>("spacer-length", 's', "Length of UMI spacer, default is 0", false, 0); //_3prime_umi_spacer
    cmd.add<int>("min-length", 'm', "Minimum length of reads after UMI processing, default is 10", false, 10); //_3prime_umi_spacer

    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }

    string file;
    file = cmd.get<string>("input");

    string out_file;
    out_file = cmd.get<string>("output");

    int _3prime_umi_len;
    _3prime_umi_len = cmd.get<int>("umi-length");

    int _3prime_umi_spacer;
    _3prime_umi_spacer = cmd.get<int>("spacer-length");

    int _min_length;
    _min_length = cmd.get<int>("min-length");

    // Check the inputs and exit if one doesn't work
    if (file == out_file) {
        cout
            << "Error: Input and output file names can not be identical"
            << endl;
        return 1;
      }
    if (file == "") {
        cout
            << "Error: Input file name cannot be blank"
            << endl;
        return 1;
      }
    if (out_file == "") {
        cout
            << "Error: Output file name cannot be blank"
            << endl;
        return 1;
      }
    if (_3prime_umi_len < 0) {
        cout
            << "Error: 3' UMI length cannot be less than 0"
            << endl;
        return 1;
      }
    if (_3prime_umi_spacer < 0) {
        cout
            << "Error: 3' UMI spacer length cannot be less than 0"
            << endl;
        return 1;
      }
    if (_min_length < 0) {
        cout
            << "Error: Minimum read length cannot be less than 0"
            << endl;
        return 1;
      }

    // Print user inputs
    cout
        << "\nfastp_lite v"
        << FASTP_LITE_VER << endl << endl;
    cout
        << "Provided input file name:                             "
        << file << endl;
    cout
        << "Provided output file name:                            "
        << out_file << endl;
    cout
        << "Selected 3' UMI length:                               "
        << _3prime_umi_len << endl;
    cout
        << "Selected 3' UMI spacer length:                        "
        << _3prime_umi_spacer << endl;
    cout
        << "Minimum allowable read length after UMI processing:   "
        << _min_length << endl << endl;
    cout
        << "Processing 3' internal UMIs from reads...\n" << endl;

    // Process file
    start(); // start elapsed time
    fstream output(out_file, std::fstream::in | std::fstream::out | std::fstream::trunc); // prepare output file
    FastqReader reader1(file); // initialize input FASTQ file
    Read* r1 = NULL;

    // Initialize counters
    int reads_keep = 0;
    int reads_toss = 0;
    int read_lengths_sum = 0;

    // Process reads from FASTQ file for 3' UMIs and spacers
    while(true){
        r1 = reader1.read();
        if (r1 == NULL) {
            break;
          }
        else {
            r1->toString();
            string umi;
            umi = r1->mSeq.mStr.substr(r1->length() - _3prime_umi_len, r1->length());
            r1->addUmiTag(umi);
            r1->trimBack(_3prime_umi_len);
            r1->trimBack(_3prime_umi_spacer);

            int read_length;
            read_length = r1->length();
            if (read_length >= _min_length) {
                output << r1->toString();
                reads_keep += 1;
                read_lengths_sum += read_length;
              }
            else {
                reads_toss += 1;
              }
          }

        delete r1;
    }

    output.close(); // close output file

    // Print results
    cout
        << "Processing complete.      "
        << endl;
    cout
        << "Total reads processed:    "
        << reads_keep + reads_toss
        << endl;
    cout
        << "Total reads kept:         "
        << reads_keep
        << endl;
    cout
        << "Total reads removed:      "
        << reads_toss
        << endl;
    if (reads_keep != 0) {
        int avg_read = static_cast<int>(read_lengths_sum / reads_keep);
        cout
            << "Average read length:      "
            << avg_read
            << endl;
      }
    stop(); // stop and print elapsed time
    cout.flush();
    return 0;
}

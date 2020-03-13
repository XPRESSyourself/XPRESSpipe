//
//  main.cpp
//  fastp_lite
//
//  Created by Jordan Berg on 3/12/20.
//  Copyright © 2020 Jordan Berg. All rights reserved.
//
/*  Adapted from fastp for 3' UMI trimming and recording
    MIT License

    Copyright (c) 2017 OpenGene - Open Source Genetics Toolbox

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
 */

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
        << "Elapsed time: "
        << (clock() - start) / (double)CLOCKS_PER_SEC << "s"
        << endl;
}

int main(int argc, char* argv[]) {

    // Parse user arguments
    if(argc == 1) {
        cerr << "fastp_lite: an ultra-fast 3' UMI preprocessor" << endl << "version " << FASTP_LITE_VER << endl;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "fastp_lite " << FASTP_LITE_VER << endl;
        return 0;
    }
    
    cmdline::parser cmd;
    cmd.add<string>("input", 'i', "Input file name", false, ""); //file
    cmd.add<string>("output", 'o', "Output file name", false, ""); //out_file
    cmd.add<int>("umi-length", 'l', "Length of UMI, default is 0", false, 0); //_3prime_umi_len
    cmd.add<int>("spacer-length", 's', "Length of UMI spacer, default is 0", false, 0); //_3prime_umi_spacer
    
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
    
    // Print user inputs
    cout
        << "\nfastp_lite v"
        << FASTP_LITE_VER;
    cout
        << "\n";
    cout
        << "\nProvided input file name:         "
        << file;
    cout
        << "\nProvided output file name:        "
        << out_file;
    cout
        << "\nSelected 3' UMI length:           "
        << _3prime_umi_len;
    cout
        << "\nSelected 3' UMI spacer length:    "
        << _3prime_umi_spacer;
    cout
        << "\n\n";
    cout
        << "Processing UMIs from SE reads...\n\n";

    // Process file
    start(); // start elapsed time
    fstream output(out_file, std::fstream::in | std::fstream::out | std::fstream::trunc); // prepare output file
    FastqReader reader1(file); // initialize input FASTQ file
    Read* r1 = NULL;
    
    // Process reads from FASTQ file for 3' UMIs and spacers
    while(true){
        r1 = reader1.read();
        if(r1 == NULL)
            break;
        else
            r1->toString();
            string umi;
            umi = r1->mSeq.mStr.substr(r1->length() - _3prime_umi_len, r1->length());
            r1->addUmiTag(umi);
            r1->trimBack(_3prime_umi_len);
            r1->trimBack(_3prime_umi_spacer);
            output << r1->toString();
            
        delete r1;
    }
    
    output.close(); // close output file
    stop(); // stop and print elapsed time
    return 0;
}
    
    
    

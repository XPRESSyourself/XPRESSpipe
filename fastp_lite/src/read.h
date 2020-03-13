//
//  read.h
//  fastp_lite
//
//  Created by Jordan Berg on 3/12/20.
//  Copyright © 2020 Jordan Berg. All rights reserved.
//
/* Adapted from fastp for 3' UMI trimming and recording
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

#ifndef READ_H
#define READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "sequence.h"
#include <vector>

using namespace std;

class Read{
public:
    Read(string name, string seq, string strand, string quality, bool phred64=false);
    Read(string name, Sequence seq, string strand, string quality, bool phred64=false);
    Read(string name, string seq, string strand);
    Read(string name, Sequence seq, string strand);
    Read(Read &r);
    void print();
    void printFile(ofstream& file);
    Read* reverseComplement();
    string firstIndex();
    string lastIndex();
    // default is Q20
    int lowQualCount(int qual=20);
    int length();
    string toString();
    string toStringWithTag(string tag);
    void resize(int len);
    void convertPhred64To33();
    void trimBack(int len);
    void addUmiTag(string umi);

public:
    static bool test();

private:


public:
    string mName;
    Sequence mSeq;
    string mStrand;
    string mQuality;
    bool mHasQuality;
};

class ReadPair{
public:
    ReadPair(Read* left, Read* right);
    ~ReadPair();

    // merge a pair, without consideration of seq error caused false INDEL
    Read* fastMerge();
public:
    Read* mLeft;
    Read* mRight;

public:
    static bool test();
};

#endif

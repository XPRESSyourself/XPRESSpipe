/* From fastp for 3' UMI trimming and recording
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

#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#include <iostream>
#include <fstream>

using namespace std;

class FastqReader{
public:
    FastqReader(string filename, bool hasQuality = true, bool phred64=false);
    ~FastqReader();

    void getBytes(size_t& bytesRead, size_t& bytesTotal);

    //this function is not thread-safe
    //do not call read() of a same FastqReader object from different threads concurrently
    Read* read();
    bool eof();
    bool hasNoLineBreakAtEnd();

public:
    static bool isFastq(string filename);
    static bool test();

private:
    void init();
    void close();
    string getLine();
    void clearLineBreaks(char* line);
    void readToBuf();

private:
    string mFilename;
    FILE* mFile;
    bool mHasQuality;
    bool mPhred64;
    char* mBuf;
    int mBufDataLen;
    int mBufUsedLen;
    bool mStdinMode;
    bool mHasNoLineBreakAtEnd;

};

class FastqReaderPair{
public:
    FastqReaderPair(FastqReader* left, FastqReader* right);
    FastqReaderPair(string leftName, string rightName, bool hasQuality = true, bool phred64 = false, bool interleaved = false);
    ~FastqReaderPair();
    ReadPair* read();
public:
    FastqReader* mLeft;
    FastqReader* mRight;
    bool mInterleaved;
};

#endif


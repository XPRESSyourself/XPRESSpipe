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

#include "fastqreader.h"
#include "util.h"
#include <string.h>

using namespace std;

#define FQ_BUF_SIZE (1<<20)

FastqReader::FastqReader(string filename, bool hasQuality, bool phred64){
    mFilename = filename;
    mFile = NULL;
    mStdinMode = false;
    mPhred64 = phred64;
    mHasQuality = hasQuality;
    mBuf = new char[FQ_BUF_SIZE];
    mBufDataLen = 0;
    mBufUsedLen = 0;
    mHasNoLineBreakAtEnd = false;
    init();
}

FastqReader::~FastqReader(){
    close();
    delete mBuf;
}

bool FastqReader::hasNoLineBreakAtEnd() {
    return mHasNoLineBreakAtEnd;
}

void FastqReader::readToBuf() {
    
    mBufDataLen = fread(mBuf, 1, FQ_BUF_SIZE, mFile);

    mBufUsedLen = 0;

    if(mBufDataLen < FQ_BUF_SIZE) {
        if(mBuf[mBufDataLen-1] != '\n')
            mHasNoLineBreakAtEnd = true;
    }
}

void FastqReader::init(){

    if(mFilename == "/dev/stdin") {
        mFile = stdin;
    }
    else
        mFile = fopen(mFilename.c_str(), "rb");
    if(mFile == NULL) {
        error_exit("Failed to open file: " + mFilename);
    }
    
    readToBuf();
}

void FastqReader::getBytes(size_t& bytesRead, size_t& bytesTotal) {

    bytesRead = ftell(mFile);//mFile.tellg();

    // use another ifstream to not affect current reader
    ifstream is(mFilename);
    is.seekg (0, is.end);
    bytesTotal = is.tellg();
}

void FastqReader::clearLineBreaks(char* line) {

    // trim \n, \r or \r\n in the tail
    int readed = strlen(line);
    if(readed >=2 ){
        if(line[readed-1] == '\n' || line[readed-1] == '\r'){
            line[readed-1] = '\0';
            if(line[readed-2] == '\r')
                line[readed-2] = '\0';
        }
    }
}

string FastqReader::getLine(){
    static int c=0;
    c++;
    int copied = 0;

    int start = mBufUsedLen;
    int end = start;

    while(end < mBufDataLen) {
        if(mBuf[end] != '\r' && mBuf[end] != '\n')
            end++;
        else
            break;
    }

    // this line well contained in this buf, or this is the last buf
    if(end < mBufDataLen || mBufDataLen < FQ_BUF_SIZE) {
        int len = end - start;
        string line(mBuf+start, len);

        // skip \n or \r
        end++;
        // handle \r\n
        if(end < mBufDataLen-1 && mBuf[end-1]=='\r' && mBuf[end] == '\n')
            end++;

        mBufUsedLen = end;

        return line;
    }

    // this line is not contained in this buf, we need to read new buf
    string str(mBuf+start, mBufDataLen - start);

    while(true) {
        readToBuf();
        start = 0;
        end = 0;
        while(end < mBufDataLen) {
            if(mBuf[end] != '\r' && mBuf[end] != '\n')
                end++;
            else
                break;
        }
        // this line well contained in this buf, we need to read new buf
        if(end < mBufDataLen || mBufDataLen < FQ_BUF_SIZE) {
            int len = end - start;
            str.append(mBuf+start, len);

            // skip \n or \r
            end++;
            // handle \r\n
            if(end < mBufDataLen-1 && mBuf[end] == '\n')
                end++;

            mBufUsedLen = end;
            return str;
        }
        // even this new buf is not enough, although impossible
        str.append(mBuf+start, mBufDataLen);
    }

    return string();
}

bool FastqReader::eof() {

    return feof(mFile);//mFile.eof();
}

Read* FastqReader::read(){

    if(mBufUsedLen >= mBufDataLen && eof()) {
        return NULL;
    }

    string name = getLine();
    // name should start with @
    while((name.empty() && !(mBufUsedLen >= mBufDataLen && eof())) || (!name.empty() && name[0]!='@')){
        name = getLine();
    }

    if(name.empty())
        return NULL;

    string sequence = getLine();
    string strand = getLine();

    // WAR for FQ with no quality
    if (!mHasQuality){
        string quality = string(sequence.length(), 'K');
        return new Read(name, sequence, strand, quality, mPhred64);
    }
    else {
        string quality = getLine();
        if(quality.length() != sequence.length()) {
            cerr << "ERROR: sequence and quality have different length:" << endl;
            cerr << name << endl;
            cerr << sequence << endl;
            cerr << strand << endl;
            cerr << quality << endl;
            return NULL;
        }
        return new Read(name, sequence, strand, quality, mPhred64);
    }

    return NULL;
}

void FastqReader::close(){

    if (mFile){
        fclose(mFile);//mFile.close();
        mFile = NULL;
    }
}

bool FastqReader::isFastq(string filename) {
    if (ends_with(filename, ".fastq"))
        return true;
    else if (ends_with(filename, ".fq"))
        return true;
    else
        return false;
}

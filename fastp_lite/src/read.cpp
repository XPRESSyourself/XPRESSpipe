//
//  read.cpp
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

#include "read.h"
#include <sstream>
#include "util.h"

Read::Read(string name, string seq, string strand, string quality, bool phred64){
    mName = name;
    mSeq = Sequence(seq);
    mStrand = strand;
    mQuality = quality;
    mHasQuality = true;
    if(phred64)
        convertPhred64To33();
}

Read::Read(string name, string seq, string strand){
    mName = name;
    mSeq = Sequence(seq);
    mStrand = strand;
    mHasQuality = false;
}

Read::Read(string name, Sequence seq, string strand, string quality, bool phred64){
    mName = name;
    mSeq = seq;
    mStrand = strand;
    mQuality = quality;
    mHasQuality = true;
    if(phred64)
        convertPhred64To33();
}

Read::Read(string name, Sequence seq, string strand){
    mName = name;
    mSeq = seq;
    mStrand = strand;
    mHasQuality = false;
}

void Read::convertPhred64To33(){
    for(int i=0; i<mQuality.length(); i++) {
        mQuality[i] = max(33, mQuality[i] - (64-33));
    }
}

Read::Read(Read &r) {
    mName = r.mName;
    mSeq = r.mSeq;
    mStrand = r.mStrand;
    mQuality = r.mQuality;
    mHasQuality = r.mHasQuality;
}

void Read::print(){
    std::cerr << mName << endl;
    std::cerr << mSeq.mStr << endl;
    std::cerr << mStrand << endl;
    if(mHasQuality)
        std::cerr << mQuality << endl;
}

void Read::printFile(ofstream& file){
    file << mName << endl;
    file << mSeq.mStr << endl;
    file << mStrand << endl;
    if(mHasQuality)
        file << mQuality << endl;
}

Read* Read::reverseComplement(){
    Sequence seq = ~mSeq;
    string qual;
    qual.assign(mQuality.rbegin(), mQuality.rend());
    string strand = (mStrand=="+") ? "-" : "+";
    return new Read(mName, seq, strand, qual);
}

void Read::resize(int len) {
    if(len > length() || len<0)
        return ;
    mSeq.mStr.resize(len);
    mQuality.resize(len);
}
   
void Read::trimBack(int len){
    len = min(length()-1, len);
    mSeq.mStr = mSeq.mStr.substr(0, mSeq.mStr.length() - len);
    mQuality = mQuality.substr(0, mQuality.length() - len);
}

string Read::lastIndex(){
    int len = mName.length();
    if(len<5)
        return "";
    for(int i=len-3;i>=0;i--){
        if(mName[i]==':' || mName[i]=='+'){
            return mName.substr(i+1, len-i);
        }
    }
    return "";
}

string Read::firstIndex(){
    int len = mName.length();
    int end = len;
    if(len<5)
        return "";
    for(int i=len-3;i>=0;i--){
        if(mName[i]=='+')
            end = i-1;
        if(mName[i]==':'){
            return mName.substr(i+1, end-i);
        }
    }
    return "";
}

int Read::lowQualCount(int qual){
    int count = 0;
    for(int q=0;q<mQuality.size();q++){
        if(mQuality[q] < qual + 33)
            count++;
    }
    return count;
}

int Read::length(){
    return mSeq.length();
}

string Read::toString() {
    return mName + "\n" + mSeq.mStr + "\n" + mStrand + "\n" + mQuality + "\n";
}

string Read::toStringWithTag(string tag) {
    return mName + " " + tag + "\n" + mSeq.mStr + "\n" + mStrand + "\n" + mQuality + "\n";
}

void Read::addUmiTag(string umi){
    
    string tag;
    tag = ":" + umi;

    int spacePos = -1;
    for(int i=0; i<mName.length(); i++) {
        if(mName[i] == ' ') {
            spacePos = i;
            break;
        }
    }
    if(spacePos == -1) {
        mName = mName + tag;
    } else {
        mName = mName.substr(0, spacePos) + tag + mName.substr(spacePos, mName.length() - spacePos);
    }

}

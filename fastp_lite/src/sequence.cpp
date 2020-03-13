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

#include "sequence.h"

Sequence::Sequence(){
}

Sequence::Sequence(string seq){
    mStr = seq;
}

void Sequence::print(){
    std::cerr << mStr;
}

int Sequence::length(){
    return mStr.length();
}

Sequence Sequence::reverseComplement(){
    string str(mStr.length(), 0);
    for(int c=0;c<mStr.length();c++){
        char base = mStr[c];
        switch(base){
            case 'A':
            case 'a':
                str[mStr.length()-c-1] = 'T';
                break;
            case 'T':
            case 't':
                str[mStr.length()-c-1] = 'A';
                break;
            case 'C':
            case 'c':
                str[mStr.length()-c-1] = 'G';
                break;
            case 'G':
            case 'g':
                str[mStr.length()-c-1] = 'C';
                break;
            default:
                str[mStr.length()-c-1] = 'N';
        }
    }
    return Sequence(str);
}

Sequence Sequence::operator~(){
    return reverseComplement();
}

bool Sequence::test(){
    Sequence s("AAAATTTTCCCCGGGG");
    Sequence rc = ~s;
    if (s.mStr != "AAAATTTTCCCCGGGG" ){
        cerr << "Failed in reverseComplement() expect AAAATTTTCCCCGGGG, but get "<< s.mStr;
        return false;
    }
    if (rc.mStr != "CCCCGGGGAAAATTTT" ){
        cerr << "Failed in reverseComplement() expect CCCCGGGGAAAATTTT, but get "<< rc.mStr;
        return false;
    }
    return true;
}


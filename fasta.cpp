#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <vector>
#include <filesystem>
#include <iostream>
using namespace std;

//http://lh3lh3.users.sourceforge.net/parsefastq.shtml
// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)

void append_clean_sequences(string& sequence_name, string& sequence, vector<string>& sequence_names, vector<string>& sequences){
    // Make the string several sequences.
    // We iterate through the string and make everything lower.
    // Whenever we find a char that has to be excluded we copy
    //          the substring to the sequence vector and
    //          the name to the name vector
    // The reason we give a reference is that we want to save memory and don't care about the string getting altered.
    char* start = &sequence[0];
    char* sequence_end = &sequence[sequence.size()];
    int counter = 0;
    char* iter;
    for(iter = &sequence[0]; iter < sequence_end; ++iter){
        // Make it a lower char
        *(iter) = tolower(*(iter));
        if(*iter == 'a' | *iter == 'c' | *iter == 'g' | *iter == 't'){
            // Do you have to do anything?
            // Idea: If we have a 'u' for example, maybe we could add both sequences.
        }
        else{
            if(iter > start) {
                sequences.push_back(sequence.substr(start-&sequence[0], iter - start));
                sequence_names.push_back(sequence_name + "%" + to_string(counter));
                ++counter;
            }
            start = iter + 1; // Is this really the way to go?
        }
    }
    // If we are finished and the last character is okay then we should add the last sequence.
    --iter;
    if(*iter == 'a' | *iter == 'c' | *iter == 'g' | *iter == 't'){
        sequences.push_back(sequence.substr(start-&sequence[0], iter - start+1));
        sequence_names.push_back(sequence_name + "%" + to_string(counter));
    }
}

void read_fasta_into_vectors(string path, vector<string>& names, vector<string>& sequences)
{
    gzFile fp;
    kseq_t *seq;
    int l;

    fp = gzopen(&path[0], "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    cout << path << endl;
    string name, sequence;
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        // --> seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s seems to be the content
        // --> seq->comment.l seems to understand whether it exists (for comment and qual)

        // -> Read the sequence; Write the stuff into a vector of strings; also keep track of the sequence ID
        name = seq -> name.s;
        sequence = seq -> seq.s;
        append_clean_sequences(name,sequence,names,sequences);
        // Ignore NNN and (not polymorphism, but sounds similar); also make everything small letters
        // -> Return the vectors
        //

    }
    // Return value: l (should be -1 probably - guess it is the value where reading stopped)
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
}


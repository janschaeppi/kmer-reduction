//
// Created by Jan Schäppi on 21.09.21.
//
#include <vector>
using namespace std;

//https://www.programcreek.com/cpp/?CodeExample=murmur+hash64
uint32_t MurmurHash (const char *key, uint32_t len, uint32_t seed) {
    uint32_t c1 = 0xcc9e2d51;
    uint32_t c2 = 0x1b873593;
    uint32_t r1 = 15;
    uint32_t r2 = 13;
    uint32_t m = 5;
    uint32_t n = 0xe6546b64;
    uint32_t h = 0;
    uint32_t k = 0;
    uint8_t *d = (uint8_t *) key; // 32 bit extract from `key'
    const uint32_t *chunks = NULL;
    const uint8_t *tail = NULL; // tail - last 8 bytes
    int i = 0;
    int l = len / 4; // chunk length

    h = seed;

    chunks = (const uint32_t *) (d + l * 4); // body
    tail = (const uint8_t *) (d + l * 4); // last 8 byte chunk of `key'

    // for each 4 byte chunk of `key'
    for (i = -l; i != 0; ++i) {
        // next 4 byte chunk of `key'
        k = chunks[i];

        // encode next 4 byte chunk of `key'
        k *= c1;
        k = (k << r1) | (k >> (32 - r1));
        k *= c2;

        // append to hash
        h ^= k;
        h = (h << r2) | (h >> (32 - r2));
        h = h * m + n;
    }

    k = 0;

    // remainder
    switch (len & 3) { // `len % 4'
        case 3: k ^= (tail[2] << 16);
        case 2: k ^= (tail[1] << 8);

        case 1:
            k ^= tail[0];
            k *= c1;
            k = (k << r1) | (k >> (32 - r1));
            k *= c2;
            h ^= k;
    }

    h ^= len;

    h ^= (h >> 16);
    h *= 0x85ebca6b;
    h ^= (h >> 13);
    h *= 0xc2b2ae35;
    h ^= (h >> 16);

    return h;
}

int LexicographicHash(char* first, int length){
    int sum = 0;
    for(char* first_iter = first; first_iter < first + length; ++first_iter) sum += int(*first_iter);
    return sum;
}

void minimizer_pos_bruteforce(char*& minimizer_ptr, int& minimizer_val, char* sequence, int sequence_len, int l){
    minimizer_val = MurmurHash(sequence,l,12345);
    minimizer_ptr = sequence;
    int candidate_minimum;
    for(char* sequence_iter = sequence+1; sequence_iter < sequence + sequence_len - l + 1; ++sequence_iter){
        candidate_minimum = MurmurHash(sequence_iter,l,12345);
        if(candidate_minimum < minimizer_val){
            minimizer_val = candidate_minimum;
            minimizer_ptr = sequence_iter;
        }
    }
}

void collect_syncmers(vector<int>& sequence_positions, vector<int>& syncmer_positions, vector<int>& sequence_ids,string& sequence, int sequence_id, int s, int k, int l){

    // s : Subsequence length
    // k : syncmer length
    // l : minimizer length
    // Make sure subsequence &sequence[0]+sequence.length()-subsequence_len+1 makes sense

    char* sequence_ptr = &sequence[0];
    char* minimum_kmer_ptr = sequence_ptr-1;
    char* last_minimum_kmer_ptr = sequence_ptr-1;
    char* minimum_lmer_ptr = sequence_ptr-1; int minimum_lmer_val;
    char* last_lmer_ptr; int last_lmer_val;
    //minimizer_pos_bruteforce(minimizer_ptr,minimizer_val, sequence_ptr,subsequence_len, l);

    for(; sequence_ptr < &sequence[0]+sequence.length()-s+1; ++sequence_ptr){
        // SET THE MINIMUM LMER POINTER
        // Case 1: minimizer_ptr < sequence_ptr, then brute force again
        if(minimum_lmer_ptr < sequence_ptr) minimizer_pos_bruteforce(minimum_lmer_ptr,minimum_lmer_val, sequence_ptr,s, l);
        // Case 2: sequence_ptr + subsequence_len - l has a smaller value than minimizer_val -> replace
        else{
            last_lmer_ptr = sequence_ptr + s - l;
            last_lmer_val = MurmurHash(last_lmer_ptr,l,12345);
            if(last_lmer_val < minimum_lmer_val){
                minimum_lmer_ptr = last_lmer_ptr;
                minimum_lmer_val = last_lmer_val;
            }
        }
        // Case 3: else just continue

        // SET THE MINIMUM KMER POINTER
        // If in the range then the syncmer starts at the minimum lmer
        // Else it starts at the last k-mer
        if(minimum_lmer_ptr + k < sequence_ptr + s) minimum_kmer_ptr = minimum_lmer_ptr;
        else minimum_kmer_ptr = sequence_ptr + s - k;

        // ADD THE POINTER POSITION TO THE VECTOR IF NECESSARY; ALSO ADD THE SEQUENCE POSITION
        if(minimum_kmer_ptr != last_minimum_kmer_ptr){
            syncmer_positions.push_back(minimum_kmer_ptr-&sequence[0]);
            sequence_positions.push_back(sequence_ptr-&sequence[0]);
            sequence_ids.push_back(sequence_id);
            last_minimum_kmer_ptr = minimum_kmer_ptr;
            //cout << "Changed Minimizer to position " <<  minimum_kmer_ptr-&sequence[0] << endl;
        }
        //else cout << "Retained Minimizer." << endl;

    }
}
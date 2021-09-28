#include <iostream>
#include "GenomeIndexer.cpp"

int main() {

    GenomeIndexer G(7,5,3);
    G.add_file("/Users/janschaeppi/Scripts/masterthesis/kmer-distribution/testfile_merged.fna");

    G.hash_all();

    vector<int> distribution(30);
    G.hash_countup(distribution);

    print_hotizontal(distribution);

    G.stoptime_journal();
    return 0;
}

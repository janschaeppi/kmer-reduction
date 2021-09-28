//
// Created by Jan Schäppi on 22.09.21.
//

#include "fasta.cpp"
#include "syncmers.cpp"
#include <filesystem>
#include "utils.cpp"
#include <iostream>
#include <chrono>
using namespace std;

struct GenomeIndexer {
    // PARAMETERS
    int s; // The window size with which the sequences will be parsed
    int k; // The syncmer length
    int l; // The l-mer length for the syncmer
    string path; // Path of a directory from which all fasta files will be included

    vector<chrono::steady_clock::time_point> stoptimes;
    vector<string> stoptime_explanations;

    // FASTA READING
    const string fasta_endings[12] = {".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa", ".fasta.gz", ".fna.gz",
                                      ".ffn.gz", ".faa.gz", ".frn.gz", ".fa.gz"};

    // INPUT MATRIX - length x
    vector<string> sequences;
    vector<string> sequence_names;

    // SYNCMER HASHED MATRIX - length y
    vector<int> sequence_ids;
    vector<int> rel_sync_sequence_start;
    vector<int> rel_sync_syncmer_start;

    // HASHTABLE
    int hashtable_size;
    int murmur_hash_seed = 352632;
    vector<vector<int>> hashtable;

    GenomeIndexer(int _s, int _k, int _l) {
        s = _s;
        k = _k;
        l = _l;
        stoptime("Genome Indexer initialized.");
    }
    void status(const string message, const int i, const int tot) {
        // We want about 100 status prints. Hence tot / 100 and go to the closest multiple of 10.
        // If int(tot / 10) == 0 then we shall print in any case.
        int mod = int(tot / 100) - int(tot / 100) % 10;
        if (mod == 0 || i % mod == 0) {
            cout << "Status: (" << i << "/" << tot << ") " << message << "\r" << flush;
        }
        if (tot == i + 1) cout << "Status: (Finished) " << message << endl;
    }
    void stoptime(const string message){
        stoptimes.push_back(chrono::high_resolution_clock::now());
        stoptime_explanations.push_back(message);
    }
    void stoptime_journal(){
        std::chrono::duration<double> diff_total, diff_to_last;

        for(int i = 1; i < stoptimes.size(); ++i){
            diff_total = stoptimes[i]-stoptimes[0];
            diff_to_last = stoptimes[i] - stoptimes[i-1];
            cout << " Total (since init): " << diff_total.count() << "s ";
            cout << " Compared to last: " << diff_to_last.count() << "s ";
            cout << stoptime_explanations[i] << endl;
        }
    }
    bool is_fasta(const string path) {
        for (const auto &fasta_ending: fasta_endings) {
            if (endswith(path, fasta_ending)) return true;
        }
        return false;
    }
    void add_directory(const string path) {
        for (const auto &entry: __fs::filesystem::directory_iterator(path)) {
            if (is_fasta(entry.path())) {
                add_file(entry.path());
            }
        }
    }
    void add_file(const string path) {
        stoptime("Adding file " + path);
        const string status_msg = "Parsing Fasta file" + path;
        const int x_before = sequences.size();
        read_fasta_into_vectors(path, sequence_names, sequences);
        for (int x = x_before; x < sequences.size(); ++x) {
            collect_syncmers(rel_sync_sequence_start, rel_sync_syncmer_start, sequence_ids, sequences[x], x, s, k, l);
            status(status_msg, sequences.size() - x, sequences.size() - x_before);
        }
        stoptime("Adding file done.");
    }
    void hash_all() {

        stoptime("Hashing Start.");
        hashtable_size = sequence_ids.size(); // Only for now
        hashtable.clear();
        hashtable.resize(hashtable_size);

        string syncmer;
        for (int index = 0; index < sequence_ids.size(); ++index) {
            // syncmer start is syncmer_positions[index]
            // syncmer is substr(syncmer_start,k)
            syncmer = sequences[sequence_ids[index]].substr(rel_sync_sequence_start[index], k);
            hashtable[MurmurHash(&syncmer[0], k, murmur_hash_seed) % hashtable_size].push_back(index);
        }
        stoptime("Hashing Done.");
    }

    void print_hastable() {
        for (int i = 0; i < hashtable.size(); ++i) {
            cout << i << ": [";
            for (int j = 0; j < hashtable[i].size(); ++j) {
                cout << hashtable[i][j];
                if (j < hashtable[i].size() - 1) cout << ", ";
            }
            cout << "]" << endl;
        }
    }
    void explain_hashtable_entry(int index){
        int sequence_id, syncmer_seqid;
        cout << "Hashtable contains indices ";
        for(int i = 0; i < hashtable[index].size(); ++i) cout << hashtable[index][i] << " ";
        cout << endl << "which are the indices in the Y-Matrix (describing the syncmer sequences), so sequence_ids, rel_sync_sequence_start and rel_sync_syncmer_start.";
        cout << "For these we have the following entries:"  << endl;
        for(int i = 0; i < hashtable[index].size(); ++i){
            cout << "---------------------------------" << endl;
            syncmer_seqid = hashtable[index][i];
            cout << "Syncmer Sequence index " << syncmer_seqid << endl;
            cout << endl;
            cout << "sequence_ids[" << syncmer_seqid << "] = " << sequence_ids[syncmer_seqid] << " (Entry position in the X-Matrix describing fasta sequences)" << endl;
            cout << "rel_sync_sequence_start[" << syncmer_seqid << "] = " << rel_sync_sequence_start[syncmer_seqid] << " (Position in the fasta sequence where the syncmer sequence starts)" << endl;
            cout << "rel_sync_syncmer_start[" << syncmer_seqid << "] = " << rel_sync_syncmer_start[syncmer_seqid] << " (Position in the fasta sequence where the syncmer starts)" << endl;
            cout << endl;
            sequence_id = sequence_ids[syncmer_seqid];
            cout << "sequence_names[" << sequence_id << "] = " << sequence_names[sequence_id] << endl;
            cout << "sequences[" << sequence_id << "] = " << sequences[sequence_id].substr(0,30) << ".." << endl;
            cout << endl;
            cout << "Associated Syncmer: " << sequences[sequence_id].substr(rel_sync_syncmer_start[syncmer_seqid],k) << endl;
            cout << "Associated Syncmer Sequence: " << sequences[sequence_id].substr(rel_sync_sequence_start[syncmer_seqid],s) << endl;
        }
    }

    void print_hashtable_minimizers(int index){
        int sequence_id;
        for(const auto& syncmer_seqid : hashtable[index]){
            sequence_id = sequence_ids[syncmer_seqid];
            cout << sequences[sequence_id].substr(rel_sync_syncmer_start[syncmer_seqid],k) << endl;
        }
    }
    void hash_countup(vector<int>& distribution){
        stoptime("Hash Countup start.");
        vector<vector<int>> hash_counter(hashtable_size);

        for(int index = 0; index < hashtable.size(); ++index){
            numberOfEqualSMers(distribution,index);
            status("Hashtable Countup",index+1,hashtable.size());
        }
        stoptime("Hash Countup end.");
    }

    int sync_sequence_length(const int sequence_id, const int syncmer_seqid){
        if(syncmer_seqid + 1 == rel_sync_sequence_start.size() || rel_sync_sequence_start[syncmer_seqid + 1] == 0){
            return (sequences[sequence_id].size() - rel_sync_sequence_start[syncmer_seqid]);
        }
        else{
            return rel_sync_sequence_start[syncmer_seqid+1] - rel_sync_sequence_start[syncmer_seqid];
        }
    }

    void numberOfEqualSMers(vector<int>& distribution, int index){

        // INITIALIZATION STEP 1: CHECK HOW MANY KMERS WE NEED TO COMPARE
        int number_of_competitors = hashtable[index].size();

        // INITIALIZATION STEP 2: POINTERS
        // (char*) first_first_of_sequence[i] <= runner[i] <= last_first_of_sequence[i] in the valid region
        // Example: (Sub)sequence ACGT, s=2 => *first_first_of_sequence = A; *last_first_of_sequence = G
        // We set the runners such that they reach the syncmers at the same time in the loop.
        // For this we calculate the maximum distance from sequence start to syncmer start (max) and set the runner to syncmer position - max
        // Example: (Sub)sequences AACCCTT and CCCTTGG with syncmer CCC, we want to align as follows
        //          AACCCTT
        //            CCCTTGG
        // so runner[0] = first_first_of_sequence[0] and runner[1] = first_first_of_sequence[1] - 2.

        vector<char*> runner(number_of_competitors);
        vector<char*> first_first_of_sequence(number_of_competitors);
        vector<char*> last_first_of_sequence(number_of_competitors);

        // INITIALIZATION 2.1: MAXIMUM DISTANCE FROM SYNCMER_SEQUENCE START TO SYNCMER START TO SYNCHRONIZE THE RUNNERS
        int max_distance_to_syncmer = 0, distance_to_syncmer;
        int syncmer_seqid, sequence_id; // will only be used in the first two for-loops
        for(int i = 0; i < number_of_competitors; ++i){
            syncmer_seqid = hashtable[index][i];
            distance_to_syncmer = rel_sync_syncmer_start[syncmer_seqid] - rel_sync_sequence_start[syncmer_seqid];
            if(distance_to_syncmer > max_distance_to_syncmer){
                max_distance_to_syncmer = distance_to_syncmer;
            }
        }

        // INITIALIZATION 2.2: SETTING RUNNER AND ITS BOUNDARIES
        // The end of a syncmer_sequence can always be determined by the beginning of the next syncmer_sequence.
        // Iff the next beginning is a zero or the nonexistent then the syncmer_sequence is the end of the entire sequence.
        for(int i = 0; i < number_of_competitors; ++i){
            syncmer_seqid = hashtable[index][i];
            sequence_id = sequence_ids[syncmer_seqid];
            first_first_of_sequence[i] = &(sequences[sequence_id][rel_sync_sequence_start[syncmer_seqid]]);
            if(syncmer_seqid + 1 == rel_sync_sequence_start.size() || rel_sync_sequence_start[syncmer_seqid + 1] == 0){
                // last: s-mer that ends at the sequence end. ACGT, s=2 we get 0 (first) + 4 (size - sync_sequence_start) - 2 (s) = 2
                // so the last first is A+2 = G (in pointer terms)
                last_first_of_sequence[i] = first_first_of_sequence[i] + (sequences[sequence_id].size() - rel_sync_sequence_start[syncmer_seqid] - s);
            }
            else{
                // last: s-mer that begins right before the next syncmer_sequence.
                // this is okay since we know that the next syncmer_sequence has size s, which means the sequence is still valid at last_first + s.
                last_first_of_sequence[i] = first_first_of_sequence[i] + (rel_sync_sequence_start[syncmer_seqid+1] - rel_sync_sequence_start[syncmer_seqid] - 1);
            }
            runner[i] = first_first_of_sequence[i] + (rel_sync_syncmer_start[syncmer_seqid]-rel_sync_sequence_start[syncmer_seqid])-max_distance_to_syncmer;
        }

        // SHORTCUT: IF THE NUMBER OF VALUES IS 0 OR 1, WE CAN SHORTCUT
        if(number_of_competitors == 0) return;
        if(number_of_competitors == 1){
            int number_of_s_mers;
            syncmer_seqid = hashtable[index][0];
            sequence_id = sequence_ids[syncmer_seqid];
            if(syncmer_seqid + 1 == rel_sync_sequence_start.size() || rel_sync_sequence_start[syncmer_seqid + 1] == 0){
                // as the sequence terminates here, we can count up everything that ENDS BEFORE the sequence end.
                // ACGT, s=2 -> 4 (size) - 0 (start) - 2 (s) + 1 = 3, which is AC, CG, GT.
                number_of_s_mers = sequences[sequence_id].size() - rel_sync_sequence_start[syncmer_seqid] - s + 1;
            }
            // as the sequence does not terminate and the algorithm has found another sequence of at least size s
            // afterwards, we can cound every s-mer starting in the between this and the next syncmer_sequence.
            else number_of_s_mers = rel_sync_sequence_start[syncmer_seqid+1] - rel_sync_sequence_start[syncmer_seqid];
            distribution[1] += number_of_s_mers;
            return;
        }

        // INITIALIZATION 3: PREPARING STATE VARIABLES
        vector<bool> match = strictLowerDiagonalMatrix(number_of_competitors,false); // strictly lower diagonal matrix; sparse matrix storing whether i and j match
        vector<int> checksums(number_of_competitors,0); // [sum(int c) for c in str]. Checksums are always bigger than 0.
        vector<bool> runner_in_valid_position(number_of_competitors); // runners can go out of their proper boundaries; in that case this will be false and i will be ignored.
        vector<bool> cluster_found(number_of_competitors,false); // important for evaluation of the match matrix.
        int cluster_size; // important for evaluation of the match matrix.

        // LOOPING OVER SEQUENCES CHARWISE UNTIL EVERY RUNNER RUNS OUT OF SCOPE
        bool has_valid_runners; // the while loop stopper.
        while(true){

            // A VALIDITY AND CHECKSUM
            // Set correctly: runner_in_valid_position, checksum, has_valid_runners
            // checksum is only set for runners in valid position,
            has_valid_runners = false;
            for(int i = 0; i < number_of_competitors; ++i){
                runner_in_valid_position[i] = first_first_of_sequence[i] <= runner[i] && runner[i] <= last_first_of_sequence[i];
                has_valid_runners = has_valid_runners || runner_in_valid_position[i];
                if(runner_in_valid_position[i]){
                    if(checksums[i] == 0) checksums[i] = LexicographicHash(runner[i],s);
                    else checksums[i] = checksums[i] - int(*runner[i]) + int(*(runner[i]+s-1));
                }
            }
            // has_valid_runners is now set correctly.
            // If we have no valid runner, all runners have run out of scope. Stop the loop.
            if(!has_valid_runners) break;

            // B FILLING THE MATCH MATRIX
            // Filling match such that after the loop it is correct wherever the runners of i and j are in valid position.
            // The goal of many of the preparatory values was to avoid having to iterate over the k-mers to check
            //      whether they are the same.
            for(int i = 0; i < number_of_competitors; ++i){
                if(runner_in_valid_position[i]){
                    for(int j = 0; j < i; ++j){
                        if(runner_in_valid_position[j]){
                            // i, j are valid runners. We check whether the s-mers at runner[i] and runner[j] match.
                            // Shortcut: If checksums do not match it is impossible that i and j match.
                            if(checksums[i] == checksums[j]){
                                // Shortcut: If checksum is the same and the s-mers have matched before, it is still a match.
                                if(match[positionInStrictLowerDiagonalMatrix(i,j)] == true); // no change in match necessary
                                // Necessary: If the checksums have not matched previously, we check brute force.
                                else if(equals(runner[i],runner[j],s)) match[positionInStrictLowerDiagonalMatrix(i,j)] = true;
                                // If none of the above cases is true, then match(i,j) was previously false and it should stay false.
                            }
                            else{
                                match[positionInStrictLowerDiagonalMatrix(i,j)] = false;
                            }
                        }
                    }
                }
            }

            // EVALUATING THE MATCH MATRIX
            // A "cluster" here is a set of equivalent s-mers.
            // At every position i (starting with the highest one as in the strict lower diagonal matrix
            //      this one has information on all its matches), we create it a new cluster (cluster_size = 1),
            //      add all matches to the cluster (++cluster_size) and
            //      do not allow them to be their own cluster anymore (cluster_found[i] = true)
            // We add every cluster to the distribution distribution[cluster_size] is the number of clusters of that size.
            // The last entry in distribution is the tail (cluster size too big).
            for(int i = 0; i < number_of_competitors; ++i) cluster_found[i] = false;
            for(int i = number_of_competitors-1; i >= 0; --i){
                cluster_size = 1;
                if(cluster_found[i] || !runner_in_valid_position[i]) continue;
                cluster_found[i] = true;
                for(int j = 0; j < i; ++j){
                    if(match[positionInStrictLowerDiagonalMatrix(i,j)] & runner_in_valid_position[j]){
                        ++cluster_size;
                        cluster_found[j] = true;
                    }
                }

                if(cluster_size >= distribution.size()-1) ++distribution[distribution.size()-1];
                else ++distribution[cluster_size];
            }
            // (while loop over)

            // D INCREMENTING THE RUNNERS
            for(int i = 0; i < number_of_competitors; ++i) ++runner[i];
        }
    }

};

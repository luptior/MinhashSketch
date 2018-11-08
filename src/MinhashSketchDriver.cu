// this is host-side app
// modified from MinhashSketch(https://github.com/daren996/MinhashSketch)

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <cstdint>
#include <random>
#include <vector>
#include <string>

#include "MinhashSketch.cu"
#include "Hash.h"
#include "Utils.h"

using namespace std;

void output_signature(vector <vector<uint64>> sig1) {
    for (int h = 0; h < sig1.size(); ++h) {
        cout << "sig[" << h << "].size(): " << sig1[h].size() << "\t";
        for (int i = 0; i < sig1[h].size(); ++i) {
            cout << hex << sig1[h][i] << dec << " ";
        }
        cout << endl;
    }
}

void usage() {
    cout << "===========================" << endl;
    cerr << "Usage: " << endl << endl;
    cerr << "    ./MinhashSketch FILE_ONE FILE_TWO MODE" << endl;
    cerr << endl;
    cerr << "    Possible MODEs are:" << endl;
    cerr << endl << bold_on;
    cerr << "        all" << endl;
    cerr << endl;
    cerr << "        minhash_parallel" << endl;
    cerr << endl;
    cerr << "Execute \"MinhashSketch help\" for an extended help section." << endl;
    cout << "===========================" << endl;
    exit(1);
}

void help() {
    cout << endl;
    cout << bold_on << "NAME" << bold_off << endl;
    cout << "    " << "MinhashSketch" << endl;
    cout << endl;
    cout << bold_on << "USAGE" << bold_off << endl;
    cout << "    " << "MinhashSketch FILE_ONE FILE_TWO " << bold_on << "MODE [PARAMETERS...]" << bold_off << endl;
    cout << endl;
    cout << "    " << "MinhashSketch calculates the similarity between two text files FILE_ONE and FILE_TWO" << endl;
    cout << "    " << "and outputs it as a number between 0 and 1, where 1 means the two files are exactly" << endl;
    cout << "    " << "the same." << endl;
    cout << endl;
    cout << bold_on << "MODE" << bold_off << endl;
    cout << "    " << "There are modes which change the way MinhashSketch computes the similarity. " << endl;
    cout << "    " << "Each may make use of different parameters, indicated as follows:" << endl;;
    cout << endl;
    cout << "    " << bold_on << "all" << bold_off << endl;
    cout << "        " << "This option executes all modes." << endl;
    cout << endl;
    cout << "    " << bold_on << "minhash" << bold_off << endl;
    cout << "        " << "Calculates the similarity by computing minhash signatures for each sequence. Used" << endl;
    cout << "        " << "parameters are." << endl;
    cout << endl;
    cout << "            " << "--k=POSITIVE_INTEGER as shingle size" << endl;
    cout << endl;
    cout << "            " << "--t=POSITIVE_INTEGER" << bold_on << " (obligatory) " << bold_off
         << "as number of hash functions used" << endl;
    cout << endl;
    cout << "            " << "--seed=INTEGER as random generator seed" << endl;
    cout << endl;
    cout << bold_on << "PARAMETERS" << bold_off << endl;
    cout << endl;
    cout << "    " << bold_on << "--k=POSITIVE_INTEGER" << bold_off << endl;
    cout << "        " << "Defaults to k=9. Indicates the size of the shingles used to calculate the simi-" << endl;
    cout << "        " << "larity between the documents." << endl;
    cout << endl;
    cout << "    " << bold_on << "--m=POSITIVE_INTEGER" << bold_off << endl;
    cout << "        " << "Defaults to m=1. Indicates the number of sketches saved in minhash modes." << endl;
    cout << endl;
    cout << "    " << bold_on << "--t=POSITIVE_INTEGER" << bold_off << endl;
    cout << "        " << "Defaults to t=1. Indicates the number of hash functions used in minhash modes." << endl;
    cout << endl;
    cout << "    " << bold_on << "--seed=INTEGER" << bold_off << endl;
    cout << "        " << "Defaults to a random value. Used by minhash modes in their random generator number." << endl;
    cout << endl;
    cout << "    " << bold_on << "-e" << bold_off << endl;
    cout << "        " << "Output in experimentation format." << endl;
    cout << endl;
    exit(0);
}

// Main function to run the program in mercator style
int run_mercator(const int k, const int m, const int t, char *dnaList, int length, uint64 *hashes_b, uint64_t* result)
{   

    MinhashSketch sketch;

    // number of blocks
    const int BLOCKS_NUM = sketch.getNBlocks();
    const int BLOCK_THREADS = 32 * 16;
    const int ITEMS_PER_THREAD = 4;
    cout << "Number of Blocks" << BLOCKS_NUM << "\n";
    int CHUNKS_NUM;
    int CHUNK_SIZE=BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD ;

    // Calculate how many chuncks does one have
    if (length % (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) == 0)
        CHUNKS_NUM = (length - k + 1) / CHUNK_SIZE;
    else
        CHUNKS_NUM = (length - k + 1) / CHUNK_SIZE + 1;

    for (int i = 0; i < (CHUNKS_NUM - 1); ++i) {
        // host_side buffers
        char *input = new char [CHUNK_SIZE];
        input = &dnaList + i*CHUNKS_SIZE ;
        uint64_t *output = new uint64_t [m];

        // begin MERCATOR usage
        Mercator::Buffer<char> ib(CHUNK_SIZE);
        Mercator::Buffer<uint64_t> ob(m);
        
        // move data into the input buffer
        ib.set(input, length);
        
            // pass in the program parameters
        sketch.getParams()->m=m;
        sketch.getParams()->k=k;
        sktech.BlockGetSketch.getParam()->thread_offset = """something_ add here""";
        sktech.BlockGetSketch.getParam()->hash_b = """something_ add here""";


        // set input and output place
        sketch.src.setSource(ib);
        sketch.snk.setSink(ob);
        
        sketch.run();
        ob.get(output, ob.size());
    
        //Print Results
    }

}

// MinhashSketch.exe ../testing_files/sequence_clip1.fasta ../testing_files/sequence_clip2.fasta all -e --k=5 --m=10 --t=10
int main(int argc, char *argv[]) {

    if (argc == 2 && string(argv[1]) == "help") help();
    if (argc < 4) usage();

    // DEFAULT VALUES
    int k, m, t, seed;
    bool e;
    k = 9;
    m = 1; // the number of sketches
    t = 1; // the number of hash functions
    seed = random_device()();
    e = false;

    // PARSE FILE_ONE FILE_TWO MODE
    string name_one = string(argv[1]);
    string name_two = string(argv[2]);
    string cal_name = string(argv[3]);
    ifstream file1(name_one);
    if (file1.fail()) {
        std::cerr << "Unable to open file 1" << name_one << std::endl;
        exit(1);
    }
    ifstream file2(name_two);
    if (file2.fail()) {
        std::cerr << "Unable to open file 2" << name_two << std::endl;
        exit(1);
    }

    // PARSE PARAMETERS
    for (int i = 4; i < argc; ++i) {
        string param(argv[i]);
        if (param == "-e") {
            e = true;
        } else {
            int param_size = (uint) param.size();
            if (param_size >= 5) {
                auto index_eq = (uint) param.find('=');
                if (index_eq + 2 <= param_size) {
                    string param_name = param.substr(0, index_eq);
                    if (param_name == "--k") {
                        k = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    } else if (param_name == "--m") {
                        m = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    } else if (param_name == "--t") {
                        t = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    } else if (param_name == "--seed") {
                        seed = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    }
                }
            }
        }
    }
    if (k < 1) {
        std::cerr << "K value too small! Minimum: 1" << std::endl;
        exit(1);
    }
    if (m < 1) {
        std::cerr << "M value too small! Minimum: 1" << std::endl;
        exit(1);
    }

    // GET TWO SEQUENCES
    string file_info1, file_info2, sequence1, sequence2, s1, s2;
    utils::file_to_string(file1, file_info1, sequence1); 
    utils::file_to_string(file2, file_info2, sequence2);
    uint64 sequence_size1 = sequence1.size(), sequence_size2 = sequence2.size();
    if (sequence1.size() < k || sequence2.size() < k) {
        cout << "k cannot be greater than the size of any document" << endl;
        exit(1);
    }
    file1.close();
    file2.close();
    cout << file_info1 << "\n" <<"sequence1 size: " << sequence1.size() << endl;
    cout << file_info2 << "\n" <<"sequence2 size: " << sequence2.size() << endl;
    // dnaLists are char array store the sequence in chars
    char dnaList1[sequence1.size()];
    char dnaList2[sequence2.size()];
    strcpy(dnaList1, sequence1.c_str());
    strcpy(dnaList2, sequence1.c_str());

    // MAIN PROGRESS
    clock_t ini_time;
    bool mode_found = false;
    double similarity, time;
    list <tuple<string, double, double>> results;

    //a list of hash randoms are calculated for t hash functions
    uint64 *hashes_b = generateHashes_b(t, seed);


    if (cal_name == "all" || cal_name == "minhash_parallel") {

        if (t < 1) {
            cerr << endl;
            cerr << "You must provide a parameter --t=POSITIVE_INTEGER parameter for minhash modes!" 
                    << endl << endl;
            exit(1);
        }

        mode_found = true;
        ini_time = clock();

        sig2  = new uint64_t[m];
        sig1  = new uint64_t[m];

        //
        // Changes should be made here to use Mercator version of codes
        run_mercator(k, m, t, dnaList1, sequence1.size(), hashes_b, sig1);
        run_mercator(k, m, t, dnaList2, sequence2.size(), hashes_b, sig2);
        
        // Change the original output type to be fixed size array
        // vector <vector<uint64>> sig1 = genSig(k, m, t, dnaList1, sequence1.size(), hashes_b);
        // vector <vector<uint64>> sig2 = genSig(k, m, t, dnaList2, sequence2.size(), hashes_b);


        cout << "sig1 size:" << sig1[0].size() << endl;
        output_signature(sig1);
        cout << "\nsig2 size:" << sig2[0].size() << endl;
        output_signature(sig2);
        cout << endl;
        similarity = computeSim(sig1, sig2);
        time = double(clock() - ini_time) / CLOCKS_PER_SEC;
        results.emplace_back("minhash_parallel", similarity, time);
    }
    if (!mode_found) usage();

    // OUTPUT RESULTS
    if (e) {
        cout << setw(12) << "cal_name" << setw(14) << "seed" << setw(9) << "k" << setw(5) << "m" << setw(7) << "t";
        cout << setw(9) << fixed << "time" << setw(13) << fixed << "similarity" << endl;
    } else {
        cout << "===========================" << endl;
        cout << "k:" << k << setw(7) << fixed << "m:" << m << setw(7) << fixed << "t:" << t << endl;
        cout << "===========================" << endl;
    }
    cout.precision(8);
    for (auto &result : results) {
        if (e) {
            cout << setw(12) << get<0>(result) << setw(14) << seed << setw(5) << k <<9setw(5) << m << setw(7) << t;
            cout << setw(9) << fixed << get<2>(result) << setw(13) << fixed << get<1>(result) << endl;
        } else {
            cout << uline_on << get<0>(result) << uline_off << endl;
            cout << "time: " << setw(21) << fixed << get<2>(result) << endl;
            cout << "similarity: " << setw(15) << fixed << get<1>(result) << endl;
            cout << "===========================" << endl;
        }
    }

}
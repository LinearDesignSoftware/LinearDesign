#include <iomanip>
#include "beam_cky_parser.h"
#include "beam_cky_parser.cc"
#include "Utils/reader.h"
#include "Utils/common.h"
#include "Utils/codon.h"

// #ifndef CODON_TABLE
// #define CODON_TABLE "./codon_usage_freq_table_human.csv"
// #endif

#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif

using namespace LinearDesign;

template <typename ScoreType, typename IndexType>
bool output_result(const DecoderResult<ScoreType, IndexType>& result, 
        const double duration, const double lambda, const bool is_verbose, 
        const Codon& codon, string& CODON_TABLE) {

    stringstream ss;
    if (is_verbose)
        ss << "Using lambda = " << (lambda / 100.) << "; Using codon frequency table = " << CODON_TABLE << endl;
    ss << "mRNA sequence:  " << result.sequence << endl;
    ss << "mRNA structure: " << result.structure << endl;
    ss << "mRNA folding free energy: " << std::setprecision(2) << fixed << result.score 
                                        << " kcal/mol; mRNA CAI: " << std::setprecision(3) 
                                        << fixed << codon.calc_cai(result.sequence) << endl;
    if (is_verbose)
        ss << "Runtime: " << duration << " seconds" << endl;
    cout << ss.str() << endl;

    return true;
}

void show_usage() {
    cerr << "echo SEQUENCE | ./lineardesign -l [LAMBDA]" << endl;
    cerr << "OR" << endl;
    cerr << "cat SEQ_FILE_OR_FASTA_FILE | ./lineardesign -l [LAMBDA]" << endl;
}


int main(int argc, char** argv) {

    // default args
    double lambda = 0.0f;
    bool is_verbose = false;
    string CODON_TABLE = "./codon_usage_freq_table_human.csv";

    // parse args
    if (argc != 4) {
        show_usage();
        return 1;
    }else{
        lambda = atof(argv[1]);
        is_verbose = atoi(argv[2]) == 1;
        if (string(argv[3]) != ""){
            CODON_TABLE = argv[3];
        }
    } 
    lambda *= 100.;
    
    // load codon table and coding wheel
    Codon codon(CODON_TABLE);
    std::unordered_map<string, Lattice<IndexType>> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    prepare_codon_unit_lattice<IndexType>(CODING_WHEEL, codon, aa_graphs_with_ln_weights, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon, lambda);

    // main loop
    string aa_seq, aa_tri_seq;
    vector<string> aa_seq_list, aa_name_list;
    // load input
    for (string seq; getline(cin, seq);){
        if (seq.empty()) continue;
        if (seq[0] == '>'){
            aa_name_list.push_back(seq); // sequence name
            if (!aa_seq.empty())
                aa_seq_list.push_back(aa_seq);
            aa_seq.clear();
            continue;
        }else{
            rtrim(seq);
            aa_seq += seq;
        }
    }
    if (!aa_seq.empty())
        aa_seq_list.push_back(aa_seq);

    // start design
    for(int i = 0; i < aa_seq_list.size(); i++){
        if (aa_name_list.size() > i)
            cout << aa_name_list[i] << endl;
        auto& aa_seq = aa_seq_list[i];
        // convert to uppercase
        transform(aa_seq.begin(), aa_seq.end(), aa_seq.begin(), ::toupper);
        aa_tri_seq.clear();
        if (is_verbose)
            cout << "Input protein: " << aa_seq << endl;
        if (!ReaderTraits<Fasta>::cvt_to_seq(aa_seq, aa_tri_seq)) 
            continue;

        // init parser
        BeamCKYParser<ScoreType, IndexType> parser(lambda, is_verbose);

        auto protein = util::split(aa_tri_seq, ' ');
        // parse
        auto system_start = chrono::system_clock::now();
        auto dfa = get_dfa<IndexType>(aa_graphs_with_ln_weights, util::split(aa_tri_seq, ' '));
        auto result = parser.parse(dfa, codon, aa_seq, protein, aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit, aa_graphs_with_ln_weights);
        auto system_diff = chrono::system_clock::now() - system_start;
        auto system_duration = chrono::duration<double>(system_diff).count();  

        // output
        output_result(result, system_duration, lambda, is_verbose, codon, CODON_TABLE);

#ifdef FINAL_CHECK
        if (codon.cvt_rna_seq_to_aa_seq(result.sequence) != aa_seq) {
            std::cerr << "Final Check Failed:" << std::endl;
            std::cerr << codon.cvt_rna_seq_to_aa_seq(result.sequence) << std::endl;
            std::cerr << aa_seq << std::endl;
            assert(false);
        }
#endif
    }
    return 0;
}

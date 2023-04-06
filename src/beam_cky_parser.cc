
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <chrono>
#include <string>

#include "beam_cky_parser.h"
#include "Utils/utility_v.h"
#include "backtrace_iter.cc"
#include "Utils/common.h"

using namespace std;

using NodeType = std::pair<LinearDesign::IndexType, LinearDesign::NumType>;

#define tetra_hex_tri -1

namespace LinearDesign {

template <typename ScoreType, typename IndexType, typename NodeType>
double BeamCKYParser<ScoreType, IndexType, NodeType>::get_broken_codon_score(
        const NodeType& start_node, const NodeType& end_node) {

    IndexType s_index = start_node.first;
    IndexType t_index = end_node.first;

    if (s_index >= t_index)
        return 0.0;

    auto aa_left = protein[s_index / 3]; // tri letter

    auto aa_right = protein[(int)(s_index / 3)];
    if (t_index / 3 < protein.size()){
        aa_right = protein[(int)(t_index / 3)];
    }
        
    auto start_node_re_index = make_pair(s_index % 3, start_node.second);
    auto end_node_re_index = make_pair(t_index % 3, end_node.second);

    double ret = 0.0;

    if (t_index - s_index < 3) {
        if (s_index / 3 == t_index / 3) {
            ret = std::get<0>(best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index,end_node_re_index)]);
        }else{
            double left_ln_cai = 0.0, right_ln_cai = 0.0;
            if (s_index % 3 != 0) 
                left_ln_cai = std::get<0>(best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index,make_pair(0, 0))]);
            if (t_index % 3 != 0) 
                right_ln_cai = std::get<0>(best_path_in_one_codon_unit[aa_right][make_tuple(make_pair(0, 0), end_node_re_index)]);
            ret = left_ln_cai + right_ln_cai;
        }
    }else{
        double left_ln_cai = 0.0, right_ln_cai = 0.0;
        if (s_index % 3 != 0) 
            left_ln_cai = std::get<0>(best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index,make_pair(0, 0))]);
        if (t_index % 3 != 0) 
            right_ln_cai = std::get<0>(best_path_in_one_codon_unit[aa_right][make_tuple(make_pair(0, 0), end_node_re_index)]);
        ret = left_ln_cai + right_ln_cai;
    }
    return ret;
}

template <typename ScoreType, typename IndexType, typename NodeType>
template <IndexType j_num>
void BeamCKYParser<ScoreType, IndexType, NodeType>::hairpin_beam(IndexType j, DFA_t& dfa) {
 
    auto j_node = make_pair(j,j_num);

    for (auto &j1_node_nucj : dfa.right_edges[j_node]) { // right_edges[j][j_num][j1_num][nuc]: false/true
        
        auto j1_node = std::get<0>(j1_node_nucj);
        auto nucj = std::get<1>(j1_node_nucj);

        auto weight_nucj = std::get<2>(j1_node_nucj);


        for (auto &j4_node : dfa.nodes[j+4]){
            const auto& jnext_list = next_pair[nucj][j4_node];

            if (jnext_list.empty())
                continue;

            for (auto &jnext_node_nucjnext : jnext_list){
                auto jnext_node = std::get<0>(jnext_node_nucjnext);
                auto nucjnext = std::get<1>(jnext_node_nucjnext);
                auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);
                auto jnext = jnext_node.first;

                IndexType hairpin_length = jnext + 1 - j; //special hairpin
                NodeNucpair temp = {j, j_num, static_cast<NucPairType>(NTP(nucj, nucjnext))};


#ifdef SPECIAL_HP
                if (hairpin_length == 5 or hairpin_length == 6 or hairpin_length == 8){
                    for(auto & seq_score_weight : hairpin_seq_score_cai[j_node][jnext_node][NTP(nucj, nucjnext)]){
                            auto seq =  get<0>(seq_score_weight);
                            auto pre_cal_score =  get<1>(seq_score_weight);
                            auto pre_cal_cai_score = get<2>(seq_score_weight);
                            update_if_better(bestH[jnext_node][temp], pre_cal_score, pre_cal_cai_score);
                    }

                    continue;
                }
#endif

                for (auto &j2_node_nucj1 : dfa.right_edges[j1_node]) {
                    auto j2_node = std::get<0>(j2_node_nucj1);
                    auto j2_num = j2_node.second;
                    auto nucj1 = std::get<1>(j2_node_nucj1);
                    auto weight_nucj1 = std::get<2>(j2_node_nucj1);

                    for (auto& jnext_1_node_list : dfa.auxiliary_left_edges[jnext_node]){

                        NodeType jnext_1_node = jnext_1_node_list.first;
                        NumType jnext_1_num = jnext_1_node.second;
                        if (jnext - j == 4 and (jnext_1_num != j2_num and dfa.nodes[j+2].size() == dfa.nodes[jnext-1].size())) continue;

                        for (auto& nucjnext_1_weight : jnext_1_node_list.second){

                            IndexType nucjnext_1 = get<0>(nucjnext_1_weight);
                            auto weight_nucjnext_1 = get<1>(nucjnext_1_weight);

                            auto newscore = - func12(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                            
                            FinalScoreType cai_score =  weight_nucj + weight_nucj1 + weight_nucjnext_1 + weight_nucjnext; //ZL need to add weight_nucjnext

                            if ((jnext_1_node.first - j2_node.first) <= SINGLE_MAX_LEN)
                                cai_score += get_broken_codon_score_map[j2_node][jnext_1_node];
                            else
                                cai_score += get_broken_codon_score(j2_node,jnext_1_node);

                            update_if_better(bestH[jnext_node][temp], newscore, cai_score); 

                        }
                    }
                }
            }
        }
    }

    // for every state h in H[j]
    //    1. extend h(i, j) to h(i, jnext)
    //    2. generate p(i, j)
    for (size_t i_node_nucpair_ = 0; i_node_nucpair_ < 16 * j; ++i_node_nucpair_) {

        if (bestH[j_node][i_node_nucpair_].score == util::value_min<ScoreType>())
            continue;

        auto i_node_nucpair = reverse_index(i_node_nucpair_);
        
        auto i = i_node_nucpair.node_first;
        auto i_num = i_node_nucpair.node_second;
        auto pair_nuc = i_node_nucpair.nucpair;
        auto i_node = make_pair(i,i_num);

        auto nuci = PTLN(pair_nuc);
        auto nucj = PTRN(pair_nuc);


        for (const auto& item : dfa.auxiliary_right_edges[j_node]){
            auto j1_node = item.first;
            auto jnext_list = next_pair[nuci][j1_node];

            if (jnext_list.empty()) continue;

            for (auto &jnext_node_nucjnext : jnext_list){
                auto jnext_node = std::get<0>(jnext_node_nucjnext);
                auto nucjnext = std::get<1>(jnext_node_nucjnext);
                auto jnext = jnext_node.first;
                auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);
                auto hairpin_length = jnext + 1 - i;

                NodeNucpair temp = {i, i_num, static_cast<NucPairType>(NTP(nuci, nucjnext))};

#ifdef SPECIAL_HP

                if (hairpin_length == 5 or hairpin_length == 6 or hairpin_length == 8){
                    for(auto & seq_score_weight : hairpin_seq_score_cai[i_node][jnext_node][NTP(nuci, nucjnext)]){
                            auto seq =  get<0>(seq_score_weight);
                            auto pre_cal_score =  get<1>(seq_score_weight);
                            auto pre_cal_cai_score = get<2>(seq_score_weight);
                            update_if_better(bestH[jnext_node][temp], pre_cal_score, pre_cal_cai_score);
                    }

                    continue;
                }

#endif          

                for (auto &i1_node_newnuci : dfa.right_edges[i_node]){
                    NucType newnuci = get<1>(i1_node_newnuci);
                    if (nuci != newnuci) continue;
                    NodeType i1_node = get<0>(i1_node_newnuci);
                    double weight_newnuci = get<2>(i1_node_newnuci);


                    for (auto &i2_node_nuci1 : dfa.right_edges[i1_node]) {
                        auto i2_node = get<0>(i2_node_nuci1);
                        auto nuci1 = get<1>(i2_node_nuci1);
                        auto weight_nuci1 = get<2>(i2_node_nuci1);

                        for (auto &jnext_1_node_nucjnext_1 : dfa.left_edges[jnext_node]) {
                            auto jnext_1_node = get<0>(jnext_1_node_nucjnext_1);
                            auto nucjnext_1 = get<1>(jnext_1_node_nucjnext_1);
                            auto weight_nucjnext_1 = get<2>(jnext_1_node_nucjnext_1);

                            auto newscore = - func12(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);

                            FinalScoreType cai_score =  weight_newnuci + weight_nuci1 + weight_nucjnext_1 + weight_nucjnext; //move weight_nucjnext from H to P to here. Since we added SH here, so it must be here.

                            if ((jnext_1_node.first - i2_node.first) <= SINGLE_MAX_LEN)
                                cai_score += get_broken_codon_score_map[i2_node][jnext_1_node];
                            else
                                cai_score += get_broken_codon_score(i2_node,jnext_1_node);

                            update_if_better(bestH[jnext_node][temp], newscore, cai_score);

                        }
                    }
                }
            }
        }

        auto& state = bestH[j_node][i_node_nucpair_];

        for (auto &j1_node_newnucj : dfa.right_edges[j_node]){
            NucType newnucj = get<1>(j1_node_newnucj);
            if (nucj != newnucj) continue;
            NodeType j1_node = get<0>(j1_node_newnucj);
            update_if_better(bestP[j1_node][i_node_nucpair_], state.score, state.cai_score);
        }
    }

}

template <typename ScoreType, typename IndexType, typename NodeType>
template <IndexType j_num>
void BeamCKYParser<ScoreType, IndexType, NodeType>::Multi_beam(IndexType j, DFA_t& dfa){
    
    NodeType j_node = make_pair(j, j_num);

    for (size_t i_node_nucpair_ = 0; i_node_nucpair_ < 16 * j; ++i_node_nucpair_){

        auto& new_state_score = bestMulti[j_node][i_node_nucpair_];

        if (new_state_score.score == util::value_min<ScoreType>())
            continue;

        auto i_node_nucpair = reverse_index(i_node_nucpair_);
        auto i = i_node_nucpair.node_first;
        auto i_num = i_node_nucpair.node_second;
        auto pair_nuc = i_node_nucpair.nucpair;
        auto nuci = PTLN(pair_nuc);
        auto nucj_1 = PTRN(pair_nuc);

        auto& jnext_list = next_pair[nuci][j_node];

        if (!jnext_list.empty()){

            for (auto &jnext_node_nucjnext : jnext_list){
                auto jnext_node = std::get<0>(jnext_node_nucjnext);
                auto nucjnext = std::get<1>(jnext_node_nucjnext);
                auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);
                auto jnext = jnext_node.first;

                for (auto &jnext1_node_newnucjnext : dfa.right_edges[jnext_node]){
                    auto jnext1_node = std::get<0>(jnext1_node_newnucjnext);
                    auto newnucjnext = std::get<1>(jnext1_node_newnucjnext);
                    if (newnucjnext == nucjnext){
                        double cai_score;

                        if ((jnext_node.first - new_state_score.pre_node.first) <= SINGLE_MAX_LEN)
                            cai_score = new_state_score.pre_left_cai + (get_broken_codon_score_map[new_state_score.pre_node][jnext_node] + weight_nucjnext);
                        else    
                            cai_score = new_state_score.pre_left_cai + (get_broken_codon_score(new_state_score.pre_node, jnext_node) + weight_nucjnext);

                        NodeNucpair temp = {i, i_num, static_cast<NucPairType>(NTP(nuci, nucjnext))};

                        update_if_better(bestMulti[jnext1_node][temp], new_state_score.score, cai_score, new_state_score.pre_node, new_state_score.pre_left_cai);
                    }
                }
            }
        }
        //  2. generate multi(i, j) -> p(i, j)
        auto newscore = new_state_score.score - func15(i, j, nuci, -1, -1, nucj_1, seq_length); // hzhang: TODO
        update_if_better(bestP[j_node][i_node_nucpair_], newscore, new_state_score.cai_score);
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
template <IndexType j_num>
void BeamCKYParser<ScoreType, IndexType, NodeType>::P_beam(IndexType j, DFA_t& dfa){
    
    auto j_node = make_pair(j, j_num);

    if (j < seq_length){
        for (size_t i_node_nucpair_ = 0; i_node_nucpair_ < 16 * j; ++i_node_nucpair_){
        
            auto& state = bestP[j_node][i_node_nucpair_];
            if (state.score == util::value_min<ScoreType>())
                continue;

            auto i_node_nucpair = reverse_index(i_node_nucpair_);
            auto i = i_node_nucpair.node_first;
            
            if (i <= 0) continue;

            auto i_num = i_node_nucpair.node_second;
            auto pair_nuc = i_node_nucpair.nucpair;
            auto i_node = make_pair(i, i_num);

            auto nuci = PTLN(pair_nuc);
            auto nucj_1 = PTRN(pair_nuc);  

            // stacking
            for (auto &j1_node_nucj : dfa.right_edges[j_node]){
                auto j1_node = std::get<0>(j1_node_nucj);
                auto nucj = std::get<1>(j1_node_nucj);
                auto weight_nucj = std::get<2>(j1_node_nucj);

                for (auto &i_1_node_nuci_1 : dfa.left_edges[i_node]){
                    auto i_1_node = std::get<0>(i_1_node_nuci_1);
                    auto nuci_1 = std::get<1>(i_1_node_nuci_1);
                    auto weight_nuci_1 = std::get<2>(i_1_node_nuci_1);
                    auto outer_pair = NTP(nuci_1, nucj);
                    if (_allowed_pairs[nuci_1][nucj]){
                        auto newscore = stacking_score[outer_pair-1][pair_nuc-1] + state.score;
                        double cai_score = state.cai_score + (weight_nuci_1 + weight_nucj);
                        NodeNucpair temp = {i_1_node.first, i_1_node.second, static_cast<NucPairType>(NTP(nuci_1, nucj))};
                        update_if_better(bestP[j1_node][temp], newscore, cai_score);
                    }
                }
            }

            // right bulge: ((...)..) 
            for (auto &j1_node_list : dfa.auxiliary_right_edges[j_node]){
                auto j1_node = j1_node_list.first;

                for (auto &i_1_node_nuci_1 : dfa.left_edges[i_node]){
                    auto i_1_node = std::get<0>(i_1_node_nuci_1);
                    auto nuci_1 = std::get<1>(i_1_node_nuci_1);
                    auto weight_nuci_1 = std::get<2>(i_1_node_nuci_1);

                    auto q_list = next_list[nuci_1][j1_node]; 

                    for (auto& q_node_nucq : q_list){

                        auto q_node = std::get<0>(q_node_nucq);

                        auto q_num = q_node.second;
                        auto q = q_node.first;

                        if (q-j > SINGLE_MAX_LEN) break;
                        
                        auto nucq = std::get<1>(q_node_nucq);
                        auto weight_nucq = std::get<2>(q_node_nucq);
                        auto outer_pair = NTP(nuci_1, nucq);
                        
                        for(auto& q1_node_list : dfa.auxiliary_right_edges[q_node]){
                            NodeType q1_node = q1_node_list.first;
                            if(dfa.nodes[q].size() == 1 and dfa.nodes[q+1].size() == 2 and ((q1_node_list.second)[0]).first != nucq) continue;
                            
                            auto newscore = bulge_score[outer_pair-1][pair_nuc-1][q-j-1]
                                            + state.score;

                            double cai_score;
                            if ((q_node.first - j_node.first) <= SINGLE_MAX_LEN)
                                cai_score = state.cai_score + (weight_nuci_1 + get_broken_codon_score_map[j_node][q_node] + weight_nucq);
                            else
                                cai_score = state.cai_score + (weight_nuci_1 + get_broken_codon_score(j_node, q_node) + weight_nucq);

                            NodeNucpair temp = {i_1_node.first, i_1_node.second, static_cast<NucPairType>(outer_pair)};

                            update_if_better(bestP[q1_node][temp], newscore, cai_score);
                            break;
                        }
                    }                        
                }
            }

            // left bulge: (..(...)) 
            for (auto &j1_node_nucj : dfa.right_edges[j_node]){
                auto j1_node = std::get<0>(j1_node_nucj);
                auto nucj = std::get<1>(j1_node_nucj);
                auto weight_nucj = std::get<2>(j1_node_nucj);

                for (auto &i_1_node_list : dfa.auxiliary_left_edges[i_node]){
                    auto i_1_node = i_1_node_list.first;
                    auto p_list = prev_list[nucj][i_1_node]; 

                    for (auto &p_node_nucp_1 : p_list){

                        auto p_node = std::get<0>(p_node_nucp_1);
                        auto p_num = p_node.second;
                        auto p = p_node.first;

                        if (i-p > SINGLE_MAX_LEN) break;

                        auto nucp_1 = std::get<1>(p_node_nucp_1);
                        auto outer_pair = NTP(nucp_1, nucj);

                        for(auto& p_1_node_new_nucp_1 : dfa.left_edges[p_node]){
                            
                            NucType new_nucp_1 = std::get<1>(p_1_node_new_nucp_1);
                            if(nucp_1 != new_nucp_1) continue;
                            NodeType p_1_node = std::get<0>(p_1_node_new_nucp_1);
                            auto weight_nucp_1 = std::get<2>(p_1_node_new_nucp_1);

                            auto newscore = bulge_score[outer_pair-1][pair_nuc-1][i-p-1]
                                            + state.score;
                            
                            double cai_score;

                            if ((i_node.first - p_node.first) <= SINGLE_MAX_LEN)
                                cai_score = state.cai_score + (weight_nucp_1 + get_broken_codon_score_map[p_node][i_node] + weight_nucj);
                            else
                                cai_score = state.cai_score + (weight_nucp_1 + get_broken_codon_score(p_node, i_node) + weight_nucj);

                            NodeNucpair temp = {p_1_node.first, (NumType)p_1_node.second, static_cast<NucPairType>(outer_pair)};
                            update_if_better(bestP[j1_node][temp], newscore, cai_score);
                        }
                    }                 
                }
            }

            // internal loop
            for (auto &j1_node_dict : dfa.auxiliary_right_edges[j_node]){
                auto j1_node = j1_node_dict.first;
                auto j1_num = j1_node.second;

                for (auto &i_1_node_nuci_1 : dfa.left_edges[i_node]){
                    auto i_1_node = std::get<0>(i_1_node_nuci_1);
                    auto i_1_num = i_1_node.second;
                    auto nuci_1 = std::get<1>(i_1_node_nuci_1);
                    auto weight_nuci_1 = std::get<2>(i_1_node_nuci_1);

                    for (IndexType p = i-1; p > max(i - SINGLE_MAX_LEN, 0); --p) {//ZL, i-(p-1)<=len => i - len < p
                        vector<pair<int, NumType>> p_node_list;
                        
                        if (p == i - 1)
                            p_node_list.push_back(i_1_node);
                        else if (p == i - 2) // hzhang: N.B. add this p, i-1, i o--o--o
                            for (auto &p_node_dict : dfa.auxiliary_left_edges[i_1_node])
                                p_node_list.push_back(p_node_dict.first);
                        else
                            p_node_list = dfa.nodes[p];
                        
                        for (auto &p_node : p_node_list){
                            for (auto &p1_node_nucp : dfa.right_edges[p_node]){

                                auto p1_node = std::get<0>(p1_node_nucp);
                                auto p1_num = p1_node.second;
                                auto nucp = std::get<1>(p1_node_nucp);
                                auto weight_nucp = std::get<2>(p1_node_nucp);

                                if (p == i - 1 and nucp != nuci_1) continue;
                                else if (p == i - 2 and p1_num != i_1_num) continue;
                                else if (p == i - 3 and p1_num != i_1_num and dfa.nodes[p+1].size() == dfa.nodes[i-1].size()) continue; 

                                for (auto &p_1_node_nucp_1 : dfa.left_edges[p_node]){
                                    auto p_1_node = std::get<0>(p_1_node_nucp_1);
                                    auto nucp_1 = std::get<1>(p_1_node_nucp_1);
                                    auto weight_nucp_1 = std::get<2>(p_1_node_nucp_1);

                                    auto q_list = next_list[nucp_1][j1_node]; 

                                    for (auto &q_node_nucq : q_list){

                                        auto q_node = std::get<0>(q_node_nucq);
                                        auto q_num = q_node.second;
                                        auto q = q_node.first;

                                        if (i-p+q-j > SINGLE_MAX_LEN) //check if q is still in the internal loop limit boundary.
                                            break;

                                        auto nucq = std::get<1>(q_node_nucq);
                                        auto weight_nucq = std::get<2>(q_node_nucq);

                                        for(auto& q1_node_list : dfa.auxiliary_right_edges[q_node]){
                                            NodeType q1_node = q1_node_list.first;
                                            if(dfa.nodes[q].size() == 1 and dfa.nodes[q+1].size() == 2 and ((q1_node_list.second)[0]).first != nucq) continue;
                                            NodeNucpair temp = {p_1_node.first, p_1_node.second, static_cast<NucPairType>(NTP(nucp_1, nucq))};
                                            auto& BestP_val = bestP[q1_node][temp];

                                            for(auto & nucj_weightj: j1_node_dict.second){
                                                auto nucj = nucj_weightj.first;
                                                auto weight_nucj = nucj_weightj.second;
                                                if (q == j+1){
                                                    auto newscore = - func14(p-1, q, i, j-1, nucp_1, nucp, nucj, nucq, nuci_1, nuci, nucj_1, nucj) + state.score;
                                                    
                                                    double weight_left;
                                                    if (p == i-1){
                                                        weight_left = weight_nucp_1 + weight_nucp;
                                                    }
                                                    else{
                                                        if (i_1_node.first - p1_node.first <= SINGLE_MAX_LEN)
                                                            weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score_map[p1_node][i_1_node] + weight_nuci_1;
                                                        else
                                                            weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score(p1_node, i_1_node) + weight_nuci_1;
                                                    }
                                                    double cai_score = state.cai_score + (weight_left + weight_nucj + weight_nucq); //j+1 == q

                                                    update_if_better(BestP_val, newscore, cai_score);
                                                }else if (q == j+2){
                                                    for(auto& q_1_node_list : dfa.auxiliary_left_edges[q_node]){
                                                        auto q_1_node = q_1_node_list.first;
                                                        NumType q_1_num = q_1_node.second;
                                                        if (q_1_num != j1_num) continue;
                                                        for(auto & nucq_1_weight : q_1_node_list.second){
                                                            auto nucq_1 = nucq_1_weight.first;
                                                            auto weight_nucq_1 = nucq_1_weight.second;
                                                            auto newscore = - func14(p-1, q, i, j-1, nucp_1, nucp, nucq_1, nucq, nuci_1, nuci, nucj_1, nucj) + state.score;
                                                            
                                                            double weight_left;
                                                            if (p == i-1){
                                                                weight_left = weight_nucp_1 + weight_nucp;
                                                            }
                                                            else{
                                                                // assert(p < i-1);
                                                                if (i_1_node.first - p1_node.first <= SINGLE_MAX_LEN)
                                                                    weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score_map[p1_node][i_1_node] + weight_nuci_1;
                                                                else
                                                                    weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score(p1_node, i_1_node) + weight_nuci_1;
                                                            }

                                                            auto cai_score = state.cai_score + (weight_left + weight_nucj + weight_nucq_1 + weight_nucq);

                                                            update_if_better(BestP_val, newscore, cai_score);
                                                        }
                                                        if(dfa.nodes[q-1].size() == 2) break;
                                                    }
                                                }else if (q == j + 3){
                                                    for(auto& q_1_node_list : dfa.auxiliary_left_edges[q_node]){
                                                        auto q_1_node = q_1_node_list.first;
                                                        NumType q_1_num = q_1_node.second;
                                                        if (q_1_num != j1_num and dfa.nodes[q-1].size() == dfa.nodes[j+1].size()) continue;
                                                        for(auto & nucq_1_weight : q_1_node_list.second){
                                                            auto nucq_1 = nucq_1_weight.first;
                                                            auto weight_nucq_1 = nucq_1_weight.second;
                                                            auto newscore = - func14(p-1, q, i, j-1, nucp_1, nucp, nucq_1, nucq, nuci_1, nuci, nucj_1, nucj) + state.score;
                                                            
                                                            double weight_left;
                                                            if (p == i-1){
                                                                weight_left = weight_nucp_1 + weight_nucp;
                                                            }
                                                            else{
                                                                // assert(p < i-1);
                                                                if (i_1_node.first - p1_node.first <= SINGLE_MAX_LEN)
                                                                    weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score_map[p1_node][i_1_node] + weight_nuci_1;
                                                                else
                                                                    weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score(p1_node, i_1_node) + weight_nuci_1;
                                                            }

                                                            double cai_score;
                                                            if (q_1_node.first - j1_node.first <= SINGLE_MAX_LEN)
                                                                cai_score = state.cai_score + (weight_left + weight_nucj + get_broken_codon_score_map[j1_node][q_1_node] + weight_nucq_1 + weight_nucq);
                                                            else
                                                                cai_score = state.cai_score + (weight_left + weight_nucj + get_broken_codon_score(j1_node, q_1_node) + weight_nucq_1 + weight_nucq);

                                                            update_if_better(BestP_val, newscore, cai_score);
                                                        }
                                                        if(dfa.nodes[q-1].size() == 2) break;
                                                    }
                                                }else{
                                                    for(auto& q_1_node_list : dfa.auxiliary_left_edges[q_node]){
                                                        auto q_1_node = q_1_node_list.first;
                                                        for(auto & nucq_1_weight : q_1_node_list.second){
                                                            auto nucq_1 = nucq_1_weight.first;
                                                            auto weight_nucq_1 = nucq_1_weight.second;
                                                            auto newscore = - func14(p-1, q, i, j-1, nucp_1, nucp, nucq_1, nucq, nuci_1, nuci, nucj_1, nucj) + state.score;


                                                            double weight_left;
                                                            if (p == i-1){
                                                                weight_left = weight_nucp_1 + weight_nucp;
                                                            }
                                                            else{
                                                                // assert(p < i-1);
                                                                if (i_1_node.first - p1_node.first <= SINGLE_MAX_LEN){
                                                                    weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score_map[p1_node][i_1_node] + weight_nuci_1;
                                                                }
                                                                else
                                                                    weight_left = weight_nucp_1 + weight_nucp + get_broken_codon_score(p1_node, i_1_node) + weight_nuci_1;
                                                            }

                                                            double cai_score;
                                                            if (q_1_node.first - j1_node.first <= SINGLE_MAX_LEN){
                                                                cai_score = state.cai_score + (weight_left + weight_nucj + get_broken_codon_score_map[j1_node][q_1_node] + weight_nucq_1 + weight_nucq);
                                                            }
                                                            else
                                                                cai_score = state.cai_score + (weight_left + weight_nucj + get_broken_codon_score(j1_node, q_1_node) + weight_nucq_1 + weight_nucq);
                                                            
                                                            update_if_better(BestP_val, newscore, cai_score);
                                                        }
                                                        if(dfa.nodes[q-1].size() == 2) break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // M = P and M_P = P
    for (size_t i_node_nucpair_ = 0; i_node_nucpair_ < 16 * j; ++i_node_nucpair_){
        auto& state = bestP[j_node][i_node_nucpair_];
        if (state.score == util::value_min<ScoreType>())
            continue;

        auto i_node_nucpair = reverse_index(i_node_nucpair_);
        auto i = i_node_nucpair.node_first;
        auto i_num = i_node_nucpair.node_second;
        auto pair_nuc = i_node_nucpair.nucpair;
        auto i_node = make_pair(i, i_num);

        auto nuci = PTLN(pair_nuc);
        auto nucj_1 = PTRN(pair_nuc);

        if (i > 0 and j < seq_length){

            auto M1_score = - func6(i, j-1, j-1, -1, nuci, nucj_1, -1, seq_length) + state.score;
            
            update_if_better(bestM[j_node][i_node], M1_score, state.cai_score);
            update_if_better(bestM_P[j_node][i_node], M1_score, state.cai_score);
        }
    }

    // M2 = M + M_P
    for (size_t i_node_ = 0; i_node_ < 2 * j; ++i_node_) {
        auto& state = bestM_P[j_node][i_node_];
        auto i_node = reverse_index2(i_node_);
        auto i = i_node.first;

        if (state.score == util::value_min<ScoreType>())
            continue;

        if (i > 0 and j < seq_length){
            for (size_t m_node = 0; m_node < 2 * i; ++m_node){
                auto& m_new_state_score = bestM[i_node][m_node];

                if (m_new_state_score.score == util::value_min<ScoreType>())
                    continue;

                auto newscore = m_new_state_score.score + state.score;
                auto cai_score = m_new_state_score.cai_score + state.cai_score;
                update_if_better(bestM2[j_node][m_node], newscore, cai_score);
            }
        }
    }

    // C = C + P
    for (size_t i_node_nucpair_ = 0; i_node_nucpair_ < 16 * j; ++i_node_nucpair_){
        auto& state = bestP[j_node][i_node_nucpair_];
        if (state.score == util::value_min<ScoreType>())
            continue;

        auto i_node_nucpair = reverse_index(i_node_nucpair_);
        auto i = i_node_nucpair.node_first;
        auto i_num = i_node_nucpair.node_second;
        auto pair_nuc = i_node_nucpair.nucpair;
        auto i_node = make_pair(i, i_num);

        auto nuci = PTLN(pair_nuc);
        auto nucj_1 = PTRN(pair_nuc); 
        
        if (i > 0){
            auto& prefix_C = bestC[i_node];

            if (prefix_C.score != util::value_min<ScoreType>()){
                auto newscore = - func3(i, j-1, nuci, nucj_1, seq_length) + prefix_C.score + state.score;

                auto cai_score = prefix_C.cai_score + state.cai_score;
                update_if_better(bestC[j_node], newscore, cai_score);
            }
        }
        else{
            auto newscore = - func3(0, j-1, nuci, nucj_1, seq_length) + state.score;

            update_if_better(bestC[j_node], newscore, state.cai_score);
        }
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
template <IndexType j_num>
void BeamCKYParser<ScoreType, IndexType, NodeType>::M2_beam(IndexType j, DFA_t& dfa){
    
    auto j_node = make_pair(j, j_num);
    for (size_t i_node_ = 0; i_node_ < 2 * j; ++i_node_) {
        auto& state = bestM2[j_node][i_node_];
        if (state.score == util::value_min<ScoreType>())
            continue;

        auto i_node = reverse_index2(i_node_);
        auto i = i_node.first;

        // 1. multi-loop
        for (IndexType p = i-1; p >= max(i - SINGLE_MAX_LEN, 0); --p){
            vector<pair<int, NumType>> p_node_list;
            if (p == i - 1)
                for(auto& p_node_dict : dfa.auxiliary_left_edges[i_node])
                    p_node_list.push_back(p_node_dict.first);
            else p_node_list = dfa.nodes[p];

            for (auto &p_node : p_node_list){
                for (auto &p1_node_nucp : dfa.right_edges[p_node]){
                    auto p1_node = std::get<0>(p1_node_nucp);
                    auto nucp = std::get<1>(p1_node_nucp);
                    auto weight_nucp = std::get<2>(p1_node_nucp);

                    if(p == i - 1 and p1_node != i_node) continue;
                    if(p == i - 2 and dfa.nodes[p+1].size() == dfa.nodes[i].size() and p1_node.second != i_node.second) continue;

                    auto q_list = next_pair[nucp][j_node];

                    for (auto &q_node_nucq : q_list){
                        auto q_node = std::get<0>(q_node_nucq);
                        auto nucq = std::get<1>(q_node_nucq);
                        auto weight_nucq = std::get<2>(q_node_nucq);
                        auto q = q_node.first;

                        if (i - p + q - j - 1 > SINGLE_MAX_LEN) continue; //ZL, i-p-1+q-j
                        auto outer_pair = NTP(nucp, nucq);
                        for (auto &q1_node_newnucq : dfa.right_edges[q_node]){
                            auto newnucq = std::get<1>(q1_node_newnucq);
                            if (newnucq == nucq) {
                                auto q1_node = std::get<0>(q1_node_newnucq);

                                double cai_score = state.cai_score + (weight_nucp + get_broken_codon_score_map[p1_node][i_node] + get_broken_codon_score_map[j_node][q_node] + weight_nucq);
                                double temp_left_cai = state.cai_score + (weight_nucp + get_broken_codon_score_map[p1_node][i_node]);
                                NodeNucpair temp = {p_node.first, p_node.second, static_cast<NucPairType>(NTP(nucp, nucq))};
                                update_if_better(bestMulti[q1_node][temp], state.score, cai_score, j_node, temp_left_cai);
                                break;
                            }
                        }
                    }
                }
            }
        }
        //  2. M = M2
        update_if_better(bestM[j_node][i_node], state.score, state.cai_score);
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
template <IndexType j_num>
void BeamCKYParser<ScoreType, IndexType, NodeType>::M_beam(IndexType j, DFA_t& dfa)
{
    
    auto j_node = make_pair(j, j_num);

    for (size_t i_node_ = 0; i_node_ < 2 * j; ++i_node_) {
        auto& state = bestM[j_node][i_node_];

        if (state.score == util::value_min<ScoreType>())
            continue;

        auto i_node = reverse_index2(i_node_);
        for (auto &j1_node_nucj : dfa.right_edges[j_node]){
            auto j1_node = std::get<0>(j1_node_nucj);
            auto nucj = std::get<1>(j1_node_nucj);
            auto weight_nucj = std::get<2>(j1_node_nucj);

            double cai_score = state.cai_score + weight_nucj;
            update_if_better(bestM[j1_node][i_node], state.score, cai_score);
        }
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
template <IndexType j_num>
void BeamCKYParser<ScoreType, IndexType, NodeType>::C_beam(IndexType j, DFA_t& dfa)
{
    //  beam of C
    //  C = C + U
    auto j_node = make_pair(j, j_num);

    auto& state = bestC[j_node];


    for (auto &j1_node_nucj : dfa.right_edges[j_node]){
        NodeType j1_node = std::get<0>(j1_node_nucj);
        IndexType nucj = std::get<1>(j1_node_nucj);
        auto weight_nucj = std::get<2>(j1_node_nucj);

        double cai_score = state.cai_score + (double)weight_nucj;
        update_if_better(bestC[j1_node], state.score, cai_score);
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
void BeamCKYParser<ScoreType, IndexType, NodeType>::get_next_pair(DFA_t& dfa) {
    vector<tuple<NodeType, NucType, double>> temp_vector; 
    for (NucType nuci = 0; nuci < NOTON; nuci++) {
        for (IndexType j = seq_length; j > 0; j--) {

            for (auto& j_node : dfa.nodes[j]) {
                for (auto& item : dfa.auxiliary_left_edges[j_node]) {
                    NodeType j_1_node = item.first;
                    temp_vector.clear();
                    for (auto& nuc_weight : item.second){
                        auto nuc = std::get<0>(nuc_weight);
                        auto weight_nuc = std::get<1>(nuc_weight);
                        if (_allowed_pairs[nuci][nuc])
                            temp_vector.push_back(make_tuple(j_1_node, nuc, weight_nuc));
                    }
                    if(temp_vector.size() == 0){
                        if (next_pair[nuci][j_1_node].size() > 0 and next_pair[nuci][j_node].size() > 0) {
                            // merge
                            IndexType index1 = std::get<0>(next_pair[nuci][j_1_node][0]).first;
                            IndexType index2 = std::get<0>(next_pair[nuci][j_node][0]).first;
                            if(index1/3 == index2/3)
                                next_pair[nuci][j_1_node].insert(next_pair[nuci][j_1_node].end(),
                                                                     next_pair[nuci][j_node].begin(),
                                                                     next_pair[nuci][j_node].end());
                            else if(index1 > index2){
                                next_pair[nuci][j_1_node].clear();
                                next_pair[nuci][j_1_node].insert(next_pair[nuci][j_1_node].end(),
                                                                     next_pair[nuci][j_node].begin(),
                                                                     next_pair[nuci][j_node].end());
                            }
                        }else if (next_pair[nuci][j_node].size() > 0)
                            next_pair[nuci][j_1_node].insert(next_pair[nuci][j_1_node].end(),
                                                                     next_pair[nuci][j_node].begin(),
                                                                     next_pair[nuci][j_node].end());
                    }
                    else
                        next_pair[nuci][j_1_node].insert(next_pair[nuci][j_1_node].end(),
                                                                     temp_vector.begin(),
                                                                     temp_vector.end());
                }
            }
        }
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
void BeamCKYParser<ScoreType, IndexType, NodeType>::get_next_pair_set() {

    for(NucType nuci=0; nuci<5; nuci++){
        for (auto& j_node_vnuc : next_pair[nuci]) {
            NodeType j_node = j_node_vnuc.first;
            next_pair_set[nuci][j_node] = set<tuple<NodeType, NucType, double>>(j_node_vnuc.second.begin(), j_node_vnuc.second.end());
        }
    }
    for(NucType nuci=0; nuci<5; nuci++){
        for (auto& j_node_vnuc : next_pair_set[nuci]) {
            NodeType j_node = j_node_vnuc.first;
            next_pair[nuci][j_node].clear();
            for(auto& item : next_pair_set[nuci][j_node]){
                next_pair[nuci][j_node].push_back(item);
            }
        }
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
void BeamCKYParser<ScoreType, IndexType, NodeType>::get_prev_pair(DFA_t& dfa) {
    vector<tuple<NodeType, NucType, double>> temp_vector;
    for (NucType nuci = 0; nuci < NOTON; nuci++) {
        for (IndexType j = 0; j < seq_length; j++) {
            for (auto& j_node : dfa.nodes[j]) {
                for (auto& item : dfa.auxiliary_right_edges[j_node]) {
                    NodeType j1_node = item.first;
                    temp_vector.clear();
                    for (auto& nuc_weight : item.second){
                        auto nuc = std::get<0>(nuc_weight);
                        auto weight_nuc = std::get<1>(nuc_weight);    
                        if (_allowed_pairs[nuci][nuc])
                            temp_vector.push_back(make_tuple(j1_node, nuc, weight_nuc));
                    }
                    if(temp_vector.size() == 0){
                        if (prev_pair[nuci][j1_node].size() > 0 and prev_pair[nuci][j_node].size() > 0) {
                            // merge
                            IndexType index1 = std::get<0>(prev_pair[nuci][j1_node][0]).first-1;
                            IndexType index2 = std::get<0>(prev_pair[nuci][j_node][0]).first-1;
                            if(index1/3 == index2/3)
                                prev_pair[nuci][j1_node].insert(prev_pair[nuci][j1_node].end(),
                                                                     prev_pair[nuci][j_node].begin(),
                                                                     prev_pair[nuci][j_node].end());
                            else if(index1 < index2){
                                prev_pair[nuci][j1_node].clear();
                                prev_pair[nuci][j1_node].insert(prev_pair[nuci][j1_node].end(),
                                                                     prev_pair[nuci][j_node].begin(),
                                                                     prev_pair[nuci][j_node].end());
                            }
                        }else if (prev_pair[nuci][j_node].size() > 0)
                            prev_pair[nuci][j1_node].insert(prev_pair[nuci][j1_node].end(),
                                                                     prev_pair[nuci][j_node].begin(),
                                                                     prev_pair[nuci][j_node].end());
                    }
                    else
                        prev_pair[nuci][j1_node].insert(prev_pair[nuci][j1_node].end(),
                                                                     temp_vector.begin(),
                                                                     temp_vector.end());
                }
            }
        }
    }
}

template <typename ScoreType, typename IndexType, typename NodeType>
void BeamCKYParser<ScoreType, IndexType, NodeType>::get_prev_pair_set() {

    for(NucType nuci=0; nuci<5; nuci++){
        for (auto& j_node_vnuc : prev_pair[nuci]) {
            NodeType j_node = j_node_vnuc.first;
            prev_pair_set[nuci][j_node] =
                set<tuple<NodeType, NucType, double>>(j_node_vnuc.second.begin(), j_node_vnuc.second.end());
        }
    }
    for(NucType nuci=0; nuci<5; nuci++){
        for (auto& j_node_vnuc : prev_pair_set[nuci]) {
            NodeType j_node = j_node_vnuc.first;
            prev_pair[nuci][j_node].clear();
            for(auto& item : prev_pair_set[nuci][j_node]){
                prev_pair[nuci][j_node].push_back(item);
            }
        }
    }
}

#ifdef SPECIAL_HP
template <typename ScoreType, typename IndexType, typename NodeType>
void BeamCKYParser<ScoreType, IndexType, NodeType>::special_hp(DFA_t& dfa, int8_t hairpin_length) {
    int8_t hairpin_type = HAIRPINTYPE(hairpin_length);
    vector<tuple<NodeType, string, double, NodeType>> queue;
    vector<tuple<NodeType, string, double, NodeType>> frontier; 
    // vector
    for(IndexType i=0; i<=seq_length - hairpin_length; i++){
        for(NodeType i_node : dfa.nodes[i]){
            int count = hairpin_length;
            queue.clear();
            queue.push_back(make_tuple(i_node, "", double(0.), i_node));
            while(count > 0){
                count --;
                frontier.clear();
                for(auto& node_str : queue){
                    NodeType cur_node = std::get<0>(node_str);
                    string cur_str = std::get<1>(node_str);
                    double cur_lncai = std::get<2>(node_str);
                    for(auto& node_nuc : dfa.right_edges[cur_node]){
                        NodeType new_node = std::get<0>(node_nuc);
                        string new_str = cur_str + GET_ACGU(std::get<1>(node_nuc));
                        double new_total_lncai = cur_lncai + std::get<2>(node_nuc);
                        frontier.push_back(make_tuple(new_node, new_str, new_total_lncai, cur_node));
                    }
                }
                queue.swap(frontier);
            }
            for(auto node_str : queue){
                auto j_node = std::get<3>(node_str);
                auto temp_seq = std::get<1>(node_str);
                auto cai_score = std::get<2>(node_str);
                auto hairpin_length = temp_seq.size();
                int8_t hairpin_type = HAIRPINTYPE(hairpin_length);
                NucType nuci = GET_ACGU_NUC(temp_seq[0]);
                NucType nucj = GET_ACGU_NUC(temp_seq[temp_seq.size() - 1]);
                auto temp_nucpair = NTP(nuci, nucj);

                ScoreType special_hairpin_score = func1(temp_seq, hairpin_type);
                if(special_hairpin_score == SPECIAL_HAIRPIN_SCORE_BASELINE){

                    auto newscore = - func12(0, hairpin_length - 1, GET_ACGU_NUC(temp_seq[0]), GET_ACGU_NUC(temp_seq[1]), GET_ACGU_NUC(temp_seq[hairpin_length-2]), GET_ACGU_NUC(temp_seq[hairpin_length-1]), tetra_hex_tri);
                    hairpin_seq_score_cai[i_node][j_node][temp_nucpair].push_back(make_tuple(temp_seq, newscore, cai_score));
                }

                else{
                    hairpin_seq_score_cai[i_node][j_node][temp_nucpair].push_back(make_tuple(temp_seq, special_hairpin_score, cai_score));
                }
            }
        }
    }
}
#endif

template <typename ScoreType, typename IndexType, typename NodeType>
void BeamCKYParser<ScoreType, IndexType, NodeType>::preprocess(DFA_t& dfa) {

    vector<tuple<NodeType, NucType, double>> new_q_list, new_p_list;
    set<NodeType> visited; 

    // next_list
    NodeType init_node = make_pair(0, 0);
    for (NucType nuci=1; nuci<NOTON; nuci++) {
        visited.clear();
        auto q_list = next_pair[nuci][init_node]; 

        while (!q_list.empty()){
            new_q_list.clear();

            for (auto& q_node_nucq : q_list){

                auto q_node = std::get<0>(q_node_nucq);
                auto q_num = q_node.second;
                auto q = q_node.first; 
                auto nucq = std::get<1>(q_node_nucq);

                // q_node
                next_list[nuci][q_node].push_back(q_node_nucq);
                // q-1 is special
                for(auto& q_1_node_dict : dfa.auxiliary_left_edges[q_node]){
                    NodeType q_1_node = q_1_node_dict.first;
                    next_list[nuci][q_1_node].push_back(q_node_nucq);
                }
                for(IndexType j=q-2; j>=max(0, q-SINGLE_MAX_LEN-1); j--)
                    for(NodeType j_node : dfa.nodes[j])
                        next_list[nuci][j_node].push_back(q_node_nucq);

                for(auto& q1_node_list : dfa.auxiliary_right_edges[q_node]){
                    NodeType q1_node = q1_node_list.first;

                    if(dfa.nodes[q].size() == 1 and dfa.nodes[q+1].size() == 2 and ((q1_node_list.second)[0]).first != nucq) continue;
                
                    if(visited.find(q1_node) == visited.end()){
                        visited.insert(q1_node);
                        new_q_list.insert(new_q_list.end(), next_pair[nuci][q1_node].cbegin(), next_pair[nuci][q1_node].cend());
                    }
                    break;
                }
            }
            q_list.swap(new_q_list);
        }
    }

    // prev_list
    init_node = make_pair(seq_length, 0);
    for (NucType nucj=1; nucj<NOTON; nucj++) {
        visited.clear();
        auto p_list = prev_pair[nucj][init_node]; 

        while (!p_list.empty()){
            new_p_list.clear();

            for (auto& p_node_nucp_1 : p_list){
                auto p_node = std::get<0>(p_node_nucp_1);
                auto p_num = p_node.second;
                auto p = p_node.first; 
                auto nucp_1 = std::get<1>(p_node_nucp_1); 

                // p_node
                prev_list[nucj][p_node].push_back(p_node_nucp_1);
                // p+1 is special
                for(auto& p1_node_dict : dfa.auxiliary_right_edges[p_node]){
                    NodeType p1_node = p1_node_dict.first;
                    prev_list[nucj][p1_node].push_back(p_node_nucp_1);
                }
                for(IndexType i=p+2; i<=min(seq_length, p+SINGLE_MAX_LEN+1); i++)
                    for(NodeType i_node : dfa.nodes[i])
                        prev_list[nucj][i_node].push_back(p_node_nucp_1);

                for(auto& p_1_node_new_nucp_1 : dfa.left_edges[p_node]){
                                    
                    NucType new_nucp_1 = std::get<1>(p_1_node_new_nucp_1);
                    if(nucp_1 != new_nucp_1) continue;
                    NodeType p_1_node = std::get<0>(p_1_node_new_nucp_1);
                
                    if(visited.find(p_1_node) == visited.end()){
                        visited.insert(p_1_node);
                        new_p_list.insert(new_p_list.end(), prev_pair[nucj][p_1_node].cbegin(), prev_pair[nucj][p_1_node].cend());
                    }
                }
            }
            p_list.swap(new_p_list);
        }
    }

    // stacking energy computation
    int newscore;
    for(int8_t outer_pair=1; outer_pair<=6; outer_pair++){
        auto nuci_1 = PTLN(outer_pair);
        auto nucq = PTRN(outer_pair);
        for(int8_t inner_pair=1; inner_pair<=6; inner_pair++){
            auto nuci = PTLN(inner_pair);
            auto nucj_1 = PTRN(inner_pair);
            newscore = - func14(0, 1, 1, 0,
                                nuci_1, nuci, nucj_1, nucq,
                                nuci_1, nuci, nucj_1, nucq);
            stacking_score[outer_pair-1][inner_pair-1] = newscore;

            for (IndexType l=0; l<=SINGLE_MAX_LEN; l++){
                newscore = - func14(0, l+2, 1, 0,
                              nuci_1, nuci, nucj_1, nucq,
                              nuci_1, nuci, nucj_1, nucq); 

                bulge_score[outer_pair-1][inner_pair-1][l] = newscore;
            }
        }   
    }

#ifdef SPECIAL_HP
    // Triloops
    special_hp(dfa, 5);
    // Tetraloop37
    special_hp(dfa, 6);
    // Hexaloops
    special_hp(dfa, 8);
#endif
}

template <typename ScoreType, typename IndexType, typename NodeType>
DecoderResult<double, IndexType> BeamCKYParser<ScoreType, IndexType, NodeType>::parse(
        DFA_t& dfa, Codon& codon, std::string& aa_seq, std::vector<std::string>& p, 
        std::unordered_map<std::string, std::string>& aa_best_in_codon,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, 
            std::hash<std::tuple<NodeType, NodeType>>>>& best_path_in_one_codon,
            std::unordered_map<string, Lattice<IndexType>>& aa_graphs_with_ln_weights) {
    

    protein = p;
    aa_graphs_with_ln_w = aa_graphs_with_ln_weights;
    aa_best_path_in_a_whole_codon = aa_best_in_codon;
    best_path_in_one_codon_unit = best_path_in_one_codon;

    seq_length = 3 * static_cast<IndexType>(aa_seq.size());
    next_pair.resize(5);
    next_pair_set.resize(5);
    get_next_pair(dfa);
    get_next_pair_set();

    prev_pair.resize(5);
    prev_pair_set.resize(5);
    get_prev_pair(dfa);
    get_prev_pair_set();

    next_list.resize(5);
    prev_list.resize(5);
    stacking_score.resize(6, vector<int>(6));
    bulge_score.resize(6, vector<vector<int>>(6, vector<int>(SINGLE_MAX_LEN+1)));

    preprocess(dfa);


    int reserved_size = (seq_length + 1) * 16; //node,nucpair
    int reserved_size2 = (seq_length + 1) * 2; //node

    bestH.resize(reserved_size2);
    bestP.resize(reserved_size2);
    bestM2.resize(reserved_size2); //slim signature, Liang Zhang
    bestMulti.resize(reserved_size2);
    bestM.resize(reserved_size2); //slim signature, Liang Zhang
    bestM_P.resize(reserved_size2); // hzhang: inter-state: P -> M

    bestC.resize(reserved_size2); //slim signature, Liang Zhang

    get_broken_codon_score_map.resize(reserved_size2);
    for (auto& e : get_broken_codon_score_map){ //slim signature, Liang Zhang
        e.resize(reserved_size2);
        for (auto& ee : e){
            ee = util::value_min<double>();
        }
    }

    for (auto& ee : bestC){
        ee.score = util::value_min<ScoreType>();
        ee.cai_score = util::value_min<double>();
    }


    for (auto& e : bestH){
        e.resize(reserved_size);
        for (auto& ee : e){
            ee.score = util::value_min<ScoreType>();
            ee.cai_score = util::value_min<double>();
        }

    }

    for (auto& e : bestP){
        e.resize(reserved_size);
        for (auto& ee : e){
            ee.score = util::value_min<ScoreType>();
            ee.cai_score = util::value_min<double>();
        }
    }

    for (auto& e : bestMulti){ //slim signature, Liang Zhang
        e.resize(reserved_size);
        for (auto& ee : e){
            ee.score = util::value_min<ScoreType>();
            ee.cai_score = util::value_min<double>();
        }        
    }

    for (auto& e : bestM2){ //slim signature, Liang Zhang
        e.resize(reserved_size2);
        for (auto& ee : e){
            ee.score = util::value_min<ScoreType>();
            ee.cai_score = util::value_min<double>();
        }        
    }

    for (auto& e : bestM){ //slim signature, Liang Zhang
        e.resize(reserved_size2);
        for (auto& ee : e){
            ee.score = util::value_min<ScoreType>();
            ee.cai_score = util::value_min<double>();
        }
    }

    for (auto& e : bestM_P){ // hzhang
        e.resize(reserved_size2);
        for (auto& ee : e){
            ee.score = util::value_min<ScoreType>();
            ee.cai_score = util::value_min<double>();
        }
    }

    for (IndexType i = 0; i <= seq_length; ++i) {
        for (auto & node_i : dfa.nodes[i]){
            for (IndexType l = 0; l <= SINGLE_MAX_LEN; ++l){
                auto j = i + l;
                if (j > seq_length)
                    break;

                for (auto & node_j : dfa.nodes[j]){
                    get_broken_codon_score_map[node_i][node_j] = get_broken_codon_score(node_i, node_j);
                }

            }
        }

    }

    bestC[make_pair(0,0)].score = 0;
    bestC[make_pair(0,0)].cai_score = double(0.);

    for(const auto& node_nue_weight : dfa.right_edges[make_pair(0,0)]) {
        auto node = std::get<0>(node_nue_weight);
        auto weight_nue = std::get<2>(node_nue_weight);
        update_if_better(bestC[node], 0, weight_nue);
    }
    
    for (IndexType j = 0; j <= seq_length; ++j) {
        cout << "j=" << j << "\r" << flush;
        
        hairpin_beam<0>(j, dfa);
        hairpin_beam<1>(j, dfa);
        
        if (j == 0)
            continue;
        
        Multi_beam<0>(j, dfa);
        Multi_beam<1>(j, dfa);
        P_beam<0>(j, dfa);
        P_beam<1>(j, dfa);
        M2_beam<0>(j, dfa);
        M2_beam<1>(j, dfa);
        
        if (j < seq_length) {
            M_beam<0>(j, dfa);
            M_beam<1>(j, dfa);
            C_beam<0>(j, dfa);
            C_beam<1>(j, dfa);
        }
    }

    auto end_node = make_pair(seq_length, 0);
    auto viterbi = bestC[end_node];

    auto backtrace_result = backtrace(dfa, viterbi, end_node);

    return DecoderResult<double, IndexType>{backtrace_result.seq, backtrace_result.structure, viterbi.score / -100.0, 0., viterbi.cai_score, seq_length};
}

template <typename ScoreType, typename IndexType, typename NodeType>
BeamCKYParser<ScoreType, IndexType, NodeType>::BeamCKYParser(const double lambda_value, const bool verbose)
        : lambda(lambda_value), is_verbose(verbose) {
        func9(0, 0);
}

}
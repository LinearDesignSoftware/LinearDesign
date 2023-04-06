
#include "beam_cky_parser.h"

using namespace std;

#define tetra_hex_tri -1

namespace LinearDesign {

template <typename ScoreType, typename IndexType, typename NodeType>
BacktraceResult BeamCKYParser<ScoreType, IndexType, NodeType>::backtrace(DFA_t& dfa, const State_t& state, NodeType end_node){

    char sequence[seq_length+1];
    memset(sequence, '.', seq_length);
    sequence[seq_length] = 0;

    char structure[seq_length+1];
    memset(structure, '.', seq_length);
    structure[seq_length] = 0;

    bool no_backpointer;

    stack<tuple<NodeType, NodeType, State_t, Beam_type, PairType>> stk;
    NodeType start_node = make_pair(0, 0);
    stk.push(make_tuple(start_node, end_node, state, Beam_type::BEAM_C, -1));

    double epsilon = 1e-8;

    while(!stk.empty()) {
        tuple<NodeType, NodeType, State_t, Beam_type, PairType> top = stk.top();
        NodeType i_node = get<0>(top), j_node = get<1>(top);
        State_t& state = get<2>(top);
        Beam_type beam_type = get<3>(top);
        PairType curr_pair_nuc = get<4>(top);
        stk.pop();

        IndexType i, j, p, q, hairpin_length;
        j = j_node.first;
        NucType nuci, nucj, nuci1, nucj_1;
        no_backpointer = true;
        
        int left_start, left_end, right_start, right_end;

        switch (beam_type) {
            case Beam_type::BEAM_C:
                if (j <= 0) continue;
                for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                    auto j_1_node = std::get<0>(j_1_node_nucj_1);
                    auto& c_state = bestC[j_1_node];
                    auto weight_nucj_1 = std::get<2>(j_1_node_nucj_1);
                    auto cai_score = c_state.cai_score + weight_nucj_1;

                    if (state.score == c_state.score && abs(state.cai_score - cai_score) < epsilon){
                        NucType nucj_1 = std::get<1>(j_1_node_nucj_1);
                        stk.push(make_tuple(i_node, j_1_node, c_state, Beam_type::BEAM_C, curr_pair_nuc));
                        sequence[j-1] = GET_ACGU(nucj_1);
                        no_backpointer = false;
                        break;
                    }
                }

                // C = C + P
                if(no_backpointer) {
                    for (size_t c_node_nucpair_ = 0; c_node_nucpair_ < 16 * seq_length; ++c_node_nucpair_){
                        auto& p_state = bestP[j_node][c_node_nucpair_];

                        if (p_state.score == util::value_min<ScoreType>()) continue;
                        auto c_node_nucpair = reverse_index(c_node_nucpair_);

                        auto c = c_node_nucpair.node_first;
                        auto c_num = c_node_nucpair.node_second;
                        auto pair_nuc = c_node_nucpair.nucpair;
                        auto c_node = make_pair(c, c_num);

                        auto nucc = PTLN(pair_nuc);
                        auto nucj_1 = PTRN(pair_nuc);                         


                        auto newscore = - func3(c, j-1, nucc, nucj_1, seq_length) + p_state.score;

                        if (c > 0){
                            auto& c_state = bestC[c_node];
                            auto cai_score = c_state.cai_score + p_state.cai_score;

                            if (state.score == c_state.score + newscore && abs(state.cai_score - cai_score) < epsilon){
                                stk.push(make_tuple(i_node, c_node, c_state, Beam_type::BEAM_C, curr_pair_nuc));
                                stk.push(make_tuple(c_node, j_node, p_state, Beam_type::BEAM_P, pair_nuc));
                                no_backpointer = false;
                                break;
                            }
                        } else{
                            if (state.score == newscore && abs(state.cai_score - p_state.cai_score) < epsilon){
                                stk.push(make_tuple(c_node, j_node, p_state, Beam_type::BEAM_P, pair_nuc));
                                no_backpointer = false;
                                break;
                            }
                        }
                    }
                    if (!no_backpointer) break;
                }
                assert(no_backpointer == false); // something wrong if no path matches
                break;

            case Beam_type::BEAM_P:
                i = i_node.first;
                j = j_node.first;
                nuci = PTLN(curr_pair_nuc);
                nucj_1 = PTRN(curr_pair_nuc);


                hairpin_length = j - i;

                for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                    NucType new_nucj_1 = std::get<1>(j_1_node_nucj_1);
                    if (new_nucj_1 != nucj_1) continue;

                    auto j_1_node = std::get<0>(j_1_node_nucj_1);
                    auto weight_nucj_1 = std::get<2>(j_1_node_nucj_1);



#ifdef SPECIAL_HP
                    if (hairpin_length == 5 or hairpin_length == 6 or hairpin_length == 8){
                        for(auto & seq_score_weight : hairpin_seq_score_cai[i_node][j_1_node][NTP(nuci, nucj_1)]){
                            auto seq =  get<0>(seq_score_weight);
                            auto pre_cal_score =  get<1>(seq_score_weight);
                            auto pre_cal_cai_score = get<2>(seq_score_weight);


                            if (state.score == pre_cal_score && abs(state.cai_score - pre_cal_cai_score) < epsilon){
                                for(int c=0; c<seq.size(); c++)
                                    sequence[i+c] = seq[c];

                                structure[i] = '(';
                                structure[j-1] = ')';

                                no_backpointer = false;
                                break;
                            }
                        }if (!no_backpointer) break;
                    }if (!no_backpointer) break;
#endif



                    for (auto& i1_node_nuci : dfa.right_edges[i_node]){
                        NucType new_nuci = std::get<1>(i1_node_nuci);
                        if (new_nuci != nuci) continue;
                        auto i1_node = std::get<0>(i1_node_nuci);
                        auto weight_nuci = std::get<2>(i1_node_nuci);

                        // helix 
                        for (auto& j_2_node_nucj_2 : dfa.left_edges[j_1_node]){
                            NucType nucj_2 = std::get<1>(j_2_node_nucj_2);

                            for (auto& i2_node_nuci1 : dfa.right_edges[i1_node]){
                                NucType nuci1 = std::get<1>(i2_node_nuci1);
                                auto pair_nuc = NTP(nuci1, nucj_2);

                                NodeNucpair temp = {i1_node.first, i1_node.second, static_cast<NucPairType>(pair_nuc)};

                                auto& p_state = bestP[j_1_node][temp];
                                auto newscore = - func14(i, j-1, i+1, j-2,
                                                                    nuci, nuci1, nucj_2, nucj_1,
                                                                    nuci, nuci1, nucj_2, nucj_1)
                                                        + p_state.score;

                                auto cai_score = p_state.cai_score + (weight_nuci + weight_nucj_1);

                                if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                    stk.push(make_tuple(i1_node, j_1_node, p_state, Beam_type::BEAM_P, pair_nuc));
                                    sequence[i] = GET_ACGU(nuci);
                                    sequence[j-1] = GET_ACGU(nucj_1);
                                    structure[i] = '(';
                                    structure[j-1] = ')';

                                    no_backpointer = false;
                                    break;
                                } 
                            }if (!no_backpointer) break;
                        }if (!no_backpointer) break;
                    }
                    if (!no_backpointer) break;

                    // hairpin
                    NodeNucpair temp = {i_node.first, LinearDesign::NumType(i_node.second), curr_pair_nuc};

                    if (state.score == bestH[j_1_node][temp].score){ //no need to check CAI score here

                        for (auto& j_2_node_nucj_2 : dfa.left_edges[j_1_node]){
                            NucType nucj_2 = std::get<1>(j_2_node_nucj_2);
                            auto j_2_node = std::get<0>(j_2_node_nucj_2);
                            auto j_2 = j_2_node.first;

                            auto weight_nucj_2 = std::get<2>(j_2_node_nucj_2);

                            for (auto& i1_node_nuci : dfa.right_edges[i_node]){
                                NucType new_nuci = std::get<1>(i1_node_nuci);
                                if (new_nuci != nuci) continue;

                                auto i1_node = std::get<0>(i1_node_nuci);
                                auto weight_nuci = get<2>(i1_node_nuci);
                                

                                for(auto& i2_node_nuci1 : dfa.right_edges[i1_node]){
                                    NucType nuci1 = std::get<1>(i2_node_nuci1);
                                    auto i2_node = std::get<0>(i2_node_nuci1);
                                    auto i2 = i2_node.first;
                                    auto weight_nuci1 = std::get<2>(i2_node_nuci1);

                                    if (j - 1 - i == 4 and (j_2_node.second != i2_node.second and dfa.nodes[i+2].size() == dfa.nodes[j-2].size())) continue;

                                    auto newscore = - func12(i, j-1, nuci, nuci1, nucj_2, nucj_1, tetra_hex_tri);
                                    auto cai_score = weight_nuci + weight_nuci1 + get_broken_codon_score(i2_node, j_2_node) + weight_nucj_2 + weight_nucj_1;

                                    if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                        sequence[i] = GET_ACGU(nuci);
                                        sequence[i+1] = GET_ACGU(nuci1);
                                        sequence[j-2] = GET_ACGU(nucj_2);
                                        sequence[j-1] = GET_ACGU(nucj_1);
                                        structure[i] = '(';
                                        structure[j-1] = ')';

                                        auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, i2_node, j_2_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                        int count = i2;
                                        for (auto & nuc : temp_string ){
                                            sequence[count] = nuc;
                                            count++;
                                        }
                                        assert(count == j_2);

                                        no_backpointer = false;
                                        break;
                                    }
                                }if (!no_backpointer) break;
                            }if (!no_backpointer) break;
                        }if (!no_backpointer) break;
                    }
                }

                // single branch
                if (no_backpointer) {
                    vector<pair<IndexType, NucType>> right_seq;
                    vector<tuple<NodeType, NodeType, NucType, NucType, NucType, vector<pair<IndexType, NucType>>, int, int, NodeType, NodeType, double, double, double, double, bool, NodeType>> q_node_nucs_list;
                    for (IndexType q = j-1; q >= std::max(j - SINGLE_MAX_LEN - 1, i + 5); --q){    
                        int right_start = -1;
                        int right_end = -1;
                        q_node_nucs_list.clear();

                        if (q == j-1){
                            for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                                if (get<1>(j_1_node_nucj_1) != nucj_1) continue;
                                auto q_node = get<0>(j_1_node_nucj_1);
                                auto weight_nucj_1 = get<2>(j_1_node_nucj_1);
                                for (auto& q_1_node_nucq_1 : dfa.left_edges[q_node]){
                                    NodeType q_1_node = get<0>(q_1_node_nucq_1);
                                    auto nucq_1 = get<1>(q_1_node_nucq_1);
                                    double weight_nucq_1 = get<2>(q_1_node_nucq_1);
                                    right_seq.push_back(make_pair(j-1, nucj_1));
                                    q_node_nucs_list.push_back(make_tuple(q_1_node, q_node, nucq_1, nucj_1, nucq_1, right_seq, right_start, right_end, make_pair(-1,0), make_pair(-1,0), weight_nucq_1, 0., 0., weight_nucj_1, true, make_pair(-1,0)));
                                    right_seq.clear();
                                }
                            }
                        }else if(q == j-2){
                            for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                                if (get<1>(j_1_node_nucj_1) != nucj_1) continue;
                                auto j_1_node = get<0>(j_1_node_nucj_1);
                                auto weight_nucj_1 = get<2>(j_1_node_nucj_1);
                                for (auto& q_node_nucq : dfa.left_edges[j_1_node]){
                                    auto q_node = get<0>(q_node_nucq);
                                    auto nucq = get<1>(q_node_nucq);
                                    auto weight_nucq = get<2>(q_node_nucq);
                                    for (auto& q_1_node_nucq_1 : dfa.left_edges[q_node]){
                                        NodeType q_1_node = get<0>(q_1_node_nucq_1);
                                        auto nucq_1 = get<1>(q_1_node_nucq_1);
                                        double weight_nucq_1 = get<2>(q_1_node_nucq_1);
                                        right_seq.push_back(make_pair(q, nucq));
                                        right_seq.push_back(make_pair(j-1, nucj_1));
                                        q_node_nucs_list.push_back(make_tuple(q_1_node, q_node, nucq_1, nucq, nucq, right_seq, right_start, right_end, make_pair(-1,0), make_pair(-1,0), weight_nucq_1, weight_nucq, 0., weight_nucj_1, false, j_1_node));
                                        right_seq.clear();
                                    }
                                }
                            }
                        }else if(q == j-3){
                            for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                                if (get<1>(j_1_node_nucj_1) != nucj_1) continue;
                                auto j_1_node = get<0>(j_1_node_nucj_1);
                                auto weight_nucj_1 = get<2>(j_1_node_nucj_1);
                                for(auto& j_2_node_nucj_2 : dfa.left_edges[j_1_node]){
                                    auto j_2_node = get<0>(j_2_node_nucj_2);
                                    auto nucj_2 = get<1>(j_2_node_nucj_2);
                                    auto weight_nucj_2 = get<2>(j_2_node_nucj_2);
                                    for (auto& q_node_nucq : dfa.left_edges[j_2_node]){
                                        auto q_node = get<0>(q_node_nucq);
                                        auto nucq = get<1>(q_node_nucq);
                                        auto weight_nucq = get<2>(q_node_nucq);
                                        for (auto& q_1_node_nucq_1 : dfa.left_edges[q_node]){
                                            NodeType q_1_node = get<0>(q_1_node_nucq_1);
                                            auto nucq_1 = get<1>(q_1_node_nucq_1);
                                            double weight_nucq_1 = get<2>(q_1_node_nucq_1);
                                            right_seq.push_back(make_pair(q, nucq));
                                            right_seq.push_back(make_pair(j-2, nucj_2));
                                            right_seq.push_back(make_pair(j-1, nucj_1));
                                            q_node_nucs_list.push_back(make_tuple(q_1_node, q_node, nucq_1, nucq, nucj_2, right_seq, right_start, right_end, make_pair(-1,0), make_pair(-1,0), weight_nucq_1, weight_nucq, weight_nucj_2, weight_nucj_1, false, j_1_node));
                                            right_seq.clear();
                                        }
                                    }
                                }
                            }
                        }
                        else if(q == j-4){
                            for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                                if (get<1>(j_1_node_nucj_1) != nucj_1) continue;
                                auto j_1_node = get<0>(j_1_node_nucj_1);
                                auto weight_nucj_1 = get<2>(j_1_node_nucj_1);
                                for(auto& j_2_node_nucj_2 : dfa.left_edges[j_1_node]){
                                    auto j_2_node = get<0>(j_2_node_nucj_2);
                                    auto nucj_2 = get<1>(j_2_node_nucj_2);
                                    auto weight_nucj_2 = get<2>(j_2_node_nucj_2);
                                    for(auto& j_3_node_nucj_3 : dfa.left_edges[j_2_node]){
                                        auto j_3_node = get<0>(j_3_node_nucj_3);
                                        for (auto& q_node_nucq : dfa.left_edges[j_3_node]){
                                            auto q_node = get<0>(q_node_nucq);
                                            auto nucq = get<1>(q_node_nucq);
                                            auto weight_nucq = get<2>(q_node_nucq);
                                            for (auto& q_1_node_nucq_1 : dfa.left_edges[q_node]){
                                                NodeType q_1_node = get<0>(q_1_node_nucq_1);
                                                auto nucq_1 = get<1>(q_1_node_nucq_1);
                                                double weight_nucq_1 = get<2>(q_1_node_nucq_1);
                                                right_seq.push_back(make_pair(q, nucq));
                                                right_seq.push_back(make_pair(j-2, nucj_2));
                                                right_seq.push_back(make_pair(j-1, nucj_1));
                                                right_start = q+1;
                                                right_end = j-2;
                                                q_node_nucs_list.push_back(make_tuple(q_1_node, q_node, nucq_1, nucq, nucj_2, right_seq, right_start, right_end, j_3_node, j_2_node, weight_nucq_1, weight_nucq, weight_nucj_2, weight_nucj_1, false, j_1_node));
                                                right_seq.clear();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else{
                            for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                                if (get<1>(j_1_node_nucj_1) != nucj_1) continue;
                                auto j_1_node = get<0>(j_1_node_nucj_1);
                                auto weight_nucj_1 = get<2>(j_1_node_nucj_1);
                                for(auto& j_2_node_nucj_2 : dfa.left_edges[j_1_node]){
                                    auto j_2_node = get<0>(j_2_node_nucj_2);
                                    auto nucj_2 = get<1>(j_2_node_nucj_2);
                                    auto weight_nucj_2 = get<2>(j_2_node_nucj_2);
                                    for (auto& q_node : dfa.nodes[q]){
                                        for (auto& q1_node_nucq : dfa.right_edges[q_node]){
                                            auto q1_node = get<0>(q1_node_nucq);
                                            auto nucq = get<1>(q1_node_nucq);
                                            auto weight_nucq = get<2>(q1_node_nucq);
                                            for (auto& q_1_node_nucq_1 : dfa.left_edges[q_node]){
                                                NodeType q_1_node = get<0>(q_1_node_nucq_1);
                                                auto nucq_1 = get<1>(q_1_node_nucq_1);
                                                double weight_nucq_1 = get<2>(q_1_node_nucq_1);
                                                right_seq.push_back(make_pair(q, nucq));
                                                right_seq.push_back(make_pair(j-2, nucj_2));
                                                right_seq.push_back(make_pair(j-1, nucj_1));
                                                right_start = q+1;
                                                right_end = j-2;
                                                q_node_nucs_list.push_back(make_tuple(q_1_node, q_node, nucq_1, nucq, nucj_2, right_seq, right_start, right_end, q1_node, j_2_node, weight_nucq_1, weight_nucq, weight_nucj_2, weight_nucj_1,false, j_1_node));
                                                right_seq.clear();
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        for (auto& q_node_nucs : q_node_nucs_list){
                            auto q_1_node = get<0>(q_node_nucs);
                            auto q_node = get<1>(q_node_nucs);
                            auto nucq_1 = get<2>(q_node_nucs);
                            auto nucq = get<3>(q_node_nucs);
                            auto nucj_2 = get<4>(q_node_nucs);
                            auto right_seq = get<5>(q_node_nucs);
                            auto right_start = get<6>(q_node_nucs);
                            auto right_end = get<7>(q_node_nucs);
                            auto right_start_node = get<8>(q_node_nucs);
                            auto right_end_node = get<9>(q_node_nucs);
                            auto weight_nucq_1 = get<10>(q_node_nucs);
                            auto weight_nucq = get<11>(q_node_nucs);
                            auto weight_nucj_2 = get<12>(q_node_nucs);
                            auto weight_nucj_1 = get<13>(q_node_nucs);
                            bool q_equ_j_1 = get<14>(q_node_nucs);
                            auto j_1_node = get<15>(q_node_nucs);

                            double weight_right = 0.0;
                            double weight_left = 0.0;

                            if(q_equ_j_1){
                                for (auto& i1_node_nuci : dfa.right_edges[i_node]){
                                    NucType new_nuci = get<1>(i1_node_nuci);
                                    if (new_nuci != nuci) continue;
                                    auto i1_node = get<0>(i1_node_nuci);
                                    auto weight_nuci = get<2>(i1_node_nuci);

                                    auto p_list = next_list[nucq_1][i1_node];
                                    for (auto &p_node_nucp : p_list){
                                        auto p_node = get<0>(p_node_nucp);
                                        auto nucp = get<1>(p_node_nucp);
                                        auto p = p_node.first; 
                                        PairType pair_nuc = NTP(nucp, nucq_1);

                                        if (p == i + 1) continue; // stack
                                        if (p - i + j - q - 2 > SINGLE_MAX_LEN) continue;
                                        
                                        NodeNucpair temp = {p_node.first, p_node.second, static_cast<NucPairType>(pair_nuc)};

                                        auto& p_state = bestP[q_node][temp];

                                        auto newscore = - func14(i, j-1, p, q-1, nuci, nucp, nucq_1, nucj_1, nuci, nucp, nucq_1, nucj_1) + p_state.score;
                                        auto weight_left = weight_nuci + get_broken_codon_score(i1_node, p_node);
                                        auto cai_score = p_state.cai_score + (weight_left + weight_nucj_1);

                                        if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                            stk.push(make_tuple(p_node, q_node, p_state, Beam_type::BEAM_P, pair_nuc));

                                            sequence[i] = GET_ACGU(nuci);

                                            auto temp_i1_to_p_nucs = get_nuc_from_dfa_cai<IndexType>(dfa, i1_node, p_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                            assert(temp_i1_to_p_nucs.size() == p - (i+1));
                                            auto count = i+1;
                                            for (auto& nuc : temp_i1_to_p_nucs){
                                                sequence[count] = nuc;
                                                count++;
                                            }
                                            assert(count == p);

                                            assert(right_seq.size() == 1);
                                            sequence[j-1] = GET_ACGU(right_seq[0].second);

                                            structure[i] = '(';
                                            structure[j-1] = ')';

                                            no_backpointer = false;
                                            break;
                                        }if (!no_backpointer) break;
                                    }if (!no_backpointer) break;
                                }
                            }else{
                                for (auto& i1_node_nuci : dfa.right_edges[i_node]){
                                    NucType new_nuci = get<1>(i1_node_nuci);
                                    if (new_nuci != nuci) continue;
                                    auto i1_node = get<0>(i1_node_nuci);
                                    auto weight_nuci = get<2>(i1_node_nuci);

                                    auto p_list = next_list[nucq_1][i1_node];
                                        
                                    for (auto &p_node_nucp : p_list){
                                        auto p_node = get<0>(p_node_nucp);
                                        auto nucp = get<1>(p_node_nucp);
                                        auto weight_nucp = get<2>(p_node_nucp);
                                        auto p = p_node.first; 
                                        PairType pair_nuc = NTP(nucp, nucq_1);

                                        if (p - i + j - q - 2 > SINGLE_MAX_LEN) continue;
                                        
                                        NodeNucpair temp = {p_node.first, p_node.second, static_cast<NucPairType>(pair_nuc)};

                                        auto& p_state = bestP[q_node][temp];
                                        auto newscore = 0;
                                        if (p == i+1){
                                            newscore = - func14(i, j-1, p, q-1, nuci, nucp, nucj_2, nucj_1, nuci, nucp, nucq_1, nucq) 
                                                            + p_state.score;

                                            weight_right = get_broken_codon_score(q_node,j_1_node) + weight_nucj_1;       

                                            auto cai_score = p_state.cai_score + (weight_nuci + weight_right);

                                            if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                                stk.push(make_tuple(p_node, q_node, p_state, Beam_type::BEAM_P, pair_nuc));

                                                sequence[i] = GET_ACGU(nuci);
                                                for(auto& idx_nucidx : right_seq){
                                                    IndexType idx = idx_nucidx.first;
                                                    NucType nucidx = idx_nucidx.second;
                                                    sequence[idx] = GET_ACGU(nucidx);
                                                }

                                                structure[i] = '(';
                                                structure[j-1] = ')';

                                                auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, q_node, j_1_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                                int count = q;
                                                for (auto & nuc : temp_string ){
                                                    sequence[count] = nuc;
                                                    count++;
                                                }

                                                assert(count == j-1);
                                                no_backpointer = false;
                                                break;
                                            }
                                        }else if(p == i+2){
                                            for (auto& i2_node_nuci1 : dfa.right_edges[i1_node]){
                                                auto i2_node = get<0>(i2_node_nuci1);
                                                if (p_node != i2_node) continue;

                                                NucType nuci1 = get<1>(i2_node_nuci1);
                                                auto weight_nuci1 = get<2>(i2_node_nuci1);
                                                newscore = - func14(i, j-1, p, q-1, nuci, nuci1, nucj_2, nucj_1, nuci1, nucp, nucq_1, nucq) 
                                                                + p_state.score;
     
                                                weight_left = weight_nuci + weight_nuci1;
                                                weight_right = weight_nucq + get_broken_codon_score(right_start_node, right_end_node) + weight_nucj_2 + weight_nucj_1;

                                                auto cai_score = p_state.cai_score + (weight_left + weight_right);

                                                if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                                    stk.push(make_tuple(p_node, q_node, p_state, Beam_type::BEAM_P, pair_nuc));

                                                    sequence[i] = GET_ACGU(nuci);
                                                    sequence[i+1] = GET_ACGU(nuci1);
                                                    for(auto& idx_nucidx : right_seq){
                                                        IndexType idx = idx_nucidx.first;
                                                        NucType nucidx = idx_nucidx.second;
                                                        sequence[idx] = GET_ACGU(nucidx);
                                                    }

                                                    structure[i] = '(';
                                                    structure[j-1] = ')';
                                                    auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, right_start_node, right_end_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                                    int count = right_start;
                                                    for (auto & nuc : temp_string ){
                                                        sequence[count] = nuc;
                                                        count++;
                                                    }
                                                    assert(count == right_end);
                                                    no_backpointer = false;
                                                    break;
                                                }        
                                            }  
                                        }else if(p == i+3){
                                            for (auto& p_1_node_nucp_1 : dfa.left_edges[p_node]){
                                                auto p_1_node = get<0>(p_1_node_nucp_1);
                                                auto nucp_1 = get<1>(p_1_node_nucp_1);
                                                auto weight_nucp_1 = get<2>(p_1_node_nucp_1);
                                                for (auto& i2_node_nuci1 : dfa.right_edges[i1_node]){
                                                    auto i2_node = get<0>(i2_node_nuci1);
                                                    if (p_1_node != i2_node) continue;

                                                    NucType nuci1 = get<1>(i2_node_nuci1);
                                                    auto weight_nuci1 = get<2>(i2_node_nuci1);
                                                    newscore = - func14(i, j-1, p, q-1, nuci, nuci1, nucj_2, nucj_1, nucp_1, nucp, nucq_1, nucq) 
                                                            + p_state.score;

                                                    weight_left = weight_nuci + weight_nuci1 + weight_nucp_1;
                                                    weight_right = weight_nucq + get_broken_codon_score(right_start_node,right_end_node) + weight_nucj_2 + weight_nucj_1;

                                                    auto cai_score = p_state.cai_score + (weight_left + weight_right);

                                                    if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                                        stk.push(make_tuple(p_node, q_node, p_state, Beam_type::BEAM_P, pair_nuc));
                                                        
                                                        sequence[i] = GET_ACGU(nuci);
                                                        sequence[i+1] = GET_ACGU(nuci1);
                                                        sequence[p-1] = GET_ACGU(nucp_1);

                                                        for(auto& idx_nucidx : right_seq){
                                                            IndexType idx = idx_nucidx.first;
                                                            NucType nucidx = idx_nucidx.second;
                                                            sequence[idx] = GET_ACGU(nucidx);
                                                        }
                                                        structure[i] = '(';
                                                        structure[j-1] = ')';
                                                        auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, right_start_node, right_end_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                                        int count = right_start;
                                                        for (auto & nuc : temp_string ){
                                                            sequence[count] = nuc;
                                                            count++;
                                                        }
                                                        assert(count == right_end);

                                                        no_backpointer = false;
                                                        break;
                                                    }          
                                                }if (!no_backpointer) break;
                                            }
                                        }else if(p == i+4){
                                            for (auto& p_1_node_nucp_1 : dfa.left_edges[p_node]){
                                                auto p_1_node = get<0>(p_1_node_nucp_1);
                                                auto nucp_1 = get<1>(p_1_node_nucp_1);
                                                auto weight_nucp_1 = get<2>(p_1_node_nucp_1);
                                                for (auto& i2_node_nuci1 : dfa.right_edges[i1_node]){
                                                    auto i2_node = get<0>(i2_node_nuci1);
                                                    NucType nuci1 = get<1>(i2_node_nuci1);
                                                    auto weight_nuci1 = get<2>(i2_node_nuci1);
                                                    for (auto& i3_node_nuci2 : dfa.right_edges[i2_node]){
                                                        auto i3_node = get<0>(i3_node_nuci2);
                                                        if (i3_node != p_1_node) continue;
                                                        auto nuci2 = get<1>(i3_node_nuci2);
                                                        auto weight_nuci2 = get<2>(i3_node_nuci2);
                                                        
                                                        newscore = - func14(i, j-1, p, q-1, nuci, nuci1, nucj_2, nucj_1, nucp_1, nucp, nucq_1, nucq) 
                                                                + p_state.score;

                                                        weight_left = weight_nuci + weight_nuci1 + weight_nuci2 + weight_nucp_1;
                                                        weight_right = weight_nucq + get_broken_codon_score(right_start_node,right_end_node) + weight_nucj_2 + weight_nucj_1;

                                                        auto cai_score = p_state.cai_score + (weight_left + weight_right);

                                                        if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                                            stk.push(make_tuple(p_node, q_node, p_state, Beam_type::BEAM_P, pair_nuc));

                                                            sequence[i] = GET_ACGU(nuci);
                                                            sequence[i+1] = GET_ACGU(nuci1);
                                                            sequence[i+2] = GET_ACGU(nuci2);
                                                            sequence[p-1] = GET_ACGU(nucp_1);

                                                            for(auto& idx_nucidx : right_seq){
                                                                IndexType idx = idx_nucidx.first;
                                                                NucType nucidx = idx_nucidx.second;
                                                                sequence[idx] = GET_ACGU(nucidx);
                                                            }
                                                            structure[i] = '(';
                                                            structure[j-1] = ')';
                                                            auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, right_start_node, right_end_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                                            auto count = right_start;
                                                            for (auto & nuc : temp_string ){
                                                                sequence[count] = nuc;
                                                                count++;
                                                            }
                                                            assert(count == right_end);

                                                            no_backpointer = false;
                                                            break;
                                                        }
                                                    }if (!no_backpointer) break;   
                                                }if (!no_backpointer) break;
                                            }
                                        }else{
                                            for (auto& i2_node_nuci1 : dfa.right_edges[i1_node]){
                                                NucType nuci1 = get<1>(i2_node_nuci1);
                                                auto i2_node = get<0>(i2_node_nuci1);
                                                auto weight_nuci1 = get<2>(i2_node_nuci1);

                                                for (auto& p_1_node_nucp_1 : dfa.left_edges[p_node]){
                                                    auto nucp_1 = get<1>(p_1_node_nucp_1);
                                                    auto p_1_node = get<0>(p_1_node_nucp_1);
                                                    auto weight_nucp_1 = get<2>(p_1_node_nucp_1);

                                                    newscore = - func14(i, j-1, p, q-1, nuci, nuci1, nucj_2, nucj_1, nucp_1, nucp, nucq_1, nucq) 
                                                                    + p_state.score;

                                                    weight_left = weight_nuci + weight_nuci1 + get_broken_codon_score(i2_node, p_1_node) + weight_nucp_1;
                                                    weight_right = weight_nucq + get_broken_codon_score(right_start_node,right_end_node) + weight_nucj_2 + weight_nucj_1;

                                                    auto cai_score = p_state.cai_score + (weight_left + weight_right);
                                                    if (state.score == newscore && abs(state.cai_score - cai_score) < epsilon){
                                                        stk.push(make_tuple(p_node, q_node, p_state, Beam_type::BEAM_P, pair_nuc));

                                                        sequence[i] = GET_ACGU(nuci);
                                                        sequence[i+1] = GET_ACGU(nuci1);
                                                        sequence[p-1] = GET_ACGU(nucp_1);
                                                        auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, i2_node, p_1_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                                        int count = i+2;
                                                        for (auto & nuc : temp_string ){
                                                            sequence[count] = nuc;
                                                            count++;
                                                        }
                                                        assert(count == p-1);

                                                        for(auto& idx_nucidx : right_seq){
                                                            IndexType idx = idx_nucidx.first;
                                                            NucType nucidx = idx_nucidx.second;
                                                            sequence[idx] = GET_ACGU(nucidx);
                                                        }
                                                        structure[i] = '(';
                                                        structure[j-1] = ')';

                                                        temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, right_start_node, right_end_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                                        count = right_start;
                                                        for (auto & nuc : temp_string ){
                                                            sequence[count] = nuc;
                                                            count++;
                                                        }
                                                        assert(count == right_end);
                                                        no_backpointer = false;
                                                        break;
                                                    }   
                                                }if (!no_backpointer) break;
                                            }
                                        }if (!no_backpointer) break;

                                    }if (!no_backpointer) break;
                                }if (!no_backpointer) break;
                            }
                        }if (!no_backpointer) break;
                    }
                }
                
                // Multi
                if (no_backpointer){

                    NodeNucpair temp = {i_node.first, i_node.second, static_cast<NucPairType>(curr_pair_nuc)};


                    auto& multi_state = bestMulti[j_node][temp];
                    auto newscore = multi_state.score - func15(i, j, nuci, -1, -1, nucj_1, seq_length);

                    if (state.score == newscore && abs(state.cai_score - multi_state.cai_score) < epsilon){
                        stk.push(make_tuple(i_node, j_node, multi_state, Beam_type::BEAM_MULTI, curr_pair_nuc));
                        
                        sequence[i] = GET_ACGU(nuci);
                        sequence[j-1] = GET_ACGU(nucj_1);
                        structure[i] = '(';
                        structure[j-1] = ')';

                        no_backpointer = false;
                    }
                }

                assert(no_backpointer == false);
                break;

            case Beam_type::BEAM_MULTI:
                nuci = PTLN(curr_pair_nuc);
                nucj_1 = PTRN(curr_pair_nuc);
                j = j_node.first;
                i = i_node.first;

                
                for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                    auto j_1_node = get<0>(j_1_node_nucj_1);
                    auto weight_nucj_1 = get<2>(j_1_node_nucj_1);
                    NodeType q_node = state.pre_node;
                    q = q_node.first;

                    if(q == j - 1 and q_node != j_1_node) continue;
                    if(q == j - 2 and dfa.nodes[q].size() == dfa.nodes[j-1].size() and q_node.second != j_1_node.second) continue;

                    

                    for (size_t p_node_ = 0; p_node_ < 2 * q; ++p_node_) {

                        auto& temp_state = bestM2[q_node][p_node_];
                        if (temp_state.score == util::value_min<ScoreType>())
                            continue;

                        auto p_node = reverse_index2(p_node_);
                        auto p = p_node.first;

                        if(p <= i) continue;

                        for (auto& i1_node_nuci : dfa.right_edges[i_node]){

                            auto i1_node = get<0>(i1_node_nuci);

                            if(p == i + 1 and p_node != i1_node) continue;
                            if(p == i + 2 and dfa.nodes[p].size() == dfa.nodes[i+1].size() and p_node.second != i1_node.second) continue;


                            double weight_nuci = double(get<2>(i1_node_nuci));

                            auto& m2_state = bestM2[q_node][p_node];

                            auto cai_score = m2_state.cai_score + (weight_nuci + get_broken_codon_score(i1_node, p_node) + get_broken_codon_score(q_node, j_1_node) + weight_nucj_1);

                            if (state.score == m2_state.score && abs(state.cai_score - cai_score) < epsilon){
                                stk.push(make_tuple(p_node, q_node, m2_state, Beam_type::BEAM_M2, -1));

                                auto temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, i1_node, p_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                auto count = i+1;
                                for (auto & nuc : temp_string ){
                                    sequence[count] = nuc;
                                    count++;
                                }
                                assert(count == p);

                                temp_string.clear();
                                temp_string = get_nuc_from_dfa_cai<IndexType>(dfa, q_node, j_1_node, protein, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon);
                                count = q;
                                for (auto & nuc : temp_string ){
                                    sequence[count] = nuc;
                                    count++;
                                }

                                assert(count == j-1);
                                
                                no_backpointer = false;
                            }if (!no_backpointer) break;
                        }if (!no_backpointer) break;
                    }if (!no_backpointer) break;
                }
                assert(no_backpointer == false);
                break;

            case Beam_type::BEAM_M2:
                // M2 = M + P
                i = i_node.first;
                j = j_node.first;

                for (size_t m_node_nucpair_ = 0; m_node_nucpair_ < 16 * j; ++m_node_nucpair_){
                
                    auto& p_state = bestP[j_node][m_node_nucpair_];
                    if (p_state.score == util::value_min<ScoreType>())
                        continue;

                    auto m_node_nucpair = reverse_index(m_node_nucpair_);
                    auto m = m_node_nucpair.node_first;
                    auto m_num = m_node_nucpair.node_second;

                    auto m_node = make_pair(m, m_num);

                    auto pair_nuc = m_node_nucpair.nucpair;

                    if (m <= i+4) continue; // no sharpturn

                    auto nucm = PTLN(pair_nuc);
                    auto nucj_1 = PTRN(pair_nuc); 
                    auto newscore = - func6(-1, -1, -1, -1, nucm, nucj_1, -1, seq_length) + p_state.score;

                    auto& m_state = bestM[m_node][i_node];
                    auto cai_score = m_state.cai_score + p_state.cai_score;

                    if (state.score == m_state.score + newscore && state.cai_score == cai_score){
                        stk.push(make_tuple(i_node, m_node, m_state, Beam_type::BEAM_M1, -1));
                        stk.push(make_tuple(m_node, j_node, p_state, Beam_type::BEAM_P, pair_nuc));

                        no_backpointer = false;
                        break;
                    }
                    if (!no_backpointer) break;
                }
                assert(no_backpointer == false);
                break;

            case Beam_type::BEAM_M1: 
                // M = M + U
                for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                    auto j_1_node = std::get<0>(j_1_node_nucj_1);
                    auto weight_nucj_1 = std::get<2>(j_1_node_nucj_1);
                    auto& m_state = bestM[j_1_node][i_node];
                    auto cai_score = m_state.cai_score + weight_nucj_1;

                    if (state.score == m_state.score && abs(state.cai_score - cai_score) < epsilon) {
                        NucType nucj_1 = std::get<1>(j_1_node_nucj_1);
                        stk.push(make_tuple(i_node, j_1_node, m_state, Beam_type::BEAM_M1, -1));

                        sequence[j-1] = GET_ACGU(nucj_1);

                        no_backpointer = false;
                        break;
                    }
                }

                // M = P
                if(no_backpointer){
                    for (auto& j_1_node_nucj_1 : dfa.left_edges[j_node]){
                        NucType nucj_1 = std::get<1>(j_1_node_nucj_1);

                        for (auto& i1_node_nuci : dfa.right_edges[i_node]){
                            NucType nuci = std::get<1>(i1_node_nuci);
                            PairType pair_nuc = NTP(nuci, nucj_1);
                            
                            NodeNucpair temp = {i_node.first, i_node.second, static_cast<NucPairType>(pair_nuc)};

                            auto& p_state = bestP[j_node][temp];
                            auto newscore = - func6(-1, -1, -1, -1, nuci, nucj_1, -1, seq_length) + p_state.score;

                            if (state.score == newscore && abs(state.cai_score - p_state.cai_score) < epsilon) {
                                stk.push(make_tuple(i_node, j_node, p_state, Beam_type::BEAM_P, pair_nuc));
                                no_backpointer = false;
                                break;
                            }
                        }if(!no_backpointer) break;
                    }
                }

                // M = M2
                if(no_backpointer){
                    auto& m2_state = bestM2[j_node][i_node];



                    if (state.score == m2_state.score && state.cai_score == m2_state.cai_score) {
                        stk.push(make_tuple(i_node, j_node, m2_state, Beam_type::BEAM_M2, -1));
                        no_backpointer = false;
                    }
                }

                assert(no_backpointer == false);
                break;
            default:  // MANNER_NONE or other cases
                
                printf("wrong beam_type at %d, %d\n", i, j); fflush(stdout);
                assert(false);
        }

    }
    assert(string(sequence).size() == string(structure).size());
    return {string(sequence), string(structure)};
}
}


#include <map>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <array>
#include "utility_v.h"
#include "common.h"
#include "codon.h"

using namespace std;

// #define is_verbose

namespace LinearDesign {

template <typename IndexType,
          typename IndexWType = tuple<IndexType, double>,
          typename NodeType = pair<IndexType, NumType>,
          typename NodeNucWType = tuple<NodeType, NucType, double>>
class Lattice {
public:
    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeNucWType>, hash_pair> left_edges;
    unordered_map<NodeType, vector<NodeNucWType>, hash_pair> right_edges;
    
    Lattice(): nodes(), left_edges(), right_edges() {};

    void add_edge(NodeType n1, NodeType n2, NucType nuc, double weight = 0.0f){
        right_edges[n1].push_back(make_tuple(n2, nuc, weight));
        left_edges[n2].push_back(make_tuple(n1, nuc, weight));
    }

    void add_node(NodeType n1){
        IndexType pos = get<0>(n1);
        nodes[pos].push_back(n1);
    }
};

template <typename IndexType,
          typename IndexWType = pair<IndexType, double>,
          typename NodeType = pair<IndexType, NumType>,
          typename NodeNucWType = tuple<NodeType, NucType, double>>
class DFA {
public:
    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeNucWType>, hash_pair> left_edges;
    unordered_map<NodeType, vector<NodeNucWType>, hash_pair> right_edges;
    unordered_map<NodeType, unordered_map<NodeType, vector<IndexWType>, hash_pair>, hash_pair> auxiliary_left_edges;
    unordered_map<NodeType, unordered_map<NodeType, vector<IndexWType>, hash_pair>, hash_pair> auxiliary_right_edges;
    unordered_map<NodeType, unordered_map<IndexType, double>, hash_pair> node_rightedge_weights;

    DFA(): nodes(), left_edges(), right_edges(), auxiliary_left_edges(), auxiliary_right_edges() {};

    void add_edge(NodeType n1, NodeType n2, IndexType nuc, double weight = 0.0f){
        right_edges[n1].push_back(make_tuple(n2, nuc, weight));
        left_edges[n2].push_back(make_tuple(n1, nuc, weight));
        auxiliary_right_edges[n1][n2].push_back(make_pair(nuc, weight));
        auxiliary_left_edges[n2][n1].push_back(make_pair(nuc, weight));
        node_rightedge_weights[n1][nuc] = weight;
    }

    void add_node(NodeType n1){
        IndexType pos = get<0>(n1);
        nodes[pos].push_back(n1);
    }    
};

template <typename IndexType,
          typename NodeType = pair<IndexType, NumType>,
          typename LatticeType = Lattice<IndexType>>
unordered_map<string, LatticeType> read_wheel(string const &filename) {
    unordered_map<string, LatticeType> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        printf("Unable to open coding_wheel file\n");
        exit(1);   // call system to stop
    }

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;

    for (string line; getline(inFile, line);) {
        stuff = util::split(line, '\t');

        aa = stuff[0];
        LatticeType graph = LatticeType();
        graph.add_node(make_pair(0,0)); // always initialize with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = util::split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_pair(2, i);
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_pair(1, i);
                graph.add_node(n1);
                graph.add_edge(make_pair(0, 0), n1, GET_ACGU_NUC(first));
            }
            else {
                n1 = make_pair(1, i-1);
            }
            last_first = first;
            graph.add_edge(n1, n2, GET_ACGU_NUC(second));
            for (auto& third : thirds) {
                graph.add_edge(n2, make_pair(0,0), GET_ACGU_NUC(third));
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
#ifdef is_verbose
        printf("-----------------Lattice------------------------\n");
        for(IndexType pos = 0; pos <= 2; pos++){
            for(auto &node : graph.nodes[pos]){
                IndexType p = get<0>(node);
                IndexType num = get<1>(node);
                printf("node, (%d, %d)\n", p, num);
                for(auto &item : graph.right_edges[node]){
                    NodeType n2 = get<0>(item);
                    IndexType p2 = get<0>(n2); IndexType num2 = get<1>(n2);
                    IndexType nuc = get<1>(item);
                    double weight = get<2>(item);
                    printf("              (%d, %d) -(%d,%lf)-> (%d, %d)\n", p, num, nuc,  weight, p2, num2);
                }
                for(auto &item : graph.left_edges[node]){
                    NodeType n1 = get<0>(item);
                    IndexType p1 = get<0>(n1); IndexType num1 = get<1>(n1);
                    IndexType nuc = get<1>(item);
                    double weight = get<2>(item);
                    printf("  (%d, %d) <-(%d,%lf)- (%d, %d)\n", p1, num1, nuc, weight, p, num);
                }
            }
        }        
#endif
    }
    inFile.close();
    return aa_graphs;
}


template <typename IndexType,
          typename NodeType = pair<IndexType, NumType>,
          typename LatticeType = Lattice<IndexType>,
          typename NucType = IndexType,
          typename NodeNucNodeType = std::tuple<NodeType, NucType, NodeType>>
unordered_map<string, LatticeType> read_wheel_with_weights(const std::string& filename,
        std::unordered_map<std::string, std::unordered_map<NodeType, double, hash_pair>>& nodes_with_best_weight,
        std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, double, std::hash<NodeNucNodeType>>>& edges_with_best_weight,
        const Codon& codon) {
    unordered_map<string, LatticeType> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) 
        throw std::runtime_error("Unable to open coding_wheel file\n");

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;

    for (string line; getline(inFile, line);) {
        stuff = util::split(line, '\t');
        aa = stuff[0];
        LatticeType graph = LatticeType();
        graph.add_node(make_pair(0,0)); // always initialize with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = util::split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_pair(2, i);
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_pair(1, i);
                graph.add_node(n1);
                auto first_num = GET_ACGU_NUC(first);

                double weight = 0.0f;
                if (nodes_with_best_weight[aa].count(make_pair(0, 0))) {
                    weight = edges_with_best_weight[aa][make_tuple(make_pair(0, 0), first_num, n1)] / nodes_with_best_weight[aa][make_pair(0, 0)];
                }

                graph.add_edge(make_pair(0, 0), n1, first_num, weight);
            }
            else {
                n1 = make_pair(1, i-1);
            }
            
            last_first = first;

            auto second_num = GET_ACGU_NUC(second);

            double weight = 0.0f;
            if (nodes_with_best_weight[aa].count(n1)) {
                weight = edges_with_best_weight[aa][make_tuple(n1, second_num, n2)] / nodes_with_best_weight[aa][n1];
            }

            graph.add_edge(n1, n2, second_num, weight);

            for (auto& third : thirds) {

                std::string three_nums = std::string(1, first) + std::string(1, second) + std::string(1, third);

                double weight = 0.0f;
                if (nodes_with_best_weight[aa].count(n2)) {
                    weight = codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2];
                } else {
                    weight = codon.get_weight(aa, three_nums);
                }

                graph.add_edge(n2, make_pair(0,0), GET_ACGU_NUC(third), weight);
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }

    inFile.close();
    return aa_graphs;
}


template <typename IndexType,
          typename NodeType = pair<IndexType, NumType>,
          typename LatticeType = Lattice<IndexType>,
          typename NucType = IndexType,
          typename NodeNucNodeType = std::tuple<NodeType, NucType, NodeType>>
unordered_map<string, LatticeType> read_wheel_with_weights_log(const std::string& filename,
        std::unordered_map<std::string, std::unordered_map<NodeType, double, hash_pair>>& nodes_with_best_weight,
        std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, double, std::hash<NodeNucNodeType>>>& edges_with_best_weight,
        const Codon& codon, double lambda_) {
    unordered_map<string, LatticeType> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) 
        throw std::runtime_error("Unable to open coding_wheel file\n");

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;

    for (string line; getline(inFile, line);) {
        stuff = util::split(line, '\t');
        aa = stuff[0];
        LatticeType graph = LatticeType();
        graph.add_node(make_pair(0,0)); // always initialize with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = util::split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_pair(2, i);
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_pair(1, i);
                graph.add_node(n1);
                auto first_num = GET_ACGU_NUC(first);

                double weight = 1.0f;
                if (nodes_with_best_weight[aa].count(make_pair(0, 0))) {
                    weight = lambda_ * log(edges_with_best_weight[aa][make_tuple(make_pair(0, 0), first_num, n1)] / nodes_with_best_weight[aa][make_pair(0, 0)]);
                }

                graph.add_edge(make_pair(0, 0), n1, first_num, weight);
            }
            else {
                n1 = make_pair(1, i-1);
            }
            
            last_first = first;

            auto second_num = GET_ACGU_NUC(second);

            double weight = 1.0f;
            if (nodes_with_best_weight[aa].count(n1)) {
                weight = lambda_ * log(edges_with_best_weight[aa][make_tuple(n1, second_num, n2)] / nodes_with_best_weight[aa][n1]);
            }

            graph.add_edge(n1, n2, second_num, weight);

            for (auto& third : thirds) {

                std::string three_nums = std::string(1, first) + std::string(1, second) + std::string(1, third);

                double weight = 1.0f;
                if (nodes_with_best_weight[aa].count(n2)) {
                    weight = lambda_ *  log(codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2]);
                } else {
                    weight = lambda_ *  log(codon.get_weight(aa, three_nums));
                }

                graph.add_edge(n2, make_pair(0,0), GET_ACGU_NUC(third), weight);
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }

    inFile.close();
    return aa_graphs;
}

template <typename IndexType,
          typename NucType = IndexType,
          typename NodeType = pair<IndexType, NumType>,
          typename NodeNucNodeType = std::tuple<NodeType, NucType, NodeType>,
          typename WeightType = double,
          typename LatticeType = Lattice<IndexType>>
void prepare_codon_unit_lattice(const std::string& wheel_path, const Codon& codon,
        std::unordered_map<string, LatticeType>& aa_graphs_with_ln_weights_ret,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>&
                best_path_in_one_codon_unit_ret,
        std::unordered_map<std::string, std::string>& aa_best_path_in_a_whole_codon_ret, double lambda_) {

    std::unordered_map<std::string, std::unordered_map<NodeType, WeightType, hash_pair>> nodes_with_best_weight;
    std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, WeightType, std::hash<NodeNucNodeType>>> edges_with_best_weight;

    unordered_map<string, LatticeType> aa_graphs_with_ln_weights;
    unordered_map<string, LatticeType> aa_graphs_with_weights = read_wheel_with_weights<IndexType>(wheel_path, nodes_with_best_weight, edges_with_best_weight, codon);

    for (auto& aa_aa_elem : aa_graphs_with_weights) {
        auto& aa = aa_aa_elem.first;
        auto& aa_elem = aa_aa_elem.second;
        for (auto& node_at_2 : aa_elem.nodes[2]) {
            for (auto& node_at_3_nuc_weight : aa_elem.right_edges[node_at_2]) {
                auto node_at_3 = std::get<0>(node_at_3_nuc_weight);
                auto nuc = std::get<1>(node_at_3_nuc_weight);
                auto weight = std::get<2>(node_at_3_nuc_weight);
                nodes_with_best_weight[aa][node_at_2] = max(nodes_with_best_weight[aa][node_at_2], weight);
                edges_with_best_weight[aa][make_tuple(node_at_2,nuc,node_at_3)] = weight;
            }
        }

        for (auto& node_at_1 : aa_elem.nodes[1]) {
            for (auto& node_at_2_nuc_weight : aa_elem.right_edges[node_at_1]) {
                auto node_at_2 = std::get<0>(node_at_2_nuc_weight);
                auto nuc = std::get<1>(node_at_2_nuc_weight);
                nodes_with_best_weight[aa][node_at_1] = max(nodes_with_best_weight[aa][node_at_1], nodes_with_best_weight[aa][node_at_2]);
                edges_with_best_weight[aa][make_tuple(node_at_1,nuc,node_at_2)] = nodes_with_best_weight[aa][node_at_2];
            }
        }

        for (auto& node_at_0 : aa_elem.nodes[0]) {
            for (auto& node_at_1_nuc_weight : aa_elem.right_edges[node_at_0]) {
                auto node_at_1 = std::get<0>(node_at_1_nuc_weight);
                auto nuc = std::get<1>(node_at_1_nuc_weight);
                nodes_with_best_weight[aa][node_at_0] = max(nodes_with_best_weight[aa][node_at_0], nodes_with_best_weight[aa][node_at_1]);
                edges_with_best_weight[aa][make_tuple(node_at_0,nuc,node_at_1)] = nodes_with_best_weight[aa][node_at_1];
            }
        }
    }

    aa_graphs_with_ln_weights = read_wheel_with_weights_log<IndexType>(wheel_path,  nodes_with_best_weight, edges_with_best_weight, codon, lambda_);

    std::unordered_map<std::string, 
                       std::unordered_map<std::tuple<NodeType, NodeType>, 
                       std::tuple<double, NucType, NucType>,
                       std::hash<std::tuple<NodeType, NodeType>>>>
                       best_path_in_one_codon_unit;


    for (auto& aa_graph : aa_graphs_with_ln_weights) {
        auto& aa = aa_graph.first;
        auto& graph = aa_graph.second;
        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_nuc_log_w : graph.right_edges[node_0]) {
                auto node_1 = std::get<0>(node_1_nuc_log_w);
                auto nuc = std::get<1>(node_1_nuc_log_w);
                auto log_weight = std::get<2>(node_1_nuc_log_w);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_1)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)] = make_tuple(log_weight,nuc,k_void_nuc);
                }
            }
        }

        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_nuc_log_w : graph.right_edges[node_1]) {
                auto node_2 = std::get<0>(node_2_nuc_log_w);
                auto nuc = std::get<1>(node_2_nuc_log_w);
                auto log_weight = std::get<2>(node_2_nuc_log_w);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_2)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)] = make_tuple(log_weight,nuc,k_void_nuc);
                }

                auto temp = best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)];
            }
        }

        for (auto& node_2 : graph.nodes[2]) {
            for (auto& node_3_nuc_log_w : graph.right_edges[node_2]) {
                auto node_3 = std::get<0>(node_3_nuc_log_w);
                auto nuc = std::get<1>(node_3_nuc_log_w);
                auto log_weight = std::get<2>(node_3_nuc_log_w);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_2,node_3)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)] = make_tuple(log_weight,nuc,k_void_nuc);
                }
            }
        }

        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_nuc_0_log_weight_0 : graph.right_edges[node_0]) {
                auto& node_1 = std::get<0>(node_1_nuc_0_log_weight_0);
                auto& nuc_0 = std::get<1>(node_1_nuc_0_log_weight_0);
                auto log_weight_0 = std::get<2>(node_1_nuc_0_log_weight_0);
                for (auto& node_2_nuc_1_log_weight_1 : graph.right_edges[node_1]) {
                    auto& node_2 = std::get<0>(node_2_nuc_1_log_weight_1);
                    auto& nuc_1 = std::get<1>(node_2_nuc_1_log_weight_1);
                    auto log_weight_1 = std::get<2>(node_2_nuc_1_log_weight_1);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_2)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)]) < log_weight_0 + log_weight_1)
                        best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)] = make_tuple(log_weight_0 + log_weight_1, nuc_0, nuc_1);
                }
            }
        }

        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_nuc_1_log_weight_1 : graph.right_edges[node_1]) {
                auto& node_2 = std::get<0>(node_2_nuc_1_log_weight_1);
                auto& nuc_1 = std::get<1>(node_2_nuc_1_log_weight_1);
                auto log_weight_1 = std::get<2>(node_2_nuc_1_log_weight_1);
                for (auto& node_3_nuc_2_log_weight_2 : graph.right_edges[node_2]) {
                    auto& node_3 = std::get<0>(node_3_nuc_2_log_weight_2);
                    auto& nuc_2 = std::get<1>(node_3_nuc_2_log_weight_2);
                    auto log_weight_2 = std::get<2>(node_3_nuc_2_log_weight_2);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_3)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)]) < log_weight_1 + log_weight_2)
                        best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)] = make_tuple(log_weight_1 + log_weight_2, nuc_1, nuc_2);
                }
            }
        }
    }

    std::unordered_map<std::string, double> max_path;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;

    for (auto& aa_path_weight : codon.aa_table_) {
        auto& aa = aa_path_weight.first; // char
        for (auto& path_weight : aa_path_weight.second) {
            if (max_path[aa] < path_weight.second) {
                max_path[aa] = path_weight.second;
                aa_best_path_in_a_whole_codon[aa] = path_weight.first;
            }
        }
    }

    aa_graphs_with_ln_weights_ret = aa_graphs_with_ln_weights;
    best_path_in_one_codon_unit_ret = best_path_in_one_codon_unit;
    aa_best_path_in_a_whole_codon_ret = aa_best_path_in_a_whole_codon;
}



template <typename IndexType,
          typename NodeType = pair<IndexType, NumType>,
          typename LatticeType = Lattice<IndexType>,
          typename DFAType = DFA<IndexType>>
DFAType get_dfa(unordered_map<string, LatticeType> aa_graphs, vector<string> aa_seq) {
    DFAType dfa = DFAType();
    NodeType newnode = make_pair(3 * static_cast<IndexType>(aa_seq.size()), 0);
    dfa.add_node(newnode);
    IndexType i = 0;
    IndexType i3;
    string aa;
    LatticeType graph;
    for(auto& item : aa_seq) {
        i3 = i * 3;
        aa = aa_seq[i];
        graph = aa_graphs[aa];
        for (IndexType pos = 0; pos <= 2; pos++) {
            for(auto& node : graph.nodes[pos]) {
                IndexType num = get<1>(node);
                newnode = make_pair(i3 + pos, num);
                dfa.add_node(newnode);
                for (auto& edge : graph.right_edges[node]) {
                    NodeType n2 = get<0>(edge);
                    IndexType nuc = get<1>(edge);
                    num = get<1>(n2);
                    NodeType newn2 = make_pair(i3 + pos + 1, num);
                    dfa.add_edge(newnode, newn2, nuc, get<2>(edge));
                }
            }
        }
        i++;
    }
#ifdef is_verbose
    printf("-----------------DFA------------------------\n");
    for(IndexType pos = 0; pos < 3 * static_cast<IndexType>(aa_seq.size()) + 1; pos++){
        for(auto& node : dfa.nodes[pos]) {
            IndexType p = get<0>(node);
            IndexType num = get<1>(node);
            printf("node, (%d, %d)\n", p, num);
            for(auto &n2 : dfa.auxiliary_right_edges[node]){
                IndexType p2 = get<0>(n2.first);
                IndexType num2 = get<1>(n2.first);
                for(auto nuc : n2.second){
                    printf("              (%d, %d) -(%d,%lf)-> (%d, %d)\n", p, num, get<0>(nuc),get<1>(nuc), p2, num2);
                }
            }
            for(auto &n1 : dfa.auxiliary_left_edges[node]){
                IndexType p1 = get<0>(n1.first); IndexType num1 = get<1>(n1.first);
                for(auto nuc : n1.second){
                    printf("  (%d, %d) <-(%d,%lf)- (%d, %d)\n", p1, num1, get<0>(nuc),get<1>(nuc), p, num);
                }
            }
        }
    }
#endif
    return dfa;
}

}

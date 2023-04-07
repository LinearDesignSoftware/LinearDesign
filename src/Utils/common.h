#ifndef common_h
#define common_h

#include <utility>
#include <functional>
#include <array>
#include <string>
#include <set>
#include <map>
#include <exception>
#include <list>
#include "base.h"

namespace LinearDesign {


using SizeType               = size_t;
using ScoreType              = int32_t;
using IndexType              = int32_t; //if less than 10000, only int16_t is needed here
using NucType                = int8_t;
using NumType                = int32_t;
using NucPairType            = int8_t;
using PairType               = int8_t;
using FinalScoreType         = double;
using NodeType               = std::pair<IndexType, NumType>;
using NodeNucType            = std::pair<NodeType, NucType>;
using NodeNucWType           = std::tuple<NodeType, NucType, double>;
using PairType               = int8_t;




enum class Manner : std::uint8_t {
    NONE = 0,              // 0: empty
    H,                     // 1: hairpin candidate
    HAIRPIN,               // 2: hairpin
    SINGLE,                // 3: single
    HELIX,                 // 4: helix
    MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
    MULTI_eq_MULTI_plus_U, // 6: multi = multi + U
    P_eq_MULTI,            // 7: P = (multi)
    M2_eq_M_plus_P,        // 8: M2 = M + P
    M_eq_M2,               // 9: M = M2
    M_eq_M_plus_U,         // 10: M = M + U
    M_eq_P,                // 11: M = P
    C_eq_C_plus_U,         // 12: C = C + U
    C_eq_C_plus_P,         // 13: C = C + P
};

enum class Beam_type : std::uint8_t {
    BEAM_C = 0,
    BEAM_P,
    BEAM_MULTI,
    BEAM_M2,
    BEAM_M1

};


template <typename ScoreType>
struct State {
    ScoreType score = util::value_min<ScoreType>();
    double cai_score = util::value_min<double>();
    NodeType pre_node;
    double pre_left_cai;
};

struct BacktraceResult {
    std::string seq;
    std::string structure;
};

template <typename ScoreType,
          typename IndexType,
          typename NodeType = std::pair<IndexType, IndexType>>
struct DecoderResult {
    std::string sequence;
    std::string structure;
    ScoreType score;
    ScoreType cai;
    ScoreType old_cai;
    IndexType num_states;
};

template <typename ScoreType, 
          typename IndexType, 
          typename NodeType = std::pair<IndexType, IndexType>>
struct ScoreInnerDate {
    ScoreType newscore;
    NodeType j_node;
    NodeType i_node;
    int nuc_pair;
};


struct NodeNucpair {
    IndexType node_first;
    NumType node_second;
    NucPairType nucpair;
};


}

#endif /* common_h */

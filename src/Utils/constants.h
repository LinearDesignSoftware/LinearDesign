#ifndef constants_h
#define constants_h

#include <map>

namespace LinearDesign {

constexpr uint8_t k_void_nuc = 127;

// static std::map<char, std::string> k_map_1_3 = {
//     {'F',"Phe"},
//     {'L',"Leu"},
//     {'S',"Ser"},
//     {'Y',"Tyr"},
//     {'*',"STOP"},
//     {'C',"Cys"},
//     {'W',"Trp"},
//     {'P',"Pro"},
//     {'H',"His"},
//     {'Q',"Gln"},
//     {'R',"Arg"},
//     {'I',"Ile"},
//     {'M',"Met"},
//     {'T',"Thr"},
//     {'N',"Asn"},
//     {'K',"Lys"},
//     {'V',"Val"},
//     {'D',"Asp"},
//     {'E',"Glu"},
//     {'G',"Gly"},
//     {'A',"Ala"}
// };

static std::map<std::string, char> k_map_3_1 = {
    {"Phe", 'F'},
    {"Leu", 'L'},
    {"Ser", 'S'},
    {"Tyr", 'Y'},
    {"STOP", '*'},
    {"Cys", 'C'},
    {"Trp", 'W'},
    {"Pro", 'P'},
    {"His", 'H'},
    {"Gln", 'Q'},
    {"Arg", 'R'},
    {"Ile", 'I'},
    {"Met", 'M'},
    {"Thr", 'T'},
    {"Asn", 'N'},
    {"Lys", 'K'},
    {"Val", 'V'},
    {"Asp", 'D'},
    {"Glu", 'E'},
    {"Gly", 'G'},
    {"Ala", 'A'}
};

}

#endif /* constants_h */

#ifndef codon_h
#define codon_h

#include <exception>
#include <fstream>
#include <string>
#include <utility>
#include <map>
#include <vector>
#include <regex>
#include <cmath>

#include <typeinfo>

#include "base.h"
#include "constants.h"

namespace LinearDesign {

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}


class Codon {
public:
	Codon(const std::string& path) : codon_table_(), aa_table_() {
        std::ifstream codon_file;
        codon_file.open(path);
        if (codon_file.is_open()) {
            int index = 0;
            for (std::string line; getline(codon_file, line);){

                rtrim(line);

                if(line.size() == 0 or line.empty())
                    continue;

                if (index++ == 0)
                    continue;
                
                const auto line_split = util::split(line, ',');
                if(line_split.size() != 3){
                    std::cerr << "Wrong format of codon frequency file!" << std::endl;
                    exit(1);
                }
                const std::string codon = line_split[0];
                const std::string aa = line_split[1];
                const float fraction = std::stof(line_split[2]);

                codon_table_[codon] = make_pair(aa, fraction);
                aa_table_[aa].push_back(make_pair(codon, fraction));
                if (!max_aa_table_.count(aa))
                    max_aa_table_[aa] = fraction;
                else
                    max_aa_table_[aa] = fmax(max_aa_table_[aa], fraction);
            }
            codon_file.close();
            if (codon_table_.size() != 64){
                std::cerr << "Codon frequency file needs to contain 64 codons!" << std::endl;
                exit(1);
            }

        } else {
            std::cerr << "The input codon frequency file does not exist!" << std::endl;
            exit(1);
        }
    }

    float calc_cai(const std::string& rna_seq) const {
        if (rna_seq.length() % 3)
            throw std::runtime_error("invalid rna seq");
        
        const int protein_length = static_cast<int>(rna_seq.length() / 3);
        float cai = 0.0f;
        
        for (int index = 3; index < rna_seq.length() + 1; index += 3) {
            const std::string tri_letter = rna_seq.substr(index - 3, 3);
            const auto f_ci_aa = codon_table_.at(tri_letter);
            const auto f_c_max = max_aa_table_.at(f_ci_aa.first);
            
            float w_i = f_ci_aa.second / f_c_max;
            cai += log2f(w_i);
        }
        
        return exp2f(cai / protein_length);
    }

    std::string find_max_codon(const char aa, 
    						   const std::string& match) const {
    	auto candidate_condons = aa_table_.at(std::string(1, aa));

    	float max_score = 0;
    	std::string max_condon;
    	for (auto& candidate : candidate_condons) {
    		if (std::regex_match(candidate.first, std::regex(match)) && 
    				candidate.second > max_score) {
    			max_condon = candidate.first;
    			max_score = candidate.second;
    		}
    	}

    	if (max_condon.empty())
    		throw std::runtime_error("invald search");

        // assert(codon_table_.at(max_condon).first == std::string(1, aa));
    	return max_condon;
    }

    std::string cvt_rna_seq_to_aa_seq(const std::string& rna_seq) const {
        if (rna_seq.length() % 3)
            throw std::runtime_error("invalid rna seq");

        std::string aa_seq;
        aa_seq.reserve(rna_seq.length());
        for (int index = 3; index < rna_seq.length() + 1; index += 3) {
            const std::string tri_letter = rna_seq.substr(index - 3, 3);
            auto aa = codon_table_.at(tri_letter).first;
            if (aa == "STOP") {
                aa_seq.append("*");
                return aa_seq;
            }
            aa_seq.append(codon_table_.at(tri_letter).first);
        }
        return aa_seq;
    }

    float get_weight(const std::string& aa_tri, const std::string& codon) const {

        if (k_map_3_1.count(aa_tri)) {
            auto codons = aa_table_.at(std::string(1, k_map_3_1[aa_tri]));
            auto it = std::find_if(codons.begin(), codons.end(), [codon](const std::pair<std::string, float>& e){
                // std::cout << typeid(e).name() << '\n';
                return e.first == codon;
            });
            if (it == codons.end()) {
                throw std::runtime_error("invalid codon");
            }
            return it->second;
        } else if (three_prime_aa_table_.count(aa_tri)) {
            return three_prime_aa_table_.at(aa_tri).second;
        }

        return 0.0f;
    }



// private:
    std::vector<std::string> aux_aa_;
    std::map<std::string, std::pair<std::string, float>> three_prime_codon_table_;
    std::map<std::string, std::pair<std::string, float>> three_prime_aa_table_;

    
    std::map<std::string, float> max_aa_table_;
	std::map<std::string, std::pair<std::string, float>> codon_table_;
	std::map<std::string, std::vector<std::pair<std::string, float>>> aa_table_;
};

}

#endif 

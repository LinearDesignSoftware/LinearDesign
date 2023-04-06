#ifndef fasta_h
#define fasta_h

#include <exception>
#include <map>
#include <fstream>
#include <string>

#include "base.h"

namespace LinearDesign {

struct Reader {
	static bool cvt_to_seq(const string& from, string& to) {
		return false;
	}
};

struct Fasta : public Reader {
	static map<char, string> map_fasta;

	static bool cvt_to_seq(const string& fasta, string& nucs) {
	    nucs.reserve(4 * fasta.length());
	    for(auto aa : fasta) {
	        if (map_fasta.count(aa)) {
	            nucs.append(map_fasta[aa] + " ");
	        } else {
	            cerr << "invalid protein sequence!\n" << endl;
	            return false;
	        }
	    }
	    nucs.pop_back();
	    return true;
	}
};

map<char, string> Fasta::map_fasta = {
	{'F',"Phe"}, 
	{'L',"Leu"}, 
	{'S',"Ser"}, 
	{'Y',"Tyr"}, 
	{'*',"STOP"}, 
	{'C',"Cys"}, 
	{'W',"Trp"}, 
	{'P',"Pro"}, 
	{'H',"His"}, 
	{'Q',"Gln"}, 
	{'R',"Arg"}, 
	{'I',"Ile"}, 
	{'M',"Met"}, 
	{'T',"Thr"}, 
	{'N',"Asn"}, 
	{'K',"Lys"}, 
	{'V',"Val"}, 
	{'D',"Asp"}, 
	{'E',"Glu"}, 
	{'G',"Gly"}, 
	{'A',"Ala"}
};

template <class T>
struct ReaderTraits {
	static bool cvt_to_seq(const string& from, string& to) {
		return T::cvt_to_seq(from, to);
	}
};

}

#endif /* fasta_h */
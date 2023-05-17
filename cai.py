import pandas as pd
from math import log2

codon_tab = pd.read_csv("./CAI_table.csv", index_col=0)


def calc_cai(transcript):
    try:
        cai = 0.0
        transcript = transcript.upper()
        for i in range(0, len(transcript), 3):
            codon = transcript[i : i + 3]  # POSIBLE ERROR: case sensitive
            w_i = codon_tab.loc[codon, "X"] / codon_tab.loc[codon, "c_max"]
            cai += log2(w_i)

        answer = 2 ** (cai / (len(transcript) / 3 - 1))
    except KeyError:
        print("ERROR: Invalid codon in transcript")
        raise
    return answer

"""
Offline substitution scores.

BLOSUM62 — NCBI BLAST matrix (ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62).
Grantham  — chemical distance, Grantham 1974 (AAindex GRAR740104).

References
----------
Grantham, R. (1974). Amino acid difference formula to help explain protein
    evolution. Science, 185(4154), 862–864.
    https://doi.org/10.1126/science.185.4154.862

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from
    protein blocks. Proceedings of the National Academy of Sciences, 89(22),
    10915–10919. https://doi.org/10.1073/pnas.89.22.10915

Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., &
    Kanehisa, M. (2008). AAindex: Amino acid index database, progress report
    2008. Nucleic Acids Research, 36(Database issue), D202–D205.
    https://doi.org/10.1093/nar/gkm998
"""
from __future__ import annotations

AA = "ARNDCQEGHILKMFPSTWYV"

_BLOSUM62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

# Lower triangle, rows/cols in AA order.
_GRANTHAM_TRI = """
0
112 0
111 86 0
126 96 23 0
195 180 139 154 0
91 43 46 61 154 0
107 54 42 45 170 29 0
60 125 80 94 159 87 98 0
86 29 68 81 174 24 40 98 0
94 97 149 168 198 109 134 135 94 0
96 102 153 172 198 113 138 138 99 5 0
106 26 94 101 202 53 56 127 32 102 107 0
84 91 142 160 196 101 126 127 87 10 15 95 0
113 97 158 177 205 116 140 153 100 21 22 102 28 0
27 103 91 108 169 76 93 42 77 95 98 103 87 114 0
99 110 46 65 112 68 80 56 89 142 145 121 135 155 74 0
58 71 65 85 149 42 65 59 47 89 92 78 81 103 38 58 0
148 101 174 181 215 130 152 184 115 61 61 110 67 40 147 177 128 0
112 77 143 160 194 99 122 147 83 33 36 85 36 22 110 144 92 37 0
64 96 133 152 192 96 121 109 84 29 32 97 21 50 68 124 69 88 55 0
"""

BLOSUM62: dict[tuple[str, str], int] = {
    (a, b): int(v)
    for a, row in zip(AA, _BLOSUM62_ROWS.split("\n")[1:])
    for b, v in zip(AA, row.split())
}

GRANTHAM: dict[tuple[str, str], int] = {}
for _i, _row in enumerate(_GRANTHAM_TRI.strip().split("\n")):
    for _j, _v in enumerate(_row.split()):
        GRANTHAM[(AA[_i], AA[_j])] = GRANTHAM[(AA[_j], AA[_i])] = int(_v)


def substitution_scores(ref: str | None, alt: str | None) -> dict:
    """BLOSUM62 and Grantham for one residue change; None where undefined."""
    key = ((ref or "").upper(), (alt or "").upper())
    return {"blosum62": BLOSUM62.get(key), "grantham": GRANTHAM.get(key)}

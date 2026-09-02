#!/usr/bin/env python3
"""
HyPR01 chloroplast genome: rotate the validated 139,646 bp genome to the
trnH-GUG / LSC boundary. Supersedes build_final_chloroplast_genome.py.

Provenance (see HyPR01_Chloroplast_Assembly_Methods_and_Audit.md, section 10):

HyPR01_Chloroplast_Genome.fasta (139,646 bp, Bandage-curated) is validated as
correct: TIPPo round 1's independent, fully-automatic Flye assembly (no
manual editing) produced a circular 139,646 bp sequence that is byte-for-byte
identical to it across the full length. An earlier attempt (see the
deprecated build_final_chloroplast_genome.py) incorrectly removed 37 bp of
real sequence based on a flawed comparison against round 2's incomplete,
2-fragment assembly. That "fix" has been retracted.

This script only rotates -- it does not remove anything. In the original
curated genome, psbA (110,377-111,438) occurs *before* trnH-GUG
(111,989-112,063), the reverse of the usual deposited-plastome convention
(trnH-GUG as gene 1, psbA immediately downstream). Rotation alone can't fix
that (it preserves relative gene order), so this reverse-complements the
whole circular molecule, then rotates so trnH-GUG starts at position 1.

Input:  HyPR01_Chloroplast_Genome.fasta       (139,646 bp)
Output: HyPR01_Chloroplast_Genome_rotated.fasta (139,646 bp, trnH-GUG at 1-75)
"""

BASE = "/path/to/your/directory"


def revcomp(s):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def read_fasta(path):
    with open(path) as f:
        lines = f.readlines()
    return "".join(l.strip() for l in lines[1:])


def main():
    curated_path = f"{BASE}/HyPR01_Chloroplast_Genome.fasta"
    out_path = f"{BASE}/HyPR01_Chloroplast_Genome_rotated.fasta"

    curated = read_fasta(curated_path)
    N = len(curated)
    assert N == 139646, f"unexpected input length {N}, refusing to proceed blind"

    # Consensus coordinates from GeSeq (ARAGORN + tRNAscan-SE agree on trnH-GUG;
    # blatX + HMMER agree on psbA), both '+' strand, 1-based inclusive.
    trnH_seq = curated[111988:112063]   # 111,989-112,063 (1-based)
    psbA_seq = curated[110376:111438]   # 110,377-111,438 (1-based)
    assert len(trnH_seq) == 75
    assert len(psbA_seq) == 1062

    idx_trnH = curated.find(trnH_seq)
    idx_psbA = curated.find(psbA_seq)
    assert idx_trnH != -1 and idx_psbA != -1

    rc = revcomp(curated)
    trnH_rc_start = N - (idx_trnH + len(trnH_seq))
    psbA_rc_start = N - (idx_psbA + len(psbA_seq))
    # verify via direct search, not just arithmetic
    assert rc.find(revcomp(trnH_seq)) == trnH_rc_start
    assert rc.find(revcomp(psbA_seq)) == psbA_rc_start
    print(f"[1] After reverse-complementing: trnH-GUG at {trnH_rc_start}, psbA at {psbA_rc_start} "
          f"({psbA_rc_start - (trnH_rc_start + len(trnH_seq))} bp downstream of trnH-GUG)")

    rotated = rc[trnH_rc_start:] + rc[:trnH_rc_start]
    assert len(rotated) == N
    assert rotated.startswith(revcomp(trnH_seq))
    new_psbA_pos = rotated.find(revcomp(psbA_seq))
    print(f"[2] Rotated: trnH-GUG at position 1 (1-based); "
          f"psbA at position {new_psbA_pos + 1} (1-based)")

    with open(out_path, "w") as f:
        f.write(f">HyPR01_Chloroplast_Genome_rotated circular length={len(rotated)}\n")
        for i in range(0, len(rotated), 70):
            f.write(rotated[i:i + 70] + "\n")
    print(f"\nWrote {out_path} ({len(rotated)} bp)")


if __name__ == "__main__":
    main()

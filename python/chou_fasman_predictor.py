"""
Chou-Fasman Secondary Structure Prediction (Helix H / Beta Strand S)
--------------------------------------------------------------------
Method used: 
1) Helix (H):
   - Nucleation: scan 6-res windows; a window is a seed if >=4 residues have Pa > 1.0
   - Extension: extend left/right by checking 4-res windows (new residue + nearest 3 in the segment);
                continue extension if sum(Pa) for those 4 residues >= 4.0

2) Strand (S):
   - Nucleation: scan 5-res windows; a window is a seed if >=3 residues have Pb > 1.0
   - Extension: extend left/right by checking 4-res windows (as above but with Pb);
                continue if sum(Pb) for those 4 residues > 4.0

3) Conflict Resolution:
   - For overlapping H and S assignments, compute average Pa and average Pb over the
     exact overlapping region and assign the higher (H or S) to those residues.

Outputs:
- Helical regions (start-end, sequence)
- Strand regions (start-end, sequence)
- Conflicts and how each was resolved
- Final assignment string (H/S/'-' where unassigned)

Indexing in outputs is 1-based
"""

from typing import List, Tuple

# Input protein sequence SEQ
SEQ = """MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYLEKPMEIARIVARCLWEESRLLQTAATAAQQGGQANHPTAAVVTEKQQMLEQHLQDVRKRVQDLEQKMKVVENLQDDFDFNYKTLKSQGDMQDLNGNNQSVTRQKMQQLEQMLTALDQMRRSIVSELAGLLSAMEYVQKTLTDEELADWKRRQQIACIGGPPNICLDRLENWITSLAESQLQTRQQIKKLEELQQKVSYKGDPIVQHRPMLEERIVELFRNLMKSAFVVERQPCMPMHPDRPLVIKTGVQFTTKVRLLVKFPELNYQLKIKVCIDKDSGDVAALRGSRKFNILGTNTKVMNMEESNNGSLSAEFKHLTLREQRCGNGGRANCDASLIVTEELHLITFETEVYHQGLKIDLETHSLPVVVISNICQMPNAWASILWYNMLTNNPKNVNFFTKPPIGTWDQVAEVLSWQFSSTTKRGLSIEQLTTLAEKLLGPGVNYSGCQITWAKFCKENMAGKGFSFWVWLDNIIDLVKKYILALWNEGYIMGFISKERERAILSTKPPGTFLLRFSESSKEGGVTFTWVEKDISGKTQIQSVEPYTKQQLNNMSFAEIIMGYKIMDATNILVSPLVYLYPDIPKEEAFGKYCRPESQEHPEADPGSAAPYLKTKFICVTPTTCSNTIDLPMSPRTLDSLMQFGNNGEGAEPSAGGQFESLTFDMELTSECATSPM
"""
SEQ = "".join([c for c in SEQ.upper() if c.isalpha()]) 
N = len(SEQ)

# Chou-Fasman parameters
P_ALPHA = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17,
    'W': 1.14, 'V': 1.14, 'F': 1.12, 'K': 1.07, 'I': 1.00, 'D': 0.98,
    'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77, 'N': 0.73, 'Y': 0.61,
    'P': 0.59, 'G': 0.53
}

P_BETA = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28,
    'Q': 1.23, 'L': 1.22, 'T': 1.20, 'W': 1.19, 'A': 0.97, 'R': 0.90,
    'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65,
    'P': 0.62, 'E': 0.26
}

def pa(res: str) -> float:
    #Return helix propensity Pa for residue; default 0.0 if unknown
    return P_ALPHA.get(res, 0.0)

def pb(res: str) -> float:
    #Return strand propensity Pb for residue; default 0.0 if unknown
    return P_BETA.get(res, 0.0)

# Helpers
def window_has_min_over_1(vals: List[float], min_count: int) -> bool:
    #True if at least min_count entries exceed 1.0 in the given list.
    return sum(v > 1.0 for v in vals) >= min_count

def extend_segment(start: int, end: int, prop_list: List[float], threshold_sum: float, seq_len: int) -> Tuple[int, int]:
    """
    Extend a segment [start, end] (inclusive, 0-based) to left/right
    At each step:
      - Try left: 4-res window = [new_left] + next three residues already in segment
      - Try right: 4-res window = last three residues in segment + [new_right]
    Continue extending while the sum over the 4-res window >= threshold_sum
    """
    changed = True
    while changed:
        changed = False
        # Left extension
        if start > 0 and end - start + 1 >= 3:
            w_left = [prop_list[start - 1]] + prop_list[start:start + 3]
            if sum(w_left) >= threshold_sum:
                start -= 1
                changed = True
        # Right extension
        if end < seq_len - 1 and end - start + 1 >= 3:
            w_right = prop_list[end - 2:end + 1] + [prop_list[end + 1]]
            if sum(w_right) >= threshold_sum:
                end += 1
                changed = True
    return start, end

def merge_segments(segments: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    #Merge overlapping/adjacent inclusive segments
    if not segments:
        return []
    segments = sorted(segments)
    merged = [segments[0]]
    for s, e in segments[1:]:
        ms, me = merged[-1]
        if s <= me + 1:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged

def segments_from_labels(labels: List[str], target: str) -> List[Tuple[int, int]]:
    #Get contiguous [start,end] (0-based, inclusive) where labels==target
    out = []
    i = 0
    n = len(labels)
    while i < n:
        if labels[i] == target:
            j = i
            while j + 1 < n and labels[j + 1] == target:
                j += 1
            out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out

def format_segments(seq: str, segments: List[Tuple[int, int]]) -> List[str]:
    #Formatted lines with 1-based indexing and the subsequence
    lines = []
    for s, e in segments:
        lines.append(f"{s+1}-{e+1}  len={e-s+1}  {seq[s:e+1]}")
    return lines

# Core prediction
Palpha_list = [pa(r) for r in SEQ]
Pbeta_list  = [pb(r) for r in SEQ]

# Helix nucleation (6-mer, >=4 residues with Pa>1)
helix_segments = []
wH = 6
for i in range(0, N - wH + 1):
    vals = Palpha_list[i:i + wH]
    if window_has_min_over_1(vals, 4):
        s, e = i, i + wH - 1
        s, e = extend_segment(s, e, Palpha_list, threshold_sum=4.0, seq_len=N)
        helix_segments.append((s, e))
helix_segments = merge_segments(helix_segments)

# Strand nucleation (5-mer, >=3 residues with Pb>1)
strand_segments = []
wS = 5
for i in range(0, N - wS + 1):
    vals = Pbeta_list[i:i + wS]
    if window_has_min_over_1(vals, 3):
        s, e = i, i + wS - 1
        # For b, strict '>' 4.0
        s, e = extend_segment(s, e, Pbeta_list, threshold_sum=4.0000001, seq_len=N)
        strand_segments.append((s, e))
strand_segments = merge_segments(strand_segments)

# Initial labels
labels_H = ['-'] * N
for s, e in helix_segments:
    for k in range(s, e + 1):
        labels_H[k] = 'H'

labels_S = ['-'] * N
for s, e in strand_segments:
    for k in range(s, e + 1):
        labels_S[k] = 'S'

# Find conflicts
conflict_mask = [(labels_H[i] == 'H' and labels_S[i] == 'S') for i in range(N)]
conflict_regions = []
i = 0
while i < N:
    if conflict_mask[i]:
        j = i
        while j + 1 < N and conflict_mask[j + 1]:
            j += 1
        conflict_regions.append((i, j))
        i = j + 1
    else:
        i += 1

# Resolve conflicts 
final_labels = ['-'] * N
for i in range(N):
    if labels_H[i] == 'H':
        final_labels[i] = 'H'
    if labels_S[i] == 'S':
        final_labels[i] = 'S' if final_labels[i] == '-' else 'X'  # mark overlap

resolved_choices = []  # (start, end, chosen, avg_Pa, avg_Pb)
for s, e in conflict_regions:
    pa_sum = sum(Palpha_list[s:e+1])
    pb_sum = sum(Pbeta_list[s:e+1])
    length = e - s + 1
    avg_pa = pa_sum / length
    avg_pb = pb_sum / length
    chosen = 'S' if avg_pb > avg_pa else 'H'
    for k in range(s, e + 1):
        final_labels[k] = chosen
    resolved_choices.append((s, e, chosen, avg_pa, avg_pb))

# Clean any accidental 'X' (should not remain after the loop above)
for i in range(N):
    if final_labels[i] == 'X':
        final_labels[i] = 'H'  # Give preference to H if any residual tie marker

# Report
helix_out  = format_segments(SEQ, helix_segments)
strand_out = format_segments(SEQ, strand_segments)

print(f"Sequence length: {N}\n")

print("HELICAL REGIONS (H):")
if helix_out:
    for line in helix_out:
        print("  ", line)
else:
    print("  None")

print("\nBETA-STRAND REGIONS (S):")
if strand_out:
    for line in strand_out:
        print("  ", line)
else:
    print("  None")

print("\nCONFLICTING REGIONS and RESOLUTION:")
if conflict_regions:
    for s, e, chosen, avg_pa, avg_pb in [
        (c[0], c[1], next(r[2] for r in resolved_choices if r[0]==c[0] and r[1]==c[1]),
         next(r[3] for r in resolved_choices if r[0]==c[0] and r[1]==c[1]),
         next(r[4] for r in resolved_choices if r[0]==c[0] and r[1]==c[1]))
        for c in conflict_regions
    ]:
        print(f"  {s+1}-{e+1}  len={e-s+1}  avg Pa={avg_pa:.3f}  avg Pb={avg_pb:.3f}  -> assign {chosen}   seg={SEQ[s:e+1]}")
else:
    print("  None")

final_str = "".join(final_labels)
print("\nFINAL ASSIGNMENT:")
print(final_str)

print("\nSequence / Assignment (blocks of 60):")
for i in range(0, N, 60):
    block_seq = SEQ[i:i+60]
    block_ann = final_str[i:i+60]
    print(f"{i+1:4d}-{i+len(block_seq):4d}  {block_seq}")
    print(" " * 12 + block_ann)


# Chou-Fasman Secondary Structure Prediction

This project implements the **Chou–Fasman algorithm** to predict secondary structure elements (α-helices and β-strands) from an input amino acid sequence. The algorithm performs nucleation, directional extension, and conflict resolution to generate a final H/S/coil assignment for every residue.

---

## Features

- Predicts:
  - **Helical regions (H)**
  - **Beta-strand regions (S)**
  - **Coil / Non-structured regions**
- Implements:
  - Helix & strand nucleation windows
  - Sliding-window extension rules
  - Overlap conflict resolution (avg Pa vs avg Pb)
- Outputs:
  - Start/end indices of each predicted region
  - Region lengths and sequences
  - Final secondary structure line annotation

---

## Algorithm Summary

| Step | Description |
|------|-------------|
| Helix Nucleation | Scan 6-residue windows; ≥4 residues must have **Pa > 1.0** |
| Helix Extension  | Extend left/right while sliding 4-residue windows have **ΣPa >= 4.0** |
| Strand Nucleation | Scan 5-residue windows; ≥3 residues must have **Pb > 1.0** |
| Strand Extension | Extend while sliding 4-residue windows have **ΣPb > 4.0** |
| Conflict Resolution | For overlapping H/S assignments, compare **avg(Pa)** vs **avg(Pb)** |

---

## ▶Running the Program

Run directly with Python:

```bash
python3 chou_fasman_predictor.py

````

Save output to a file:

```bash
python3 chou_fasman_predictor.py > output.txt
```

---

## Example Output (Excerpt)

```
HELICAL REGIONS (H):
  1-23   MAQWNQLQQLDTRYLEQLHQLYS
  25-67  SFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQ
  ...

STRAND REGIONS (S):
  52-95  HATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYL
  ...
```

---

## License

Released under the **MIT License**
Free for academic, personal, and commercial use.

---


**Evani Menon**
B.Tech — IIIT Delhi

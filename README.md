# Chou–Fasman Secondary Structure Prediction

This project implements the **Chou–Fasman algorithm** to predict secondary structure elements —
**alpha-helices (H)** and **beta-strands (S)** — directly from an input protein sequence.

The algorithm performs:

* **Nucleation** (seed detection)
* **Directional extension**
* **Conflict resolution** between helix/strand predictions

It outputs both **region-level predictions** and a **final residue-level annotation**.

---

## Features

* Predicts secondary structure regions:

  * **Helices (H)**
  * **Beta-strands (S)**
  * **Coil / unassigned ('-')**
* Implements full Chou–Fasman rules:

  * Helix & strand nucleation windows
  * Sliding-window extension in both directions
  * Final conflict resolution (avg **Pa** vs **Pb**)
* Outputs:

  * Start–end positions of each predicted H/S region
  * Region lengths and residue sequences
  * Full final structure line
    (e.g., `HHHSSSHH---HH...`)

---

## Algorithm Summary

| Step                    | Logic                                                              |
| ----------------------- | ------------------------------------------------------------------ |
| **Helix Nucleation**    | Scan **6-residue** windows. Valid if **>=4 residues have Pa > 1.0** |
| **Helix Extension**     | Extend left/right if **ΣPa (4-res window) >= 4.0**                  |
| **Strand Nucleation**   | Scan **5-residue** windows. Valid if **>=3 residues have Pb > 1.0** |
| **Strand Extension**    | Extend if **ΣPb (4-res window) > 4.0**                             |
| **Conflict Resolution** | Compare **avg(Pa)** vs **avg(Pb)** over overlap → choose H or S    |

---

## Running the Python Program

Run normally:
```bash
python3 chou_fasman_predictor.py
```

Save results to a file:
```bash
python3 chou_fasman_predictor.py > output.txt
```

---

## Web Application (GitHub Pages)
To open the website:
```
https://evanimenon.github.io/secondary-structure-predictor/
```

---

## License

Released under the **MIT License**.
Free for academic, research, and personal use.

---

## Author

**Evani Menon**
B.Tech – IIIT Delhi

---

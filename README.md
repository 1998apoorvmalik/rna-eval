# 🧬 RNA Evaluation Toolkit

This repository provides Python scripts and utilities for analyzing RNA structure and alignment data. It includes tools for evaluating secondary structure similarity, ensemble defect from base pairing probabilities, and sequence identity in multiple sequence alignments. Supported input formats include dot-bracket notation, CT files, BPP matrices, and FASTA-based MSAs.

The repository also includes a script (`ct2db.py`) to convert CT files into dot-bracket structures with up to four levels of pseudoknots. It uses a Nussinov-style dynamic programming algorithm to optimally (maximize base pairs on each page/level) assign base pairs to different bracket types: `()`, `[]`, `{}`, and `<>`.

---

## 🛠️ Requirements

- Python 3.7+
- `hopcroftkarp` (install via pip):

```bash
pip install hopcroftkarp
```

---

## 📂 Repository Structure

- `utility.py` – Core library with parsers and evaluation functions
- `eval_sd.py` – Structural distance, precision, sensitivity, F1 score
- `eval_ed.py` – Ensemble defect evaluation from BPP matrices
- `eval_seq_identity.py` – Pairwise sequence identity among aligned RNA sequences
- `ct2db.py` – Extract dot-bracket structure from `.ct` file
- `ct2seq.py` – Extract sequence from `.ct` file

---

## Accepted Input Formats

The tools in this repository accept the following types of input formats across various scripts:

**Dot-bracket string**
```
(((..(((...)))..)))
```

**CT File (first few lines)**
```
19 example
1 A 0 2 0 1
2 U 1 3 0 2
3 G 2 4 19 3
...
```

**BPP File**
```
1 19 0.88
2 18 0.84
3 17 0.91
...
```

**FASTA File (MSA)**
```
>seq1
AUGC-AUCGU
>seq2
AUGC-AA-GU
```

---

## 🔧 Environment Setup

You can quickly set up and activate the required Python environment using the provided script:

```bash
source setup_env.sh
```

This will:
- Create a virtual environment in `venv/` if it doesn't exist
- Upgrade `pip`
- Install the required dependencies (e.g., `hopcroftkarp`)
- Activate the environment

Once setup is complete, you can activate it again in the future with:

```bash
source venv/bin/activate
```

---

### 📌 Usage Instructions

#### 1. Evaluate Structural Similarity

Compare two secondary structures (dot-bracket or CT files).

```bash
./eval_sd.py struc1 struc2
```

- **Options:**
  - `--slip` / `-s` : Allow 1-nt slip in paired positions

- **Accepted formats:**
  - Dot-bracket strings: `"((..))"`
  - CT files: `*.ct`

- **Examples:**

```bash
./eval_sd.py "((..))." "((..))."
./eval_sd.py example1.ct example2.ct
./eval_sd.py example1.ct "((..))." --slip
```

---

#### 2. Evaluate Ensemble Defect

Compare a base pairing probability matrix to a target structure to evaluate the quality of a RNA ensemble.

```bash
./eval_ed.py target_structure bpp_file
```

- **Options:**
  - `--zero-index` / `-z`: If the BPP file uses 0-based indices
  - `--threshold` / `-t`: Minimum probability to consider (default: 0.01)

- **Accepted structure formats:**
  - Dot-bracket string: `"((..))"`
  - CT file: `target.ct`

- **BPP file format:**

```
i j prob   # where i < j, 1-based or 0-based (controlled with --zero-index)
```

- **Examples:**

```bash
./eval_ed.py "((..))" bpp.txt
./eval_ed.py target.ct bpp.txt -z -t 0.05
```

---

#### 3. Compute Sequence Identity

Compare aligned sequences and compute pairwise % identity.

```bash
./eval_seq_identity.py input1 input2 ... [options]
```

- **How to provide inputs:**
  - Directly via command-line:

  ```bash
  ./eval_seq_identity.py "AUGC-AUCGU" "AUGC-AA-GU"
  ```

  - As files:

  ```bash
  ./eval_seq_identity.py seq1.fasta seq2.fasta
  ```

  - From a file containing multiple sequences (MSA):

  ```bash
  ./eval_seq_identity.py msa.fasta
  ```

- **Flags:**
  - `--avg` or `-a` : Print only average sequence identity
  - `--mat` or `-m` : Print only pairwise matrix
  - Default (no flags): Prints both

- **Examples:**

```bash
./eval_seq_identity.py "AUGC-AUCGU" "AUGC-AA-GU"
./eval_seq_identity.py seq1.txt seq2.txt -a
./eval_seq_identity.py msa.fasta -m
```

---

#### 4. CT to Dot-Bracket Structure (with Pseudoknots)

Generate dot-bracket structure using base-pairing information from a .ct file.

```bash
./ct2db.py input.ct
```

- **Option:**
  - `--seq` / `-s` : Also print sequence

- **Example:**

```bash
./ct2db.py example.ct -s
```

When converting a CT file to dot-bracket format using `ct2db.py`, pseudoknots are supported via the `get_db_structure` function in `db_lib.py`.

### How Pseudoknots Are Handled

The `get_db_structure` function assigns **descending priority** to bracket types:  
`()` > `[]` > `{}` > `<>`

This means:
- It first tries to fit as many base pairs as possible into the primary (most nested) structure using `()`.
- Then it tries to assign the remaining pairs to the next available page using `[]`, and so on.

For each page (i.e. bracket level), the algorithm:
1. Selects a **maximum subset of non-crossing base pairs** from the remaining pool. This is done using a dynamic programming approach similar to the Nussinov algorithm.
2. Assigns those base pairs a unique bracket symbol depending on the current page level.
3. Removes those pairs from the pool and repeats the process with the next available bracket type.

This continues until all four pages are filled or no more pairs can be added.

### Limitations

- Only up to 4 pseudoknot levels are supported. Any additional pairs beyond those will be excluded with a warning.

### Example with Pseudoknot Levels

For the input file `./test_data/pk_struc_01.ct`, running:

```bash
python3 ct2db.py ./test_data/pk_struc_01.ct -s
```

Produces:

```
.((..[[[.))({]])}]
AUGCAGCAUGACGUAGCU
```
---

#### 5. CT to Sequence

Extract RNA sequence from a .ct file.

```bash
./ct2seq.py input.ct
```
---

## ⚙️ Utility Functions in `utility.py`

Import and use them in your Python code:

```python
from utility import (
    parse_ct_file,
    parse_secondary_structure,
    evaluate,
    calculate_pair_sequence_identity,
    calculate_msa_seq_identity,
    get_seq_from_file,
    get_ensemble_defect,
)
```

---

## 📜 License

MIT License

---

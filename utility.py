from hopcroftkarp import HopcroftKarp
from collections import defaultdict

VALID_PAIRS = set(["AU", "UA", "CG", "GC", "GU", "UG"])


def parse_secondary_structure(struc, get_pair_info=False):
    """
    Parse an RNA secondary structure string and return the indices of paired and unpaired nucleotides.

    Args:
        struc (str): A string representing RNA secondary structure, where '(' and ')' denote paired nucleotides
                    and '.' denotes unpaired nucleotides.
        get_pair_info (bool): If True, return the type of pairing (e.g., ('(', ')')) along with the indices.
                            If False, only return the indices of paired nucleotides.
    Returns:
        tuple: A tuple containing two lists:
            - A list of tuples representing the indices of paired nucleotides in the structure string.
            - A list of indices representing the indices of unpaired nucleotides in the structure string.

    Example:
        >>> parse_secondary_structure('((..((...)).))')
        ([(0, 11), (4, 9), (5, 8)], [2, 3, 6, 7, 12, 13])

    If the input string contains unbalanced parentheses, the function returns None and prints an error message.
    """
    stack1, stack2, stack3, stack4 = [], [], [], []
    paired_indices = []  # list of tuples: [(i1, j1), (i2, j2), ...]
    unpaired_indices = []  # list of indices: [i1, i2, ...]

    try:
        for i, x in enumerate(struc):
            if x == "(":
                stack1.append(i)
            elif x == "[":
                stack2.append(i)
            elif x == "{":
                stack3.append(i)
            elif x == "<":
                stack4.append(i)
            elif x == ")":
                u, v = stack1.pop(), i
                assert struc[u] == "(" and struc[v] == ")", "Invalid pair"
                paired_indices.append((u, v, ("(", ")")) if get_pair_info else (u, v))
            elif x == "]":
                # paired_indices.append((stack2.pop(), i, ('[', ']')) if get_pair_info else (stack2.pop(), i))
                u, v = stack2.pop(), i
                assert struc[u] == "[" and struc[v] == "]", "Invalid pair"
                paired_indices.append((u, v, ("[", "]")) if get_pair_info else (u, v))
            elif x == "}":
                # paired_indices.append((stack3.pop(), i, ('{', '}')) if get_pair_info else (stack3.pop(), i))
                u, v = stack3.pop(), i
                assert struc[u] == "{" and struc[v] == "}", "Invalid pair"
                paired_indices.append((u, v, ("{", "}")) if get_pair_info else (u, v))
            elif x == ">":
                # paired_indices.append((stack4.pop(), i, ('<', '>')) if get_pair_info else (stack4.pop(), i))
                u, v = stack4.pop(), i
                assert struc[u] == "<" and struc[v] == ">", "Invalid pair"
                paired_indices.append((u, v, ("<", ">")) if get_pair_info else (u, v))
            elif x == ".":
                unpaired_indices.append(i)
    except Exception as _:
        print("[Error] Unbalanced parenthesis in structure string")
        return None

    if stack1 or stack2 or stack3 or stack4:
        print("[Error] Unbalanced parenthesis in structure string")

    return paired_indices, unpaired_indices


def evaluate(
    struc1="",
    struc2="",
    struc1_paired_pos_tuple=None,
    struc1_unpaired_pos=None,
    struc2_paired_pos_tuple=None,
    struc2_unpaired_pos=None,
    allow_slip=False,
):
    """
    Evaluates one RNA secondary structure against another.

    Args:
        struc1 (str): A string representing the first RNA secondary structure.
        struc2 (str): A string representing the second RNA secondary structure.
        struc1_paired_pos_tuple (set of tuple): Optional set of paired positions for struc1.
        struc1_unpaired_pos (set of int): Optional set of unpaired positions for struc1.
        struc2_paired_pos_tuple (set of tuple): Optional set of paired positions for struc2.
        struc2_unpaired_pos (set of int): Optional set of unpaired positions for struc2.
        allow_slip (bool): Whether to allow one-nucleotide slips in pairing comparison.

    Returns:
        tuple: (precision, sensitivity, f1, structural_distance)
    """

    if struc1_paired_pos_tuple is None or struc1_unpaired_pos is None:
        res = parse_secondary_structure(struc1)
        if res is None:
            return 0, 0, 0, float("inf")
        struc1_paired_pos_tuple, struc1_unpaired_pos = res
    if struc2_paired_pos_tuple is None or struc2_unpaired_pos is None:
        res = parse_secondary_structure(struc2)
        if res is None:
            return 0, 0, 0, float("inf")
        struc2_paired_pos_tuple, struc2_unpaired_pos = res

    struc1_len = (
        2 * len(struc1_paired_pos_tuple) + len(struc1_unpaired_pos)
        if struc1 == ""
        else len(struc1)
    )

    struc2_len = (
        2 * len(struc2_paired_pos_tuple) + len(struc2_unpaired_pos)
        if struc2 == ""
        else len(struc2)
    )
    assert (
        struc1_len == struc2_len
    ), f"Length of structures must be equal\nstruc1 Length: {struc1_len}\nstruc2 Length: {struc2_len}"

    if allow_slip:
        graph = defaultdict(set)
        for i, j in struc1_paired_pos_tuple:
            for x, y in [(i, j), (i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                if (x, y) in struc2_paired_pos_tuple:
                    graph[(i, j)].add((str(x), str(y)))

        matching = HopcroftKarp(graph).maximum_matching()
        # only select values that are tuple of string
        common_paired = set(
            [(int(i), int(j)) for (i, j) in matching.values() if isinstance(i, str)]
        )

        all_paired_pos = set()
        new_unpaired_pos = []
        for i, j in common_paired:
            all_paired_pos.update([i, j])
        for i, j in struc1_paired_pos_tuple:
            if (i, j) not in matching:
                all_paired_pos.update([i, j])

        # get new unpaired pos
        for i in range(struc1_len):
            # if i not in all_paired_pos and struc1[i] != "*": [NOTE]: Removed this condition
            if i not in all_paired_pos:
                new_unpaired_pos.append(i)
        # get common unpaired pos
        common_unpaired = set(new_unpaired_pos).intersection(struc2_unpaired_pos)
    else:
        common_paired = set(struc1_paired_pos_tuple).intersection(
            struc2_paired_pos_tuple
        )
        common_unpaired = set(struc1_unpaired_pos).intersection(struc2_unpaired_pos)

    precision = len(common_paired) / (len(struc1_paired_pos_tuple) + 1e-10)
    sensitivity = len(common_paired) / (len(struc2_paired_pos_tuple) + 1e-10)
    f1 = 2 * precision * sensitivity / (precision + sensitivity + 1e-10)
    structural_distance = struc2_len - (2 * len(common_paired) + len(common_unpaired))
    return precision, sensitivity, f1, structural_distance


def get_alignment_to_sequence_mapping(aligned_sequence):
    """
    Returns a dictionary mapping the indices of the aligned sequence to the indices of the unaligned sequence.

    Args:
        aligned_sequence (str): A string representing the aligned sequence, where '-' denotes a gap.

    Returns:
        dict: A dictionary mapping the indices of the aligned sequence to the indices of the unaligned sequence.

    Example:
        >>> get_alignment_to_sequence_mapping('AUCG-AUCG')
        {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 5, 7: 6}
    """
    mapping = {}
    j = 0
    for i, x in enumerate(aligned_sequence):
        if x != "-":
            mapping[i] = j
            j += 1
        else:
            mapping[i] = j - 1
    return mapping


def pairs_to_struc(pairs, seq_length):  # need to rewrite this function
    struc = ["." for _ in range(seq_length)]
    for i, j in pairs:
        struc[i] = "("
        struc[j] = ")"
    return "".join(struc)


def map_consns_struc_to_aln_seq(consns_struc, aln_seq):
    a2s_map = get_alignment_to_sequence_mapping(aln_seq)
    seq = aln_seq.replace("-", "")  # ungapped sequence

    # get the structure corresponding to the unaligned sequence
    struc = ["."] * len(seq)

    for i, j, (b1, b2) in parse_secondary_structure(consns_struc, True)[0]:
        if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
            continue
        if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
            continue

        struc[a2s_map[i]] = b1
        struc[a2s_map[j]] = b2

    return "".join(struc), seq


def map_consns_bpp_to_aln_seq(consns_bpp, aln_seq, threshold=0.001):
    a2s_map = get_alignment_to_sequence_mapping(aln_seq)
    seq = aln_seq.replace("-", "")  # ungapped sequence

    # get the structure corresponding to the unaligned sequence
    bpp = defaultdict(lambda: defaultdict(float))

    for i in range(len(aln_seq)):
        for j in consns_bpp[i].keys():
            if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
                continue
            if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
                continue
            if consns_bpp[i][j] >= threshold:
                bpp[a2s_map[i]][a2s_map[j]] = consns_bpp[i][j]

    return bpp, seq


def get_struc_from_file(file_path, backward_search=False):
    """
    Get the structure string from a file.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
        for line in lines[::-1] if backward_search else lines:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] in [".", "("]:
                return line


def get_seq_from_file(file_path, backward_search=False):
    """
    Get the sequence string from a file.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
        for i, line in enumerate((lines[::-1] if backward_search else lines)):
            header = lines[i - 1].strip() if i > 0 else ""
            line = line.strip()
            if line[0] in set(["A", "U", "C", "G", "T", "-"]):
                return header, line


def get_bpp_from_file(file_path, zero_index=False, threshold=0.01):
    """
    Get the bpp matrix from a file.
    """
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File {file_path} not found.")
    except Exception as e:
        raise Exception(f"Error reading file {file_path}: {e}")
    if len(lines) == 0:
        raise ValueError(f"File {file_path} is empty.")

    bpp_matrix = defaultdict(lambda: defaultdict(float))
    seq_length = -1
    for line in lines:
        split = line.split(" ")
        if len(split) != 3:
            raise ValueError(
                f"Invalid BPP format in line: {line.strip()}. Expected format: 'i j prob'."
            )
        elif len(split[0]) < 100 and float(split[2]) >= threshold:
            i, j, prob = int(split[0]), int(split[1]), float(split[2])
            if i > j:
                i, j = j, i
                print(
                    "[Warning] BPP file contains i > j. Swapping indices i={i} and j={j}."
                )
            elif i == j:
                raise ValueError("BPP file contains self-pairing.")
            if not zero_index:
                i -= 1
                j -= 1
            bpp_matrix[i][j] = prob
            seq_length = max(seq_length, j + 1)
    return bpp_matrix, seq_length


def convert_bpp_to_list(bpp_matrix, seq_length, threshold=0.01):
    consns_bpp = []
    for i in range(seq_length):
        for j in bpp_matrix[i].keys():
            if bpp_matrix[i][j] >= threshold:
                consns_bpp.append((i, j, bpp_matrix[i][j]))
    return consns_bpp


def get_unpaired_probs(paired_probs, seq_length):
    unpaired_probs = {}

    for i in range(seq_length):
        unpaired_probs[i] = 1.00

    for l in paired_probs.keys():
        for r, p in paired_probs[l].items():
            unpaired_probs[l] -= p
            unpaired_probs[r] -= p

    return unpaired_probs


def parse_ct_file(ct_file_path):
    """
    Parse a CT file and return paired and unpaired positions.
    """
    paired_pos_tuple = set()
    paired_pos = set()
    unpaired_pos = set()
    seq = []

    with open(ct_file_path, "r") as file:
        # Skip header lines if necessary
        next(file)

        for line in file.readlines():
            parts = line.strip().split()
            if len(parts) < 6:
                continue  # skip malformed lines

            # get the nucleotide
            nuc = parts[1]
            seq.append(nuc)

            # get the pairing information
            i = int(parts[0]) - 1  # Adjust index for 0-based Python indexing
            j = int(parts[4]) - 1  # Adjust index for 0-based Python indexing

            if j != -1 and i < j:
                if (
                    i in paired_pos
                    or j in paired_pos
                    or i in unpaired_pos
                    or j in unpaired_pos
                ):
                    raise ValueError("CT file contains overlapping pairs.")
                paired_pos_tuple.add((i, j))
                paired_pos.add(i)
                paired_pos.add(j)
            elif j == -1:
                if i in paired_pos:
                    raise ValueError("CT file contains overlapping pairs.")
                unpaired_pos.add(i)

    return "".join(seq), paired_pos_tuple, unpaired_pos


def calculate_pair_sequence_identity(seq1, seq2):
    """
    Calculate the pairwise sequence identity between two aligned sequences.
    Skips positions where both sequences have gaps.

    Args:
    seq1 (str): First aligned sequence with possible gaps ('-')
    seq2 (str): Second aligned sequence with possible gaps ('-')

    Returns:
    float: Pairwise sequence identity as a percentage
    """

    # Ensure the sequences are of the same length
    if len(seq1) != len(seq2):
        raise ValueError(
            "The two sequences must be of the same length for comparison, but found {} {} {}".format(
                seq1, "\n\n", seq2
            )
        )

    # Initialize counters
    identical_positions = 0
    aligned_positions = 0

    # Loop through each position to compare the two sequences
    for base1, base2 in zip(seq1, seq2):
        # Skip positions where both sequences have gaps
        if base1 == "-" and base2 == "-":
            continue

        # Increment the counter for aligned positions
        aligned_positions += 1

        # Count identical positions
        if base1 == base2:
            identical_positions += 1

    # If no aligned positions are present, return 0 to avoid division by zero
    if aligned_positions == 0:
        return 0.0

    # Calculate the pairwise sequence identity as a percentage
    return (identical_positions / aligned_positions) * 100


def calculate_msa_seq_identity(msa):
    """
    Calculate the average pairwise sequence identity for a multiple sequence alignment (MSA).

    Args:
    msa (list of str): List containing aligned sequences.

    Returns:
    float: Average pairwise sequence identity.
    """
    num_sequences = len(msa)
    if num_sequences < 2:
        raise ValueError(
            "At least two sequences are required to compute pairwise sequence identity"
        )

    total_identity = 0
    num_comparisons = 0

    # Calculate pairwise identity for each unique pair of sequences
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            identity = calculate_pair_sequence_identity(msa[i], msa[j])
            total_identity += identity
            num_comparisons += 1

    # Calculate the average identity
    average_identity = total_identity / num_comparisons if num_comparisons > 0 else 0
    return average_identity


# ------------------------------ Ensemble defect calculations ------------------------------
def get_unpaired_probs(pair_probs, length):
    unpaired_probs = {}
    for i in range(length):
        unpaired_probs[i] = 1.00
    for l in pair_probs.keys():
        for r, p in pair_probs[l].items():
            unpaired_probs[l] -= p
            unpaired_probs[r] -= p
    return unpaired_probs


def compute_ensemble_defect(pairs, unpairs, length, pair_probs):
    """
    Compute the ensemble defect based on paired and unpaired positions.
    Args:
        pairs (list): List of paired positions.
        unpairs (list): List of unpaired positions.
        length (int): Length of the sequence.
        pair_probs (dict): Dictionary of pair probabilities.
        unpair_probs (dict): Dictionary of unpair probabilities.
    Returns:
        float: Ensemble defect value.
    """
    unpair_probs = get_unpaired_probs(pair_probs, length)
    value = length
    for pair in pairs:
        value -= 2 * pair_probs[pair[0]][pair[1]]
    for i in unpairs:
        value -= unpair_probs[i]
    return value


def get_ensemble_defect(struc_data, bpp_file, zero_index=False, threshold=0.01):
    """
    Calculate the ensemble defect of a predicted BPP matrix against a target structure.

    Args:
        struc (str): Target structure (dot-bracket or CT file).
        bpp (str): Predicted BPP file.
        threshold (float): Probability threshold for BPPs.

    Returns:
        float: Ensemble defect value.
    """
    try:
        pairs, unpairs = None, None
        if struc_data.endswith(".ct"):
            _, pairs, unpairs = parse_ct_file(struc_data)
        else:
            pairs, unpairs = parse_secondary_structure(struc_data) 
    except Exception as e:
        print(f"[Error] Failed to parse structure file: {e}")
        return float('inf')
    try:
        bpp, _ = get_bpp_from_file(bpp_file, zero_index, threshold)
    except Exception as e:
        print(f"[Error] Failed to parse BPP file: {e}")
        return float('inf')
    return compute_ensemble_defect(
        pairs,
        unpairs,
        2 * len(pairs) + len(unpairs),
        bpp,
    )

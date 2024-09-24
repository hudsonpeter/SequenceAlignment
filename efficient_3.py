import psutil
import sys
import time

delta = 30
alpha = {
    "A": {"A": 0, "C": 110, "G": 48, "T": 94},
    "C": {"A": 110, "C": 0, "G": 118, "T": 48},
    "G": {"A": 48, "C": 118, "G": 0, "T": 110},
    "T": {"A": 94, "C": 48, "G": 110, "T": 0},
}


def generate_string(base_string, insertions):
    generated_string = base_string
    for index in insertions:
        generated_string = (
            base_string[: index + 1] + base_string + base_string[index + 1 :]
        )
        base_string = generated_string
    return generated_string


def read_input_file(filename):
    base_str1, base_str2 = None, None
    str1_insertions, str2_insertions = [], []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                if line.isdigit():
                    if base_str2 is not None:
                        str2_insertions.append(int(line))
                    else:
                        str1_insertions.append(int(line))
                else:
                    if base_str1 is None:
                        base_str1 = line
                    else:
                        base_str2 = line

    return base_str1, str1_insertions, base_str2, str2_insertions


def space_efficient_alignment(gene_seq1, gene_seq2):
    gene_seq1_len = len(gene_seq1)
    gene_seq2_len = len(gene_seq2)
    dp = []
    for i in range(2):
        dp.append([0] * (gene_seq2_len + 1))
    for i in range(gene_seq2_len + 1):
        dp[0][i] = delta * i
    for i in range(1, gene_seq1_len + 1):
        dp[1][0] = i * delta
        for j in range(1, gene_seq2_len + 1):
            dp[1][j] = min(
                dp[0][j - 1] + alpha[gene_seq1[i - 1]][gene_seq2[j - 1]],
                dp[0][j] + delta,
                dp[1][j - 1] + delta,
            )
        for j in range(gene_seq2_len + 1):
            dp[0][j] = dp[1][j]

    penalty = dp[1]
    return penalty


def get_minimum_penalty(gene_seq1, gene_seq2):

    gene_seq1_len, gene_seq2_len = len(gene_seq1), len(gene_seq2)

    dp = [[0] * (gene_seq2_len + 1) for _ in range(gene_seq1_len + 1)]

    for i in range(gene_seq2_len + 1):
        dp[0][i] = i * delta
    for i in range(gene_seq1_len + 1):
        dp[i][0] = i * delta
    for i in range(1, gene_seq1_len + 1):
        for j in range(1, gene_seq2_len + 1):
            dp[i][j] = min(
                dp[i - 1][j - 1] + alpha[gene_seq1[i - 1]][gene_seq2[j - 1]],
                dp[i][j - 1] + delta,
                dp[i - 1][j] + delta,
            )

    gene1_aligned, gene2_aligned = "", ""

    while gene_seq1_len and gene_seq2_len:
        if (
            dp[gene_seq1_len][gene_seq2_len]
            == dp[gene_seq1_len - 1][gene_seq2_len - 1]
            + alpha[gene_seq1[gene_seq1_len - 1]][gene_seq2[gene_seq2_len - 1]]
        ):
            gene1_aligned = gene_seq1[gene_seq1_len - 1] + gene1_aligned
            gene2_aligned = gene_seq2[gene_seq2_len - 1] + gene2_aligned
            gene_seq1_len -= 1
            gene_seq2_len -= 1
        elif (
            dp[gene_seq1_len][gene_seq2_len]
            == dp[gene_seq1_len - 1][gene_seq2_len] + delta
        ):
            gene1_aligned = gene_seq1[gene_seq1_len - 1] + gene1_aligned
            gene2_aligned = "_" + gene2_aligned
            gene_seq1_len -= 1
        elif (
            dp[gene_seq1_len][gene_seq2_len]
            == dp[gene_seq1_len][gene_seq2_len - 1] + delta
        ):
            gene1_aligned = "_" + gene1_aligned
            gene2_aligned = gene_seq2[gene_seq2_len - 1] + gene2_aligned
            gene_seq2_len -= 1
    while gene_seq1_len:
        gene1_aligned = gene_seq1[gene_seq1_len - 1] + gene1_aligned
        gene2_aligned = "_" + gene2_aligned
        gene_seq1_len -= 1
    while gene_seq2_len:
        gene1_aligned = "_" + gene1_aligned
        gene2_aligned = gene_seq2[gene_seq2_len - 1] + gene2_aligned
        gene_seq2_len -= 1
    return [gene1_aligned, gene2_aligned, dp[len(gene_seq1)][len(gene_seq2)]]


def divide_and_conquer(gene_seq1, gene_seq2):
    gene_seq1_len = len(gene_seq1)
    gene_seq2_len = len(gene_seq2)
    gene_seq1_mid = gene_seq1_len // 2
    if gene_seq1_len < 2 or gene_seq2_len < 2:
        return get_minimum_penalty(gene_seq1, gene_seq2)
    else:
        left_half = space_efficient_alignment(gene_seq1[:gene_seq1_mid], gene_seq2)
        right_half = space_efficient_alignment(
            gene_seq1[gene_seq1_mid:][::-1], gene_seq2[::-1]
        )
        alignment_costs = [
            left_half[j] + right_half[gene_seq2_len - j]
            for j in range(gene_seq2_len + 1)
        ]
        # print("alignment_costs:", alignment_costs)
        optimal_index = alignment_costs.index(min(alignment_costs))
        # print("optimal_index:", optimal_index)
        callLeft = divide_and_conquer(
            gene_seq1[:gene_seq1_mid], gene_seq2[:optimal_index]
        )
        callRight = divide_and_conquer(
            gene_seq1[gene_seq1_mid:], gene_seq2[optimal_index:]
        )
        result = [callLeft[r] + callRight[r] for r in range(3)]
        # print(result)
    return result[0], result[1], result[2]


def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        return

    input_file = sys.argv[1]
    output_filename = sys.argv[2]
    base_str1, base_str1_indexes, base_str2, base2_str_indexes = read_input_file(
        input_file
    )

    gene_seq1 = generate_string(base_str1, base_str1_indexes)
    gene_Seq2 = generate_string(base_str2, base2_str_indexes)

    start_time = time.time()
    str1_aligned, str2_aligned, penalty = divide_and_conquer(gene_seq1, gene_Seq2)
    end_time = time.time()
    time_taken = (end_time - start_time) * 1000

    process = psutil.Process()
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss / 1024)

    with open(output_filename, "w") as output_file:
        output_file.write(str(penalty) + "\n")
        output_file.write(str1_aligned + "\n")
        output_file.write(str2_aligned + "\n")
        output_file.write(str(time_taken) + "\n")
        output_file.write(str(memory_consumed) + "\n")
        output_file.close()


if __name__ == "__main__":
    main()

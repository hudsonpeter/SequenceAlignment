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


def get_minimum_penalty(gene1: str, gene2: str):

    process = psutil.Process()
    gene1_ptr, gene2_ptr = 0, 0
    gene1_len, gene2_len = len(gene1), len(gene2)

    dp = [[0] * (gene2_len + 1) for _ in range(gene1_len + 1)]
    for i in range(gene1_len + 1):
        dp[i][0] = i * delta
    for j in range(gene2_len + 1):
        dp[0][j] = j * delta

    gene1_ptr = 1
    while gene1_ptr <= gene1_len:
        gene2_ptr = 1
        while gene2_ptr <= gene2_len:
            if gene1[gene1_ptr - 1] == gene2[gene2_ptr - 1]:
                dp[gene1_ptr][gene2_ptr] = dp[gene1_ptr - 1][gene2_ptr - 1]
            else:
                dp[gene1_ptr][gene2_ptr] = min(
                    dp[gene1_ptr - 1][gene2_ptr - 1]
                    + alpha[gene1[gene1_ptr - 1]][gene2[gene2_ptr - 1]],
                    dp[gene1_ptr - 1][gene2_ptr] + delta,
                    dp[gene1_ptr][gene2_ptr - 1] + delta,
                )
            gene2_ptr += 1
        gene1_ptr += 1

    max_gene_len = gene2_len + gene1_len
    gene1_ptr = gene1_len
    gene2_ptr = gene2_len

    gene1_pos = max_gene_len
    gene2_pos = max_gene_len

    gene1_matched = [0] * (max_gene_len + 1)
    gene2_matched = [0] * (max_gene_len + 1)

    while not (gene1_ptr == 0 or gene2_ptr == 0):
        if gene1[gene1_ptr - 1] == gene2[gene2_ptr - 1]:
            gene1_matched[gene1_pos] = gene1[gene1_ptr - 1]
            gene2_matched[gene2_pos] = gene2[gene2_ptr - 1]
            gene1_pos -= 1
            gene2_pos -= 1
            gene1_ptr -= 1
            gene2_ptr -= 1
        elif (
            dp[gene1_ptr - 1][gene2_ptr - 1]
            + alpha[gene1[gene1_ptr - 1]][gene2[gene2_ptr - 1]]
        ) == dp[gene1_ptr][gene2_ptr]:

            gene1_matched[gene1_pos] = gene1[gene1_ptr - 1]
            gene2_matched[gene2_pos] = gene2[gene2_ptr - 1]
            gene1_pos -= 1
            gene2_pos -= 1
            gene1_ptr -= 1
            gene2_ptr -= 1

        elif (dp[gene1_ptr - 1][gene2_ptr] + delta) == dp[gene1_ptr][gene2_ptr]:
            gene1_matched[gene1_pos] = gene1[gene1_ptr - 1]
            gene2_matched[gene2_pos] = "_"
            gene1_pos -= 1
            gene2_pos -= 1
            gene1_ptr -= 1

        elif (dp[gene1_ptr][gene2_ptr - 1] + delta) == dp[gene1_ptr][gene2_ptr]:
            gene1_matched[gene1_pos] = "_"
            gene2_matched[gene2_pos] = gene2[gene2_ptr - 1]
            gene1_pos -= 1
            gene2_pos -= 1
            gene2_ptr -= 1

    while gene1_pos > 0:
        if gene1_ptr > 0:
            gene1_ptr -= 1
            gene1_matched[gene1_pos] = gene1[gene1_ptr]
            gene1_pos -= 1
        else:
            gene1_matched[gene1_pos] = "_"
            gene1_pos -= 1

    while gene2_pos > 0:
        if gene2_ptr > 0:
            gene2_ptr -= 1
            gene2_matched[gene2_pos] = gene2[gene2_ptr]
            gene2_pos -= 1
        else:
            gene2_matched[gene2_pos] = "_"
            gene2_pos -= 1

    index = 1
    gene1_ptr = max_gene_len
    while gene1_ptr >= 1:
        if (gene2_matched[gene1_ptr] == "_") and gene1_matched[gene1_ptr] == "_":
            index = gene1_ptr + 1
            break

        gene1_ptr -= 1

    min_penalty = dp[gene1_len][gene2_len]

    gene1_ptr = index
    gene1_aligned = ""
    while gene1_ptr <= max_gene_len:
        gene1_aligned += gene1_matched[gene1_ptr]
        gene1_ptr += 1

    gene1_ptr = index
    gene2_aligned = ""
    while gene1_ptr <= max_gene_len:
        gene2_aligned += gene2_matched[gene1_ptr]
        gene1_ptr += 1

    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss / 1024)
    return min_penalty, gene1_aligned, gene2_aligned, memory_consumed


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
    penalty, str1_aligned, str2_aligned, memory_used = get_minimum_penalty(
        gene_seq1, gene_Seq2
    )
    end_time = time.time()
    time_taken = (end_time - start_time) * 1000

    with open(output_filename, "w") as output_file:
        output_file.write(str(penalty) + "\n")
        output_file.write(str1_aligned + "\n")
        output_file.write(str2_aligned + "\n")
        output_file.write(str(time_taken) + "\n")
        output_file.write(str(memory_used) + "\n")
        output_file.close()


if __name__ == "__main__":
    main()

import sys

def generate_string(base_string, insertions):
    for index in insertions:
        # print(index)
        generated_string = base_string[:index+1] + base_string + base_string[index+1:]
        # print(generated_string)
        base_string = generated_string
    return generated_string

def read_input_file(filename):
    base_str1 = None
    base_str1_indexes = []
    base_str2 = None
    base2_str_indexes = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                if line.isdigit():
                    if base_str2 is not None:
                        base2_str_indexes.append(int(line))
                    else:
                        base_str1_indexes.append(int(line))
                else:
                    if base_str1 is None:
                        base_str1 = line
                    else:
                        base_str2 = line

    return base_str1, base_str1_indexes, base_str2, base2_str_indexes

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file")
        return

    input_file = sys.argv[1]
    base_str1, base_str1_indexes, base_str2, base2_str_indexes = read_input_file(input_file)
    print("Base string 1:", base_str1)
    print("Base str1 indexes:", base_str1_indexes)
    print("Base string 2", base_str2)
    print("Base str2 indexes", base2_str_indexes)
    
    base1_generated = generate_string(base_str1, base_str1_indexes)
    base2_generated = generate_string(base_str2, base2_str_indexes)

    print("generated str1: ", base1_generated)
    print("generated str2: ", base2_generated)

if __name__ == "__main__":
    main()

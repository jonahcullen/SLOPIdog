import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--input")
parser.add_argument('-o', "--output")
parser.add_argument('-r', "--region")
args = parser.parse_args()

with open(args.input, "rt") as infile:
    with open(args.output, "wt") as outfile:
        print("pos\tchr\tcM", file=outfile)
        for line in infile:
            line = line.strip().split('\t')
            if line[0] == args.region:
                print(f"{line[3]}\t{line[0]}\t{line[2]}", file=outfile)

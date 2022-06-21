# This python takes a .cov.gz file and writes a .txt file with the desired coverage.
# Usage example : python3 smart2_cov.py HFD-M16-1.cov.gz 10

import gzip
import sys

# The arguments passed in the shell are saved in variables.
command_line = sys.argv
file_gz = command_line[1]
coverage = int(command_line[2])

# Unzips the file to be able to read it. (contents) is a list in which element is a line.
file_coverage = gzip.GzipFile(file_gz, 'rb')
contents = file_coverage.readlines()

# Gives the name of the file with the coverage as a prefix, writes a .txt file.
with open(str(coverage) + "X-" + file_gz[:-6] + "txt", "w") as f:

    # Iterates over each line in the file. Decodes it, since it's in byte form. Splits it over tabulations, to be able
    # to analyze each column individually.
    for line in contents:
        line_ = line.decode().split("\t")

        # Keeping only the lines that have at least the coverage wanted. Adds 'chr' in front of the chromosome number
        # and subtracts 1 from the start position in column 2 to convert it in 0 based coordinates. Columns 5 and 6 are
        # not needed, and column 3 is divided by 100. We want 3 decimals, hence the '.3f' .
        if int(line_[4]) + int(line_[5]) >= coverage:
            f.write("chr{0}\t{1}\t{2}\t{3}\n".format(line_[0], int(line_[1]) - 1, int(line_[2]),
                                                     format((float(line_[3]) / 100), '.3f')))

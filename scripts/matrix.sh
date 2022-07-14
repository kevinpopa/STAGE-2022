#!/bin/sh
module load python/3.9.6

# If doesn't exist, makes a matrix directory. coverage_files to clear the directory after the operation
mkdir -p matrix
mkdir -p coverage_files

# Variable that contains the coverage wanted
COVERAGE=10

# Loop that looks for all the .cov.gz files. Runs the smart2_cov.py script on it with the desired coverage.
# Then sorts the file, changes the extension and moves it in the ./matrix directory.
for file in *.cov.gz; do
  python3 ~/bin/smart2_cov.py "$file" $COVERAGE
  sort -k1,1 -V -s ${COVERAGE}*.txt > "${file%??????}"bg
  mv "${file%??????}"bg matrix
  mv ${COVERAGE}*.txt coverage_files
done
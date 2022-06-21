# python3 matrix.py 10X-HFD-M16-1.txt

import os.path
import sys
import pandas as pd

command_line = sys.argv
input1 = command_line[1]


if os.path.exists("matrix.txt"):
    df1 = pd.read_csv(input1, sep='\t', lineterminator='\n', names=["chrom", "start", "end", input1[:-4]],
                      dtype={'chrom': 'str'})
    df1 = df1.set_index(['chrom', 'start', "end"])

    df2 = pd.read_csv("matrix.txt", sep='\t', lineterminator='\n', dtype={'chrom': 'str'})
    df2 = df2.set_index(['chrom', 'start', "end"])

    df_3 = pd.merge(df2, df1, how="outer", on=["chrom", "start", "end"])
    df_3 = df_3.sort_values(by=['chrom', 'start', "end"])
    df_3 = df_3.fillna(value=-1)
    df_3.to_csv("matrix.txt", sep="\t", float_format='%.3f')
else:
    df = pd.read_csv(input1, sep='\t', lineterminator='\n', names=["chrom", "start", "end", input1[:-4]],
                     dtype={'chrom': 'str'})
    df = df.set_index(['chrom', 'start', "end"])
    df = df.sort_values(by=['chrom', 'start', "end"])
    df.to_csv("matrix.txt", sep="\t", float_format='%.3f')

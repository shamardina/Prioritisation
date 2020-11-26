#!/usr/bin/python3


import pandas as pd
import sys
import os.path


INPUT = sys.argv[1]
data = pd.read_excel(INPUT, comment="#", na_filter=False)


data.to_csv(os.path.splitext(INPUT)[0] + ".csv", sep="\t", index=False)

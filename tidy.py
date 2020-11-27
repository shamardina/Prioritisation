#!/usr/bin/python3


import csv
import sys


# Tidy the result:
#   * split VEP_HGVSc (e.g. "ENST00000373103.1:c.340C>T") and VEP_HGVSp (e.g. "ENSP00000362195.1:p.Gln114Ter") on ":"
#   * in CLNV_CLNSIGCONF replace %3B with semicolon
#   * in HGMD_PHEN remove %2C (comma, probably an artifact)


INPUT = sys.argv[1]
SPLIT = ["VEP_HGVSc", "VEP_HGVSp"]


with open(INPUT, newline='') as f:
    fields = f.readline().rstrip().split("\t")

    new_fields = fields[:]
    for s in SPLIT:
        new_fields.insert(new_fields.index(s), "%s_ID" % s)
    print("\t".join(new_fields))

    reader = csv.DictReader(f, fieldnames=fields, delimiter="\t")
    for row in reader:
        for s in SPLIT:
            s_ID = "%s_ID" % s
            if ":" in row[s]:
                row[s_ID], row[s] = row[s].split(":")
            else:
                row[s_ID] = ""

        row["CLNV_CLNSIGCONF"] = row["CLNV_CLNSIGCONF"].replace("%3B", "; ")
        row["HGMD_PHEN"] = row["HGMD_PHEN"].replace("%2C", "")

        print("\t".join([row[i] for i in new_fields]))

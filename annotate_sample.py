#!/usr/bin/python3


import csv
import collections
import sys
import subprocess



CONF = sys.argv[1]
# from config, export variables:
VARS = ["DOMAIN", "INCLUDE", "FLAGSHIP_INCLUDE", "SEX", "ETHN", "CONSANG", "RELATED"]

# read config:
echo_line = "; ".join(["echo ${%s}" % v for v in VARS])
configuration_from_bash = subprocess.check_output(["bash", "-c", "source {}; {}".format(CONF, echo_line)], encoding="utf8").split("\n")[:-1]
PARAMS = dict(zip(VARS, configuration_from_bash))

# annotate samples with this info:
FIELDS = ["WGS_ID",
          "Sex", "Sex_karyotype", "Sex_flag",
          "Ethnicity",
          "Consanguinity", "Consanguinity_f_score",
          "Network", "Relatedness (Sample_ID:in_flagship:domain, if different)",
          "In_flagship"]

wgs_to_sample = {}
anno = {}  # with sampleID as key


def how_related(Z0, Z1, Z2, PI_HAT):
    """Based on assign_relatedness.R from Salih Tuna"""
    z0, z1, z2, pi_hat = [round(float(i)/0.25)*0.25 for i in [Z0, Z1, Z2, PI_HAT]]

    # Check relationship
    # Idenical:
    if z0 == 0 and z1 == 0 and z2 == 1 and pi_hat == 1:
        relationship = "Duplicate or MZ twin"
    # Parent:
    elif z0 == 0 and z1 == 1 and z2 == 0 and pi_hat==0.5:
        relationship = "Parent/Offspring"
    # Sibling:
    elif z0 == 0.25 and z1 == 0.5 and z2 == 0.25 and pi_hat == 0.5:
        relationship = "Full sibling"
    # Grandparent or half sibling or aunt or nephew:
    elif z0 == 0.5 and z1 == 0.5 and z2 == 0 and pi_hat == 0.25:
        relationship = "Half sibling/Grandparent/Aunt/Nephew"
    # Cousin:
    elif z0 == 0.75 and z1 == 0.25 and z2==0:
        relationship = "First cousin"
    # Unknown:
    else:
        relationship = "Unknown"

    return relationship


def related_domain_annotation(sample, domain, relatedness):
    """Create a string for the summary annotation of related samples"""
    if sample in relatedness:
        if relatedness[sample] == domain:
            rv = ":Y"
        else:
            rv = ":Y:%s" % relatedness[sample]
    else:
        rv = ""
    return rv


###### Create sample annotation table ######

# Get all samples from domain, create ID mapping, set defaults
with open(PARAMS["INCLUDE"], newline='') as f:                   # IlluminaID  wgsID  sampleID  Domain
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if row[3] == PARAMS["DOMAIN"]:
            wgs_to_sample[row[1]] = row[2]
            anno[row[2]]={"WGS_ID": row[1], "In_flagship": "N"}


# Flag samples from the flagship release
flagship_sample_domain = {}  # will be used to annotate related samples
with open(PARAMS["FLAGSHIP_INCLUDE"], newline='') as f:          # IlluminaID  wgsID  sampleID  Domain
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        flagship_sample_domain[row[2]] = row[3]
        # only if domain and wgsID are the same, we don't need to reanalise:
        if (row[3] == PARAMS["DOMAIN"]) and (row[2] in anno) and (anno[row[2]]["WGS_ID"] == row[1]):
                anno[row[2]]["In_flagship"] = "Y"


# Add genomic sex
with open(PARAMS["SEX"], newline='') as f:                       # sample(wgsID) XAutoRatio YAutoRatio Hratio declared_gender genotype_gender gender_mismatch discordant_Hratio outside_thresholds flag karyotype
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["sample"] in wgs_to_sample:
            anno[wgs_to_sample[row["sample"]]].update({"Sex": row["genotype_gender"],
                                                       "Sex_karyotype": row["karyotype"],
                                                       "Sex_flag": row["flag"]})


# Add genomic ethnicity
with open(PARAMS["ETHN"], newline='') as f:                      # Sample(sampleID) megapop
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["Sample"] in anno:
            anno[row["Sample"]]["Ethnicity"] = row["megapop"]


# Add consanguinity information
with open(PARAMS["CONSANG"], newline='') as f:                   # ID(sampleID) nsnp f
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["ID"] in anno:
            if float(row["f"]) > 0.03:
                consang_status = "Y"
            else:
                consang_status = "N"
            anno[row["ID"]].update({"Consanguinity": consang_status, "Consanguinity_f_score": row["f"]})


# Add relatedness information; for related samples indicate if they were in the flagship release
sample_network = {}
relatedness_type = {}
with open(PARAMS["RELATED"], newline='') as f:                   # Network FID1 IID1(sampleID) FID2 IID2 NA UN Z0 Z1 Z2 PI_HAT
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["IID1"] in anno or row["IID2"] in anno:
            R = how_related(row["Z0"], row["Z1"], row["Z2"], row["PI_HAT"])
            relatedness_type[(row["IID1"], row["IID2"])] = R
            relatedness_type[(row["IID2"], row["IID1"])] = R
            for i in ["IID1", "IID2"]:
                if row[i] not in sample_network:
                    sample_network[row[i]] = {"networks": set(), "samples": set()}
                sample_network[row[i]]["networks"].add(row["Network"])
            sample_network[row["IID1"]]["samples"].add(row["IID2"])
            sample_network[row["IID2"]]["samples"].add(row["IID1"])

for sample in anno:
    if sample in sample_network:
        relatedness = ", ".join(sorted(["{} ({}{})".format(relatedness_type[(sample, rel)], rel, related_domain_annotation(rel, PARAMS["DOMAIN"], flagship_sample_domain)) for rel in sample_network[sample]["samples"]]))
        network = ", ".join(sample_network[sample]["networks"])
    else:
        relatedness = "N"
        network = "NA"
    anno[sample].update({"Relatedness (Sample_ID:in_flagship:domain, if different)": relatedness,
                         "Network": network})


###### Output sample annotation table ######
print("\t".join(["SAMPLE"] + FIELDS))
for sample in anno:
    print("\t".join([sample] + [anno[sample][a] for a in FIELDS]))

#!/usr/bin/python3


import argparse
import csv
import sys
import os.path


# ./merge.py -a file1 -b file2 -1 name1 -2 name2 -c name_after


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--file1",
        required=True,
        help="provide filename for file 1",
        metavar="FILE1",
        dest="filename1")
    parser.add_argument(
        "-b",
        "--file2",
        required=True,
        help="provide filename for file 2",
        metavar="FILE2",
        dest="filename2")
    parser.add_argument(
        "-1",
        "--column1",
        required=True,
        help="provide name of the column in file 1 to join on",
        metavar="COLUMN1",
        dest="column1")
    parser.add_argument(
        "-2",
        "--column2",
        required=True,
        help="provide name of the column in file 2 to join on",
        metavar="COLUMN2",
        dest="column2")
    parser.add_argument(
        "-c",
        "--after",
        required=True,
        help="provide name of the column in file 1 after which the fields from file 2 will be inserted",
        metavar="AFTER",
        dest="after")
    parser.add_argument(
        "-p",
        "--prefix",
        help="provide prefix to add to the field names from file 2 (optional)",
        default="",
        metavar="PREFIX",
        dest="prefix")
    parser.add_argument(
        "-r",
        "--remove",
        help="remove merging column from file 1 or 2",
        choices=["1", "2"],
        metavar="REMOVE",
        dest="remove"
    )

    return parser.parse_args()


def check_args(args):
    """Check that files exit; check that fields are present in the headers"""
    for f in [args.filename1, args.filename2]:
        if not os.path.exists(f):
            sys.exit("File {} not found".format(f))

    header1 = get_header(args.filename1)
    for c in [args.column1, args.after]:
        if c not in header1:
            sys.exit("File {} doesn't have field {}".format(args.filename1, c))

    header2 = get_header(args.filename2)
    if args.column2 not in header2:
        sys.exit("File {} doesn't have field {}".format(args.filename2, args.column2))


def read_annotation(filename, column):
    """Read annotation data, output dictionary with values form column as key"""
    data = {}
    with open(filename, newline='') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            data[row[column]] = row
    return data


def get_header(filename):
    with open(filename) as f:
        header = f.readline().rstrip().split("\t")
    return header


def main():
    args = parse_arguments()
    check_args(args)
    annotation = read_annotation(args.filename2, args.column2)
    anno_header = get_header(args.filename2)
    header = get_header(args.filename1)
    i1 = header.index(args.column1)
    a1 = header.index(args.after) + 1
    if args.remove:
        if args.remove == "2":
            anno_header.remove(args.column2)
        else:    # i.e. is "1"
            header.pop(i1)
            if a1 > i1:
                a1 -= 1

    merged_header = header[0:a1] + [args.prefix + a for a in anno_header] + header[a1:]

    print("\t".join(merged_header))
    with open(args.filename1, newline='') as f:
        f.readline()
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            key = row[i1]
            if key not in annotation:
                annotation[key] = {s: "NA" for s in anno_header}
            if args.remove and args.remove == "1":
                row.pop(i1)
            print("\t".join(row[0:a1] + [annotation[key][s] for s in anno_header] + row[a1:]))


if __name__ == "__main__":
    main()

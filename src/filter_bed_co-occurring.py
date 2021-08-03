#!/usr/bin/env python3

import sys
import os
import argparse
from random import choice

# Read input options
# Options
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-i", dest="input", type=str,
                    help=("Sorted BED file (can be gzipped). "
                          "Use - for stdin")
                    )
parser.add_argument("-d", dest="dist", type=int, default=500,
                    help=("Minimum distance required between BED entries "
                          "in order to output both")
                    )
parser.add_argument("-c", dest="col", type=int,
                    help=("Which column should be used for picking entries "
                          "(1-based)")
                    )
parser.add_argument("-t", dest="criteria",
                    type=str, help=("Keep [min,max] in a cluster "
                                    "or none [strict]")
                    )

args = parser.parse_args()


# Define auxiliary functions
def open_maybe_gzipped(filename):
    if filename != "-":
        with open(filename, 'rb') as test_read:
            byte1, byte2 = test_read.read(1), test_read.read(1)
            if byte1 and ord(byte1) == 0x1f and byte2 and ord(byte2) == 0x8b:
                f = gzip.open(filename, mode='rt')
            else:
                f = open(filename, 'rt')
        return f
    else:
        return sys.stdin


if args.criteria == "max":
    def keep_previous(current, previous):
        if current == previous:
            return choice([True, False])
        else:
            keep = previous > current
            return keep
elif args.criteria == "min":
    def keep_previous(current, previous):
        if current == previous:
            return choice([True, False])
        else:
            keep = previous < current
            return keep
elif args.criteria == "strict":
    pass
else:
    sys.exit('Comparison should be either "min", "max", or "strict"!')


# Start main script
# Read user input
if args.col is None or args.criteria is None:
    sys.exit("Arguments -c and -t are mandatory!")

infile = args.input
min_dist = args.dist
val_index = args.col - 4

# Start counters
total_entries = 0
entries_outputted = 0

# Parse file
with open_maybe_gzipped(infile) as input_stream:
    # Read first line separately
    # assign values directly to previous_line variable
    previous_line = input_stream.readline().strip()
    total_entries += 1
    previous_bed_info = previous_line.split("\t")
    previous_chrom = previous_bed_info.pop(0)
    previous_start = int(previous_bed_info.pop(0))
    previous_end = int(previous_bed_info.pop(0))
    previous_value = previous_bed_info[val_index]
    try:
        previous_value = float(previous_value)
    except ValueError:
        sys.exit("Line {} in column {} is not numeric. Check user input!"
                 .format(total_entries, args.col))

    # Read rest of the file
    ignore_last = False
    for line in input_stream:
        # import pdb; pdb.set_trace()  # debug
        total_entries += 1
        line = line.strip()
        bed_info = line.split("\t")

        chrom = bed_info.pop(0)
        start = int(bed_info.pop(0))
        end = int(bed_info.pop(0))
        value = bed_info[val_index]
        try:
            value = float(value)
        except ValueError:
            sys.exit("Line {} in column {} is not numeric. "
                     "Check user input!".format(total_entries, args.col))
        # We are in the same chromosome
        if chrom == previous_chrom:
            dist_from_last = abs(start - previous_end)
            # Entries are far enough apart -
            # print previous entry and update values
            if dist_from_last > min_dist:
                if not ignore_last:
                    entries_outputted += 1
                    print(previous_line)

                ignore_last = False
                previous_chrom = chrom
                previous_start = start
                previous_end = start
                previous_value = value
                previous_line = line
            # Entries are too close
            else:
                # We are not being strict -
                # decide which we want to keep
                if args.criteria != "strict":
                    keep_last = keep_previous(value, previous_value)
                    # Current value is desirable - update everything
                    if not keep_last:
                        previous_line = line
                        previous_chrom = chrom
                        previous_start = start
                        previous_end = start
                        previous_value = value
                    # Current value is not desirable - discard it
                    # else: None
                # We are being strict -
                # ignore this and previous entry
                # (but update read buffer)
                else:
                    ignore_last = True
                    previous_line = line
                    previous_chrom = chrom
                    previous_start = start
                    previous_end = start
                    previous_value = value
        # We changed chromosomes -
        # print previous entry and update values
        else:
            if not ignore_last:
                entries_outputted += 1
                print(previous_line)

            ignore_last = False
            previous_chrom = chrom
            previous_start = start
            previous_end = start
            previous_value = value
            previous_line = line
    # Print last line
    if not ignore_last:
        entries_outputted += 1
        print(previous_line)
# We're done

# Print log at the end
entries_removed = (total_entries
                   - entries_outputted)
print("Done!\nProcessed {} BED entries".format(total_entries), file=sys.stderr)
print("Outputted {} entries".format(entries_outputted), file=sys.stderr)
print("Discarded {} entries".format(entries_removed), file=sys.stderr)

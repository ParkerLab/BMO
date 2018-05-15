#!/usr/bin/env python

## measureRawSignal.py ##
## This script will get the number of fragments mapping to the vicinity of a motif.
# Usage: measureRawSignal.py -b input.bam -m motifs.bed(.gz)

from __future__ import print_function
import gzip
import multiprocessing
from optparse import OptionParser
import signal
import sys
import time

import pysam


def worker_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    global samfile
    samfile = pysam.AlignmentFile(input_bam, 'rb')
    if not samfile.check_index():
        sys.exit("Please supply an indexed BAM file.")

# Check if input is gzipped and open files accordingly (John's voodoo here)
def open_maybe_gzipped(filename):
    with open(filename, 'rb') as test_read:
        byte1, byte2 = test_read.read(1), test_read.read(1)
        if byte1 and ord(byte1) == 0x1f and byte2 and ord(byte2) == 0x8b:
            f = gzip.open(filename, mode='rt')
        else:
            f = open(filename, 'rt')
    return f

def measure_motif_signal(motif_line):

    # Process entry and define extended region to look
    bedData = motif_line.strip().split("\t")
    chrom, start, end = bedData[0:3]
    start, end        = int(start), int(end) - 1
    exStart, exEnd    = start - extension, end + extension

    # Assert we are within the chromosomal space
    if exStart < 0:
        exStart = 0
    if exEnd > chromSizes[chrom]:
        exEnd = chromSizes[chrom]

    # Fetch fragments mapping to extended region
    mappingFrags = []
    toExclude    = []
    for read in samfile.fetch(chrom, exStart, exEnd):

        # Correct for the ATAC offset
        #    Note: I'm shifting the entire information in the BAM entry,
        # including the other pair (I use this information on the next step to
        # define the order of the pairs).
        if read.is_reverse:
            integration = read.reference_end - 5
        else:
            integration  = read.reference_start + 4

        fragName = read.query_name

        # Check if integration is within the required coords (pysam can fetch from outside)
        # Action: skip this read altogether
        if integration < exStart or integration > exEnd:
            continue

        # Create mapping and exclusion sets based on motif integration
        if not noIgnOvlp:
            if integration < start or integration > end:
                mappingFrags.append(fragName)
            else:
                toExclude.append(fragName)
        else:
            mappingFrags.append(fragName)

    mappingFrags = set(mappingFrags)
    toExclude    = set(toExclude)
    nValidFrags  = len(mappingFrags.difference(toExclude))

    # append the number of fragments mapping to the motif's extended
    # region, but not within the motif, to the original BED line data
    bedData.append(nValidFrags)

    return bedData


# Start the clock
start_time = time.time()

# Set options
parser = OptionParser()
parser.add_option("-b", dest = "bam", type = str, help = "Sorted and indexed bam file)")
parser.add_option("-m", dest = "motif", type = str, help = "Gzipped motif bed file")
parser.add_option("-e", dest = "ext", type = int, default = 100,
                  help = "Size of the extended region (+- {int} bp). Defaults to 100")
parser.add_option("-p", dest = "parallel", type = int, default = 1,
                  help = "Number of parallel measuring processes to use")
parser.add_option("--no_ignore", dest = "overlap", action = "store_true",
                  help = "Don't ignore integration events at motif")

(options, args) = parser.parse_args()

# Main parameters
input_bam = options.bam
input_bed = options.motif
extension = options.ext
noIgnOvlp = options.overlap


# Get chromosome sizes
samfile = pysam.AlignmentFile(input_bam, 'rb')
chromSizes = {}
for i in range(0, len(samfile.lengths)):
    chromName, chromLen = samfile.references[i], samfile.lengths[i]
    chromSizes[chromName] = chromLen


# Open BED file
try:
    print("Opening motifs file and processing...", file = sys.stderr)
    with open_maybe_gzipped(input_bed) as motifs:
        if options.parallel > 1:
            print("Using {} measuring processes.".format(options.parallel), file=sys.stderr)
            pool = multiprocessing.Pool(processes=options.parallel, initializer=worker_init)
            measured_motifs = pool.imap(measure_motif_signal, motifs, options.parallel)
        else:
            measured_motifs = (measure_motif_signal(motif) for motif in motifs)

        nLines = 0
        for bed_line in measured_motifs:
            print(*bed_line, sep = "\t")
            nLines += 1

    duration = round(time.time() - start_time, 1)
    print("Done!\nProcessed {} entries in {} seconds.".format(nLines, duration), file = sys.stderr)

except KeyboardInterrupt:
    logger.info('Keyboard interrupt received in process {}.'.format(os.getpid()))

    if pool:
        if sys.version_info.major < 3:
            logger.warn('Interrupt handling in multiprocessing programs is not reliable before Python 3, so I may have to exit without cleaning up all worker processes. Sorry.')
            signal.signal(signal.SIGALRM, atactk.util.exit_forcefully)
            signal.alarm(10)

        logger.info('Telling worker processes to terminate...')
        pool.terminate()
        logger.info('Waiting for worker processes to terminate...')
        pool.join()
    logger.info('Exiting.')
    sys.exit(1)

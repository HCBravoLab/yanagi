
import logging
from argparse import ArgumentParser, RawTextHelpFormatter

from lib.ASQuantifier import quantifyEvents

import os, glob

description = \
    "Description:\n\n" + \
    "This subcommand calculates PSI values for a set of splicing events in samples from " + \
    "its underlying segment counts."
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-es", "--events2segs", help="specify file name mapping events to segments. (Prepared through 'segment' subcommand)")
parser.add_argument("-s", "--segs-meta", help="specify file name of the segments meta file. (Prepared through 'segment' subcommand)")
parser.add_argument("-i", "--segCounts-dir", help="specify directory where the samples segment counts files exist. " + \
                    "Each sample should have a .tsv file")
parser.add_argument("-o", "--out-dir", help="specify output directory")
parser.add_argument("-opf", "--out-prefix", help="specify prefix appended to output file names", required=False)


def main():
    args = parser.parse_args()

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    samplesFnames = glob.glob(os.path.join(args.segCounts_dir, "*.tsv"))
    if not args.out_prefix:
        args.out_prefix = os.path.basename(args.events2segs).split(".")[0]

    segIDs = []
    segLens = {}
    with open(args.segs_meta) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            tokens = line.strip().split("\t")
            segIDs.append(tokens[0])
            segLens[tokens[0]] = int(tokens[-1])
    quantifyEvents(args.events2segs, samplesFnames, args.out_dir, args.out_prefix, segIDs, segLens)

    logger.info("Done")
    
if __name__ == '__main__':
    main()

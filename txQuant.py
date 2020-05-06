
import logging
from argparse import ArgumentParser, RawTextHelpFormatter

from lib.txQuantifier import quantifyTxExperiment

import os, glob

description = \
    "Description:\n\n" + \
    "This subcommand runs EM to estimate transcripts abundance in samples given " + \
    "its underlying segment counts."
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-wd", "--work-dir", help="specify directory where the preprocessed annotation files " + \
                    "exist (same output directory used in the preprocess subcommand)")
parser.add_argument("-s", "--segs-meta", help="specify file name of the segments meta file. (Prepared through 'segment' subcommand)")
parser.add_argument("-i", "--segCounts-dir", help="specify directory where the samples segment counts files exist. " + \
                    "Each sample should have a .tsv file")
parser.add_argument("-o", "--out-dir", help="specify output directory")
parser.add_argument("-l", "--l-eff", help="Average Fragment Length")
parser.add_argument("--discover", action='store_true', help="Run correction procedure for unannotated bias correction (default: no)")

#parser.add_argument("-opf", "--out-prefix", help="specify prefix appended to output file names", required=False)


def main():
    args = parser.parse_args()

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    samplesFnames = glob.glob(os.path.join(args.segCounts_dir, "*.tsv"))
    
    quantifyTxExperiment(args.work_dir, samplesFnames, args.out_dir, args.segs_meta, int(args.l_eff), args.discover)

    logger.info("Done")
    
if __name__ == '__main__':
    main()

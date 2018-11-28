
import sys
import logging
from argparse import ArgumentParser, RawTextHelpFormatter

from lib.SegGraphCreator import createSegments



description = \
    "Description:\n\n" + \
    "The genSegs subcommand breaks the preprocessed transcriptome into a set of maximal L-disjoint" + \
    "segments (prepares the segments library)."
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-wd", "--work-dir", help="specify work directory where the preprocessed files " + \
                    "exist (same output directory used in the prep subcommand)", required=True)
parser.add_argument("-l", "--L", help="to create L-disjoint segments (typically equals to read length)",
                    required=True)

parser.add_argument("-o", "--output-prefix", help="specify name prefix of output files", required=False)


def main():
    args = parser.parse_args()

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    if not args.output_prefix:
        args.output_prefix = "segs_"+args.L+".fa"
    createSegments(int(args.L), args.work_dir, args.output_prefix)

    logger.info("Done")
    
if __name__ == '__main__':
    main()

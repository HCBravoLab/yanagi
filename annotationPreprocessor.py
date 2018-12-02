
import sys
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from subprocess import check_call


description = \
    "Description:\n\n" + \
    "This subcommand preprocesses the transcriptome by breaking exons into " + \
    "disjoint exonic bins and find their transcripts mapping"
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-gtf", "--annotation-file", help="specify transcriptome annotation GTF file",
                    required=True)
parser.add_argument("-fa", "--sequence-file", help="specify transcriptome sequences FASTA file",
                    required=True)
parser.add_argument("-o", "--output-dir", help="specify output path", required=True)


def main():
    args = parser.parse_args()

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    # Run Preprocessing Step
    ret = check_call(["Rscript", "R/preprocessing_main.R", args.annotation_file, args.sequence_file, args.output_dir])
    
    if ret != 0:
        logger.info("Preprocessing Failed!")
    else:
        logger.info("Done")
    
if __name__ == '__main__':
    main()

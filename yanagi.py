

import annotationPreprocessor as preprocessAnn
import segmentsGenerator as generateSegments
import segmentsAligner as alignSegments
import psiCalculator as calculatePSI

import logging
import argparse 
import sys

description = "Description:\n\n" + \
              "Yanagi allows you to slice the transcriptome graph into set of maximal L-disjoint segments to be used as reference library for pseudo-alignment of short RNA-seq reads, \n" \
              "calculates segment counts (or segment-pair counts) as local statistics without the need for transcript quantification, \n" \
              "which can be used for gene-level and alternative splicing analysis \n" \
              "For further information, see the help of each subcommand."

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers()

# AnnotationPreprocessor Parser
preprocessAnnSubparser = subparsers.add_parser(
    "preprocess", parents=[preprocessAnn.parser],
    help="Preprocesses transcriptome annotation by breaking exons into disjoint exonic bins and find their transcript mapping.")
preprocessAnnSubparser.set_defaults(which="preprocess")

# SegmentsGenerator Parser
segmentsGeneratorSubparser = subparsers.add_parser(
    "segment", parents=[generateSegments.parser],
    help="Generates a set of maximal L-disjoint segments from the preprocessed transcriptome annotation.")
segmentsGeneratorSubparser.set_defaults(which="segment")

# SegmentsAligner Parser
segmentsAlignerSubparser = subparsers.add_parser(
    "align", parents=[alignSegments.parser],
    help="Pseudo aligns reads (single or paired-end) into the segments and obtain segment counts (single segment or segment pair counts).")
segmentsAlignerSubparser.set_defaults(which="align")

# PSICalculator Parser
psiCalculatorSubparser = subparsers.add_parser(
    "psiCalc", parents=[calculatePSI.parser],
    help="Calculates PSI values for a set of splicing events in samples from its underlying segment counts.")
psiCalculatorSubparser.set_defaults(which="psiCalc")

# Setting logging preferences
logger = logging.getLogger(__name__)

def main():
    try:
        args = parser.parse_args()
        if args.which == "preprocess":
            preprocessAnn.parser = parser  # Setting the module aparser
            preprocessAnn.main()
        elif args.which == "segment":
            generateSegments.parser = parser  # Setting the module aparser
            generateSegments.main()
        elif args.which == "align":
            alignSegments.parser = parser  # Setting the module aparser
            alignSegments.main()
        elif args.which == "psiCalc":
            calculatePSI.parser = parser  # Setting the module aparser
            calculatePSI.main()
            
    except Exception:
        logger.error("Unknown error: {}".format(sys.exc_info()))
        sys.exit(1)
        
if __name__ == '__main__':
    main()

#!/usr/bin/env python

import os
import sys
import logging
from argparse import ArgumentParser

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

def checkDir(Dirname):
    logging.info("Checking folder: '%s'" % Dirname)
    dirname = os.path.abspath(Dirname)
    if not os.path.isdir(dirname):
        logging.error("Oops! Folder: '%s' does not exit. Please check!" % Dirname)
        sys.exit(-1)
    if not os.access(dirname, os.W_OK):
        logging.error("Oops! Folder: '%s' is not writable. Please check!" % Dirname)
        sys.exit(-1)

def checkFile(Filename):
    logging.info("Checking file: '%s'" % Filename)
    filename =  os.path.abspath(Filename)
    if not os.path.isfile(filename):
        logging.error("Oops! File: '%s' does not exit. Please check!" % Filename)
        sys.exit(-1)
    if not os.access(filename, os.R_OK):
        logging.error("Oops! File: '%s' is not readable. Please check!" % Filename)
        sys.exit(-1)


def renameFastaHeader(fasta, sampleName, delimiter, outDir):
    checkFile(fasta)
    checkDir(outDir)
    name = fasta.rsplit('.', 1)[0]
    fasta = os.path.abspath(fasta)
    outDir = os.path.abspath(outDir)
    outputfile = os.path.join(outDir, '%s.headerModified.fa' % name)
    logging.info("Start to rename the header ...")
    with open (outputfile, 'w') as fd:
        with open (fasta, 'r') as fa:
            for line in fa: 
                line = line.strip()
                if line.startswith(">"):
                    line=">%s%s%s" % (sampleName, delimiter, line[1:])
                fd.write('%s\n' % line)
    logging.info("Complete rename the header ...")


if __name__=="__main__":
    parser = ArgumentParser(description='rename the header of a given fasta file')
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument('-f', dest='fasta', help='a fasta format file', type = str)
    parser.add_argument('-n', dest='name', help='name of the sample', type = str)
    parser.add_argument('-d', dest='delim', help="delimiter. Default: '||'", type = str)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()

    if None in [args.fasta, args.name, args.output]:
        print(parser.print_help())
        exit(-1)
    if args.delim == None:
        args.delim = "||"
    renameFastaHeader(args.fasta, args.name, args.delim, args.output)
        


#!/usr/bin/env python

import argparse, sys, os, argparse, shutil, subprocess, gzip, bz2

from pispino.seqtools import *
from pispino.logger import *
from pispino.runcmd import *

__version__    = 1.1

__author__     = "Hyun Soon Gweon"
__copyright__  = "Copyright 2015, Originally from the PIPITS Project"
__credits__    = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__    = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__      = "h.s.gweon@reading.ac.uk"

PEAR                       = "pear"
VSEARCH                    = "vsearch"
FASTQJOIN                  = "fastq-join"
FASTX_FASTQ_QUALITY_FILTER = "fastq_quality_filter"
FASTX_FASTQ_TO_FASTA       = "fastq_to_fasta"


def count_sequences(
    input_dir, 
    filenames_list, 
    logging_file, 
    summary_file, 
    verbose):

    filenameextensions = []
    for filename in filenames_list:
        filenameextensions.append(filename.split(".")[-1].rstrip())
    if len(set(filenameextensions)) > 1:
        logger("Error: More than two types of extensions", logging_file)
        exit(1)
    extensionType = next(iter(filenameextensions))

    # Count
    numberofsequences = 0
    for filename in filenames_list:
        numberofsequences += int(getFileLineCount(input_dir  + "/" + filename, extensionType) / 4)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences in the rawdata!", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of reads: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of reads: " + str(numberofsequences) + "\n")


def reindex_fastq(
    input_dir,
    output_dir,
    sampleids_list,
    filenames_list,
    logging_file,
    summary_file,
    verbose):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:                                                                                                                                                                                                                                   
        shutil.rmtree(output_dir)
        os.mkdir(output_dir) 

    filenameextensions = []
    for filename in filenames_list:
        filenameextensions.append(filename.split(".")[-1].rstrip())
    if len(set(filenameextensions)) > 1:
        logger("Error: More than two types of extensions", logging_file, display = True)
        exit(1)
    extensionType = next(iter(filenameextensions))

    for i in range(len(filenames_list)):
        
        if extensionType == "gz":
            f = gzip.open(input_dir + "/" + filenames_list[i], 'rt')
        elif extensionType == "bz2":
            f = bz2.open(input_dir + "/" + filenames_list[i], 'rt')
        elif extensionType == "fastq":
            f = open(input_dir + "/" + filenames_list[i], 'rt')
        else:
            logger_error("Unknown extension found.", loggint_file)
            exit(1)
        
        outfile = open(output_dir + "/" + sampleids_list[i] + ".fastq", "w")

        logger("Reindexing and saving: " + output_dir + "/" + sampleids_list[i] + ".fastq", logging_file, display = False)

        line_number = 1
        sequence_number = 1
        for line in f:
            if line_number % 4 == 1:
                outfile.write("@" + sampleids_list[i] + "_" + str(sequence_number) + "\n")
                sequence_number += 1
            else:
                outfile.write(line.rstrip() + "\n")
            line_number += 1
        outfile.close()


def skipjoin(
    input_dir_f,
    input_dir_r,
    output_dir,
    sampleids_list,
    base_phred_quality_score,
    joiner_method,
    threads,
    PEAR_parameters,
    logging_file,
    summary_file,
    verbose):

    logger("Skip joining (just using forward reads)", logging_file, display = True)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    for i in range(len(sampleids_list)):

        # If empty, then create empty outputs
        if os.stat(input_dir_f + "/" + sampleids_list[i] + ".fastq").st_size == 0:
            open(output_dir + "/" + sampleids_list[i] + ".fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".discarded.fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".unassembled.forward.fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".unassembled.reverse.fastq", 'a').close()
            continue

        # Simply copy the files                                                                                                                                                                                                                                                                                                                        
        cmd = " ".join(["cp",
                        input_dir_f + "/" + sampleids_list[i] + ".fastq",
                        output_dir + "/" + sampleids_list[i] + ".fastq"])
        run_cmd(cmd, logging_file, verbose)

    # Count
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fastq"
        numberofsequences += int(getFileLineCount(filename) / 4)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences!", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of forward reads: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of forward reads: " + str(numberofsequences) + "\n")


# Join paired-end reads
def join(
    input_dir_f,
    input_dir_r,
    output_dir,
    sampleids_list,
    base_phred_quality_score,
    joiner_method,
    threads,
    PEAR_parameters,
    logging_file,
    summary_file,
    verbose):

    logger("Joining paired-end reads" + " " + "[" + joiner_method + "]", logging_file, display = True)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    for i in range(len(sampleids_list)):

        # If empty, then create empty outputs
        if os.stat(input_dir_f + "/" + sampleids_list[i] + ".fastq").st_size == 0 or os.stat(input_dir_r + "/" + sampleids_list[i] + ".fastq").st_size == 0:
            open(output_dir + "/" + sampleids_list[i] + ".fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".discarded.fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".unassembled.forward.fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".unassembled.reverse.fastq", 'a').close()
            continue
        
        if joiner_method == "PEAR":

            cmd = " ".join([
                PEAR,
                "-f", input_dir_f + "/" + sampleids_list[i] + ".fastq",
                "-r", input_dir_r + "/" + sampleids_list[i] + ".fastq",
                "-o", output_dir  + "/" + sampleids_list[i],
                "-j", threads,
                "-b", base_phred_quality_score,
                "-q 30",
                "-p 0.0001",
                PEAR_parameters])
            run_cmd(cmd, logging_file, verbose)

            # Change name from .assembled.fastq to .joined.fastq
            cmd = " ".join(["mv -f", 
                            output_dir + "/" + sampleids_list[i] + ".assembled.fastq", 
                            output_dir + "/" + sampleids_list[i] + ".fastq"])
            run_cmd(cmd, logging_file, verbose)

        elif joiner_method == "FASTQJOIN":

            cmd = " ".join([FASTQJOIN,
                            input_dir_f + "/" + sampleids_list[i] + ".fastq",
                            input_dir_r + "/" + sampleids_list[i] + ".fastq",
                            "-o",
                            output_dir  + "/" + sampleids_list[i]])
            run_cmd(cmd, logging_file, verbose)

            cmd = " ".join(["mv -f",
                            output_dir  + "/" + sampleids_list[i] + ".joined.fastqjoin",
                            output_dir  + "/" + sampleids_list[i]+ ".fastq"])

            run_cmd(cmd, logging_file, verbose)

        elif joiner_method == "VSEARCH":

            logger("Joining with VSEARCH.", logging_file, display = False)

            cmd = " ".join([
                VSEARCH,
                "--fastq_mergepairs",       input_dir_f + "/" + sampleids_list[i] + ".fastq",
                "--reverse",                input_dir_r + "/" + sampleids_list[i] + ".fastq",
                "--fastqout",               output_dir  + "/" + sampleids_list[i]+ ".fastq",
                "--threads",                threads,
                "--fastq_allowmergestagger",
                "--fastq_maxdiffs 500",
                "--fastq_minovlen 20",
                "--fastq_minmergelen 100"
                ])

            run_cmd(cmd, logging_file, verbose)

    # Count
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fastq"
        numberofsequences += int(getFileLineCount(filename) / 4)
    
    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences after joining. Something is not right.", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of joined reads: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of joined reads: " + str(numberofsequences) + "\n")


def qualityfilter(
    input_dir,
    output_dir,
    sampleids_list,
    base_phred_quality_score,
    FASTX_fastq_quality_filter_q,
    FASTX_fastq_quality_filter_p,
    logging_file,
    summary_file,
    verbose):

    logger("Quality filtering [FASTX]", logging_file, display = True)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    for i in range(len(sampleids_list)):

        if os.stat(input_dir + "/" + sampleids_list[i] + ".fastq").st_size == 0:
            open(output_dir + "/" + sampleids_list[i] + ".fastq", "a").close()
            continue

        cmd = " ".join([
            FASTX_FASTQ_QUALITY_FILTER,
            "-i",  input_dir + "/" + sampleids_list[i] + ".fastq", 
            "-o", output_dir + "/" + sampleids_list[i] + ".fastq", 
            "-q", FASTX_fastq_quality_filter_q,
            "-p", FASTX_fastq_quality_filter_p,
            "-Q" + base_phred_quality_score])
        run_cmd(cmd, logging_file, verbose)

    # Count
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fastq"
        numberofsequences += int(getFileLineCount(filename) / 4)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences after quality filtering...", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of quality filtered reads: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of quality filtered reads: " + str(numberofsequences) + "\n")


def convert(
    input_dir,
    output_dir,
    sampleids_list,
    base_phred_quality_score,
    FASTX_fastq_to_fasta_n,
    logging_file,
    summary_file,
    verbose):

    # Removing reads with \"N\" and FASTA conversion
    logger("Converting FASTQ to FASTA [FASTX] (also removing reads with \"N\" nucleotide if specified with \"--FASTX-n\")", logging_file, display = True)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    fastq_to_fasta_n = ""
    if FASTX_fastq_to_fasta_n:
        pass
    else:
        fastq_to_fasta_n = "-n"

    for i in range(len(sampleids_list)):

        if os.stat(input_dir + "/" + sampleids_list[i] + ".fastq").st_size == 0:
            open(output_dir + "/" + sampleids_list[i] + ".fasta", "a").close()
            continue

        cmd = " ".join([FASTX_FASTQ_TO_FASTA, 
                        "-i",  input_dir + "/" + sampleids_list[i] + ".fastq", 
                        "-o", output_dir + "/" + sampleids_list[i] + ".fasta", 
                        "-Q" + base_phred_quality_score,
                        fastq_to_fasta_n])
        run_cmd(cmd, logging_file, verbose)

    # Count
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fasta"
        numberofsequences += int(getFileLineCount(filename) / 2)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences after converting FASTQ to FASTA...", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of prepped sequences: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of prepped sequences: " + str(numberofsequences) + "\n")


def merge(
    input_dir,
    output_dir,
    sampleids_list,
    logging_file,
    verbose):

    # Merge all into a file
    logger("Merging into a single file", logging_file, display = True)

    outfile = open(output_dir + "/prepped.fasta", "w")
    for i in range(len(sampleids_list)):
        line_index = 1
        logger("Reading " + input_dir + "/" + sampleids_list[i] + ".fasta", logging_file, display = False)

        infile_fasta = open(input_dir + "/" + sampleids_list[i] + ".fasta")
        for line in infile_fasta:
            outfile.write(line.rstrip() + "\n")
    outfile.close()

    logger(BLUE + "... done" + ENDC, logging_file, display = True)


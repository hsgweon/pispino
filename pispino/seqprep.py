#!/usr/bin/env python
import argparse, sys, os, argparse, shutil, subprocess, gzip, bz2

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

HEADER = "\033[95m"
BLUE = "\033[94m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
CYAN = "\033[96m"
RED = "\033[91m"
ENDC = "\033[0m"

# For summary
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def getFileLineCount(filename, extensionType = "uncompressed"):
    if extensionType == "gz":
        f = gzip.open(filename, "r")
        return sum(bl.count(b"\n") for bl in blocks(f))
    elif extensionType == "bz2":
        f = bz2.open(filename, "r")
        return sum(bl.count(b"\n") for bl in blocks(f))
    else:
        f = open(filename, "r")
        return sum(bl.count("\n") for bl in blocks(f))


# Modified to be compatible with PYTHON3 - need to watch out though...
def run_cmd(command, log_file, verbose):
    logger_verbose(command, log_file, verbose)
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

    stdout = p.communicate()[0].decode('UTF-8')
    log_file.write(stdout)

    if verbose:
        sys.stdout.write(stdout)

    p.wait()
    FNULL.close()

    if p.returncode != 0:
        logger_error("None zero returncode: " + command, log_file)
        exit(1)


# VSEARCH outputs copyright info and licence to STDOUT; and the running outputs to STDERRDATA
def run_cmd_VSEARCH(command, log_file, verbose):
    logger_verbose(command, log_file, verbose)
    FNULL = open(os.devnull, 'w')
    if verbose:
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    else:
        p = subprocess.Popen(command, shell=True, stdout=FNULL, stderr=FNULL)

    p.wait()
    FNULL.close()

    if p.returncode != 0:
        logger_error("None zero returncode: " + command, log_file)
        exit(1)


def logger_info(string, log_file):
    import time
    output = time.strftime("%Y-%m-%d %H:%M:%S") + " " + string + "\n"
    sys.stdout.write(output)
    log_file.write(output)


def logger_verbose(string, log_file, verbose):
    import time
    output = time.strftime("%Y-%m-%d %H:%M:%S") + " " + string + "\n"
    if verbose:
        sys.stdout.write(output)
    log_file.write(output)


def logger_error(string, log_file):
    import time
    output = time.strftime("%Y-%m-%d %H:%M:%S") + " ERROR: " + string + "\n"
    sys.stdout.write(output)
    log_file.write(output)


def count_sequences(input_dir, 
                    filenames_list, 
                    logging_file, 
                    summary_file, 
                    verbose):
    

    filenameextensions = []
    for filename in filenames_list:
        filenameextensions.append(filename.split(".")[-1].rstrip())
    if len(set(filenameextensions)) > 1:
        logger_error("More than two types of extensions", logging_file)
        exit(1)
    extensionType = next(iter(filenameextensions))

    numberofsequences = 0
    for filename in filenames_list:
        numberofsequences += int(getFileLineCount(input_dir  + "/" + filename, extensionType) / 4)

    logger_info(GREEN + "... number of reads: " + str(numberofsequences) + ENDC, logging_file)

    # Write to summary_file
    summary_file.write("Number of reads: " + str(numberofsequences) + "\n")


def reindex_fastq(input_dir,
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
        logger_error("More than two types of extensions", logging_file)
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

        logger_verbose("Reindexing and saving: " + output_dir + "/" + sampleids_list[i] + ".fastq", logging_file, verbose)

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


# Just forward reads
def justForwardReads(input_dir_f,
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

    logger_info("!!! N.B. Just using forward reads !!!", logging_file)

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


    # For Summary:
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fastq"
        numberofsequences += int(getFileLineCount(filename) / 4)

    logger_info("... number of forward reads: " + str(numberofsequences), logging_file)
    summary_file.write("Number of forward reads: " + str(numberofsequences) + "\n")


# Join paired-end reads
def join(input_dir_f,
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

    logger_info(HEADER + "Joining paired-end reads" + " " + "[" + joiner_method + "]" + ENDC, logging_file)

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

            cmd = " ".join([PEAR,
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

            logger_verbose("Joining with VSEARCH.", logging_file, verbose)

            cmd = " ".join([VSEARCH,
                            "--fastq_mergepairs",
                            input_dir_f + "/" + sampleids_list[i] + ".fastq",
                            "--reverse",
                            input_dir_r + "/" + sampleids_list[i] + ".fastq",
                            "--fastqout",
                            output_dir  + "/" + sampleids_list[i]+ ".fastq",
                            "--threads", 
                            threads,
                            "--fastq_allowmergestagger",
                            "--fastq_maxdiffs 500",
                            "--fastq_minovlen 20",
                            "--fastq_minmergelen 100"])

            run_cmd_VSEARCH(cmd, logging_file, verbose)

    # For Summary:
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fastq"
        numberofsequences += int(getFileLineCount(filename) / 4)
    
    logger_info(GREEN + "... number of joined reads: " + str(numberofsequences) + ENDC, logging_file)
    summary_file.write("Number of joined reads: " + str(numberofsequences) + "\n")


def qualityfilter(input_dir,
                  output_dir,
                  sampleids_list,
                  base_phred_quality_score,
                  FASTX_fastq_quality_filter_q,
                  FASTX_fastq_quality_filter_p,
                  logging_file,
                  summary_file,
                  verbose):

    logger_info(HEADER + "Quality filtering [FASTX]" + ENDC, logging_file)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    for i in range(len(sampleids_list)):

        if os.stat(input_dir + "/" + sampleids_list[i] + ".fastq").st_size == 0:
            open(output_dir + "/" + sampleids_list[i] + ".fastq", "a").close()
            continue

        cmd = " ".join([FASTX_FASTQ_QUALITY_FILTER,
                        "-i",  input_dir + "/" + sampleids_list[i] + ".fastq", 
                        "-o", output_dir + "/" + sampleids_list[i] + ".fastq", 
                        "-q", FASTX_fastq_quality_filter_q,
                        "-p", FASTX_fastq_quality_filter_p,
                        "-Q" + base_phred_quality_score])
        run_cmd(cmd, logging_file, verbose)


    # For Summary:
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fastq"
        numberofsequences += int(getFileLineCount(filename) / 4)

    logger_info(GREEN + "... number of quality filtered reads: " + str(numberofsequences) + ENDC, logging_file)
    summary_file.write("Number of quality filtered reads: " + str(numberofsequences) + "\n")


def convert(input_dir,
            output_dir,
            sampleids_list,
            base_phred_quality_score,
            FASTX_fastq_to_fasta_n,
            logging_file,
            summary_file,
            verbose):

    # Removing reads with \"N\" and FASTA conversion
    logger_info(HEADER + "Converting FASTQ to FASTA [FASTX] (also removing reads with \"N\" nucleotide if specified with \"--FASTX-n\")" + ENDC, logging_file)

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



    # For Summary:
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + sampleids_list[i]+ ".fasta"
        numberofsequences += int(getFileLineCount(filename) / 2)

    logger_info(GREEN + "... number of prepped sequences: " + str(numberofsequences) + ENDC, logging_file)
    summary_file.write("Number of prepped sequences: " + str(numberofsequences) + "\n")


def merge(input_dir,
          output_dir,
          sampleids_list,
          logging_file,
          verbose):

    # Merge all into a file
    logger_info(HEADER + "Merging into a single file" + ENDC, logging_file)

    outfile = open(output_dir + "/prepped.fasta", "w")
    for i in range(len(sampleids_list)):
        line_index = 1
        logger_verbose("Reading " + input_dir + "/" + sampleids_list[i] + ".fasta", logging_file, verbose)

        infile_fasta = open(input_dir + "/" + sampleids_list[i] + ".fasta")
        for line in infile_fasta:
            outfile.write(line.rstrip() + "\n")
    outfile.close()

    logger_info(GREEN + "... done" + ENDC, logging_file)


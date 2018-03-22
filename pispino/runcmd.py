#!/usr/bin/env python

import os, subprocess, sys
from pispino.logger import *
from pispino.colours import *

__version__    = 2.0

__author__      = "Hyun Soon Gweon"
__copyright__   = "Copyright 2015, The PIPITS Project"
__credits__     = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__     = "GPL"
__maintainer__  = "Hyun Soon Gweon"
__email__       = "h.s.gweon@reading.ac.uk"

def run_cmd(command, log_file, verbose):
    
    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    for l in p.stdout:
        if verbose:
            logger(str(l, 'utf-8').rstrip(), log_file, display = True, timestamp = False)
        else:
            logger(str(l, 'utf-8').rstrip(), log_file, display = False, timestamp = False)

    p.wait()
    FNULL.close()

    if p.returncode != 0:
        logger("Error: None zero returncode: " + command, log_file, display = True)
        exit(1)

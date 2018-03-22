#!/usr/bin/env python

import sys, time
from pispino.colours import *

def logger(string, log_file, display, timestamp = True):

	if timestamp:
		output = RED + time.strftime("%Y-%m-%d %H:%M:%S") + ENDC + " " + string + "\n"
	else:
		output = string + "\n"

	log_file.write(output)

	if display:
		sys.stdout.write(output)

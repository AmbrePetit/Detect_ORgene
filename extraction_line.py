#!/usr/local/bin/python

import sys
import re

gff_file=sys.argv[1]
output=sys.argv[2]

fg= open(gff_file,"r")
sg= open(output, "w")

for line in fg:
	if "gene" in line:
		sg.write(line)

fg.close()
sg.close()

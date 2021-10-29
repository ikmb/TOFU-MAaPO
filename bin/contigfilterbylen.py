#!/usr/bin/env python3

#Removes all contigs with less than, in this case, 1500 bp.

import sys

a_file=open(sys.argv[2])
limit=int(sys.argv[1])
for i, line in enumerate(a_file):
    b=a_file.readline()
    if not line.startswith(">"): print(line.strip())
    else:
        lencheck=line.split()[3]
        lencheck=int(lencheck.replace('len=',''))
        if lencheck>limit:
            print(line, b, end='',sep='')
#!/usr/bin/python
#this script merge idxstat file into one
#usage: python get_count_table.py *.idxstat.txt > count.txt
#result: contig_id, length, counts

import sys
flag = 0
ids = []
length = []
count = []
fname = []
for f in sys.argv[1:]:
    fname.append(f)
    co = []
    for n,line in enumerate(open(f)):
        #print line
        if(line[:1] == '*'):
            continue
        spl = line.strip().split('\t')
        if(flag == 0):
            ids.append(spl[0])
            length.append(spl[1])
        su = int(spl[2])+int(spl[3])
        co.append(str(su))
    count.append(co)
    co = []
    flag = 1
names = '\t'.join(fname)
print ('\t'.join(['contig','length',names]))
for i in range(0,len(ids)):
    tco = []
    for j in range(0,len(count)):
        tco.append(count[j][i])
    result =[ids[i],length[i],'\t'.join(tco)]
    print ('\t'.join(result))

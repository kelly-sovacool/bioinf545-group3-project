"""
countKegg.py input output

**Inputs and outputs**

The Input is the .gene file from the geneList.sh tool.

The Output file is a kegg count file. (.txt)

"""

import sys

KeggDict = {}
keggCount = {}

with open("data/metagenome/bwa_DB/IGC.kegg") as file:
    for line in file:
        (key, val) = line.split()
        KeggDict[key] = val

geneList = []

with open(sys.argv[1]) as genes:
    for i in genes.readlines():
        line = i.split()
        if line[1] != "*":
            geneList.append((eval(line[0]), line[1]))

for read in geneList:
    KEGG = KeggDict[read[1]]
    if KEGG in keggCount:
        keggCount[KEGG] = keggCount[KEGG] + read[0]
    else:
        keggCount[KEGG] = read[0]

# print keggCount["K13730"]

Klist = keggCount.items()


def getKey(item):
    return item[1]


Ksort = sorted(Klist, key=getKey, reverse=True)

outfile = sys.argv[2]

with open(outfile, "w") as out:
    for item in Ksort:
        out.write("".join(str(item)) + "\n")

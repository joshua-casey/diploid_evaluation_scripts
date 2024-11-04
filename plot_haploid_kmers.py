import pysam
import sys

from Bio import SeqIO
from Bio.Seq import Seq                                                 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.rcParams["figure.figsize"] = (20, 20)
import matplotlib.patches as mpatches

import re

f = open(sys.argv[1])
occs = []
for line in f:
    a = line.strip().split("\t")
    occs.append((int(a[3]), int(a[5])))
    
ratios = [0 for _ in range(10)]
xs = range(0, 100, 10)
for a, b in occs:
    if (a + b) > 0:
        c = int(10 * (a / (a+b)))
        if c == 10:
            c = 9
        ratios[c] += 1
plt.bar(xs, ratios, width=10, align='edge')
plt.savefig("kmer_ratio.png", bbox_inches="tight")

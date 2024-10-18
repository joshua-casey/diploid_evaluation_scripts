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

records = SeqIO.parse(sys.argv[1], "fasta")
bin_size = 1000
bin_plot_multiplier = 100
kmer_forward = "TTAGGG"
kmer_reverse = str(Seq(kmer_forward).reverse_complement())

seq = {}
for record in records:
    seq[record.id] = str(record.seq).upper()
    
occ_fwd = {}
occ_rev = {}

for rn in seq:
    st = str(seq[rn]).upper()
    
    print(rn, len(st))
    pos_fwd = [m.start() for m in re.finditer(kmer_forward, st)]
    pos_rev = [m.start() for m in re.finditer(kmer_reverse, st)]
    
    occ_fwd[rn] = [0 for _ in range(len(st) // bin_size + 2)]
    occ_rev[rn] = [0 for _ in range(len(st) // bin_size + 2)]
    
    for i in pos_fwd:
        occ_fwd[rn][i // bin_size] += 1
    for i in pos_rev:
        occ_rev[rn][i // bin_size] += 1
        
values = {}
for rn in occ_fwd:
    x = []
    y = []
    j = 0
    for i in range(0, len(occ_fwd[rn]), bin_plot_multiplier):
        left = sum(occ_fwd[rn][i:i+bin_plot_multiplier])
        right = sum(occ_rev[rn][i:i+bin_plot_multiplier])
        value = left + right
        y.append(value)
        x.append(j)
        j += 1
    values[rn] = y

plt.rcParams["figure.figsize"] = (20, 5)
for rn in values:
    x = range(0, len(occ_fwd[rn]), bin_plot_multiplier)
    y = values[rn]
    plt.figure()
    plt.text(0.5, 0.8, rn, fontsize=30, transform=plt.gcf().transFigure)
    plt.bar(x, y, width=100)
    plt.savefig("telomere/%s.png" % rn, bbox_inches="tight")

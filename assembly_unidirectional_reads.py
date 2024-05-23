import pysam
import sys

from Bio import SeqIO
from Bio.Seq import Seq                                                 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 22})
import matplotlib.patches as mpatches

if len(sys.argv) < 2:
    print("Usage: %s <bam file>" % (sys.argv[0]))
    exit 1

align_file = sys.argv[1]

samfile = pysam.AlignmentFile(align_file, "rb")

direction_rev = {}
direction_fwd = {}

for i in range(len(samfile.references)):
    direction_rev[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 2)]
    direction_fwd[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 2)]

current_rn = "a"
current_rs = 0

samfile = pysam.AlignmentFile(align_file, "rb")

for align in samfile:
    if not align.is_secondary:
        get_read_length = align.infer_read_length()
        
        # get left-most position
        
        rs = align.reference_start
        re = align.reference_end
        rn = align.reference_name
        
        if current_rn != rn or rs - current_rs >= 1000000:
            current_rn = rn
            current_rs = rs
            print("Processing", rn, rs)
        
        if re is not None:
            fr = rs // 1000
            to = re // 1000

            if align.is_reverse:
                for i in range(fr, to+1):
                    direction_rev[rn][i] += 1
            else:
                for i in range(fr, to+1):
                    direction_fwd[rn][i] += 1
                    
values = {}
for rn in direction_rev:
    x = []
    y = []
    j = 0
    for i in range(0, len(direction_rev[rn]), 10):
        left = sum(direction_rev[rn][i:i+10])
        right = sum(direction_fwd[rn][i:i+10])
        if left + right == 0:
            value = random.randint(50, 60)
        else:
            value = 100 * max(left, right) / (left + right)
        y.append(value)
        x.append(j)
        j += 1
    values[rn] = y
    
plt.rcParams["figure.figsize"] = (20, 5)

for rn in values:
    x = range(0, len(direction_rev[rn]), 10)
    y = values[rn]
    plt.text(0.5, 0.8, rn, fontsize=30, transform=plt.gcf().transFigure)
    plt.bar(x, y, width=100)
    print(y[:5], y[-5:])
    print(rn)
    plt.ylim(45, 100)
    plt.show()

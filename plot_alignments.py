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

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
bin_size = 1000
bin_plot_multiplier = 10
softclip_percentage = 0.1
softclip_hardmin = 500


direction_rev = {}
direction_fwd = {}

count_softclips_left = {}
count_softclips_right = {}
count_coverage_left = {}
count_coverage_right = {}
count_softclips = {}

for i in range(len(samfile.references)):
    count_softclips_left[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 1)]
    count_softclips_right[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 1)]
    count_coverage_left[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 1)]
    count_coverage_right[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 1)]
    count_softclips[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 1)]
    direction_rev[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 2)]
    direction_fwd[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // bin_size + 2)]

current_rn = None
current_rs = 0

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

for align in samfile:
    if not align.is_secondary:
        get_read_length = align.infer_read_length()
        
        # get left-most position
        
        rs = align.reference_start
        re = align.reference_end
        rn = align.reference_name
        
        if current_rn is None or current_rn != rn or rs - current_rs >= 1000000:
            current_rn = rn
            current_rs = rs
            print("Processing", rn, rs)
        
        if re is not None:
            fr = rs // bin_size
            to = re // bin_size

            if align.is_reverse:
                for i in range(fr, to+1):
                    direction_rev[rn][i] += 1
            else:
                for i in range(fr, to+1):
                    direction_fwd[rn][i] += 1
        
        
            get_cigar = align.cigartuples
            
            if get_cigar is not None:
                # left clips
                left_clip = 0
                start = 0
                while get_cigar[start][0] == 4 or get_cigar[start][0] == 5:
                    left_clip += get_cigar[start][1]
                    start += 1
                    
                    if start >= len(get_cigar):
                        break
                
                # right clips
                right_clip = 0
                start = len(get_cigar) - 1
                while get_cigar[start][0] == 4 or get_cigar[start][0] == 5:
                    right_clip += get_cigar[start][1]
                    start -= 1
                    
                    if start < 0:
                        break
                
                if left_clip >= softclip_percentage * get_read_length or left_clip >= softclip_hardmin:
                    count_softclips_left[rn][rs // bin_size] += 1
                    count_softclips[rn][rs // bin_size] += 1
                count_coverage_left[rn][rs // bin_size] += 1
                
                if right_clip >= softclip_percentage * get_read_length or right_clip >= softclip_hardmin:
                    count_softclips_right[rn][re // bin_size] += 1
                    count_softclips[rn][rs // bin_size] += 1
                count_coverage_right[rn][re // bin_size] += 1

# Reads Directionality

values = {}
for rn in direction_rev:
    x = []
    y = []
    j = 0
    for i in range(0, len(direction_rev[rn]), bin_plot_multiplier):
        left = sum(direction_rev[rn][i:i+bin_plot_multiplier])
        right = sum(direction_fwd[rn][i:i+bin_plot_multiplier])
        if left + right == 0:
            value = 50
        else:
            value = 100 * max(left, right) / (left + right)
        y.append(value)
        x.append(j)
        j += 1
    values[rn] = y

plt.rcParams["figure.figsize"] = (20, 5)
for rn in values:
    x = range(0, len(direction_rev[rn]), bin_plot_multiplier)
    y = values[rn]
    plt.figure()
    plt.text(0.5, 0.8, rn, fontsize=30, transform=plt.gcf().transFigure)
    plt.bar(x, y, width=bin_plot_multiplier)
    plt.ylim(45, 100)
    plt.savefig("direction/%s.png" % rn, bbox_inches="tight")
    
# Reads Coverage

values = {}
for rn in direction_rev:
    x = []
    y = []
    j = 0
    for i in range(0, len(direction_rev[rn]), bin_plot_multiplier):
        left = sum(direction_rev[rn][i:i+bin_plot_multiplier])
        right = sum(direction_fwd[rn][i:i+bin_plot_multiplier])
        
        value = (left + right) / 10
        
        if value > 200:
            value = 200
        
        y.append(value)
        x.append(j)
        j += 1
    values[rn] = y
    
for rn in values:
    x = range(0, len(direction_rev[rn]), 10)
    y = values[rn]
    plt.figure()
    plt.text(0.5, 0.8, rn, fontsize=30, transform=plt.gcf().transFigure)
    plt.bar(x, y, width=bin_plot_multiplier)
    plt.ylim(45, 100)
    plt.savefig("coverage/%s.png" % rn, bbox_inches="tight")
    
# Reads Clipping

values3 = {}
for rn in direction_rev:
    x = []
    y = []
    
    for i in range(0, len(direction_rev[rn]), bin_plot_multiplier):
        left = sum(count_softclips_left[rn][i:i+bin_plot_multiplier])
        right = sum(count_softclips_right[rn][i:i+bin_plot_multiplier])
        
        value = (left + right) / bin_plot_multiplier
        
        y.append(value)
        x.append(j)
        j += 1
    
    values3[rn] = y
    
for rn in values3:
    x = range(0, len(direction_rev[rn]), bin_plot_multiplier)
    y = values[rn]
    plt.figure()
    plt.text(0.5, 0.8, rn, fontsize=30, transform=plt.gcf().transFigure)
    plt.bar(x, y, width=bin_plot_multiplier)
    plt.ylim(45, 100)
    plt.savefig("clips/%s.png" % rn, bbox_inches="tight")
    
get_covs = []
get_clips = []

for rn in values:
    get_covs.append(values[rn])
    get_clips.append(values3[rn])
    
plt.rcParams.update({'font.size': 22})
plt.rcParams["figure.figsize"] = (50, 20)
pos = range(1, len(get_covs)*2+1, 2)
pos2 = range(2, len(get_covs)*2+1, 2)

c = 'black'
bp1 = plt.boxplot(get_covs, positions=pos, boxprops=dict(color=c),
            capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c))
c = 'red'
bp4 = plt.boxplot(get_clips, positions=pos2, boxprops=dict(color=c),
            capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c))

xs = [1.5]
while len(xs) < 20:
    xs.append(xs[-1] + 2)
xticks = []
y = 1
while y < 20:
    xticks.append("chr%d" % y)
    y += 1
    
xticks.append("chrX")

a = mpatches.Patch(color='black', label='Coverage')
b = mpatches.Patch(color='red', label='Clipped reads')

plt.legend(handles=[a, b], loc='upper right')

plt.ylim(0, 150)
plt.xticks(xs, xticks)
plt.savefig("coverage/boxplot.png", bbox_inches="tight")


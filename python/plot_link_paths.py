from collections import defaultdict
from scipy import interpolate
import numpy as np
import matplotlib.pylab as plt
import parser_wrapper
import glob
import itertools
import matplotlib.patches as mpatches
import argparse

parser = argparse.ArgumentParser(description='Plots link occupancy')
parser.add_argument('--file', type=str, help='Metric file')
parser.add_argument('--sofile', type=str, help='Library file for parser', default='libmetrics_parser.dylib')
parser.add_argument('--metric', type=str, help='Metric id', default='path_bytes')
parser.add_argument('--x_min', type=float, default=0)
parser.add_argument('--x_max', type=float, default=2000)
args = parser.parse_args()

INTERESTING_LINKS = ['N0->N1', 'N4->N5', 'N8->N9', 'N12->N13']

p = parser_wrapper.MetricsParser(args.file, args.sofile)
data = p.Parse(args.metric, '.*', deltas=True)
print data
exit

ax_f, axarr = plt.subplots(len(INTERESTING_LINKS), sharex=True, sharey=True)

def SrcDstLabel(src, dst):
    s = str(src) + u'\u2192' + str(dst)
    return s.replace('N', '')

def AggFromPath(path):
    path = path.split('[')[1].split(']')[0]
    pieces = path.split('->')
    return SrcDstLabel(pieces[0], pieces[-1])

cm = plt.get_cmap('hot')

NUM_COLORS=5
colors = itertools.cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
color_map = {}

def GetColor(label):
    if label in color_map:
        return color_map[label]
    return color_map.setdefault(label, colors.next())

GetColor(SrcDstLabel(0, 1))
GetColor(SrcDstLabel(2, 3))
GetColor(SrcDstLabel(6, 7))
GetColor(SrcDstLabel(10, 11))

for i, link in enumerate(INTERESTING_LINKS):
    ax = axarr[i]

    xs = []
    fs = []
    labels = []

    for key, value in data.items():
        assert(key[0] == args.metric)
        path = key[1]
        if link in path:
            x, y = value
            x = np.array(x, dtype=np.float64) * 0.000000000001
            y = np.array(y, dtype=np.float64) * (100.0 / 1000.0 / 1000.0 / 1000.0) * 8
            x, y = parser_wrapper.Bin(x, y, 10)
            xs.append(x)
            fs.append(interpolate.interp1d(x,y, bounds_error=False, fill_value=0))
            labels.append(AggFromPath(path))

    if len(xs) == 0:
        continue
            
    max_x = max(len(i) for i in xs)
    x = None
    for xi in xs:
        if len(xi) == max_x:
            x = xi
    
    ys = [f(x) for f in fs]
    colors_list = [GetColor(i) for i in labels]
    ngons = ax.stackplot(x, ys, labels=labels, colors=colors_list)
#    ngons[0].set_hatch('//')

    ax.set_ylabel('Gbps')
    ax.legend(loc=1)

#color_items = color_map.items()
#ax.legend([plt.Rectangle((0, 0), 1, 1, fc=v) for _, v in color_items],
#            [k for k, _ in color_items], ncol=2, loc=2)

ax_f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in ax_f.axes[:-1]], visible=False)
plt.xlim([args.x_min, args.x_max])
plt.ylim([0,0.999])
plt.xlabel('seconds')
plt.show()

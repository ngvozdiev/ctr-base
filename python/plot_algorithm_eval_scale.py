from collections import defaultdict
from scipy import interpolate
import numpy as np
import matplotlib.pylab as plt
import parser_wrapper
import glob
import itertools
import matplotlib.patches as mpatches
import argparse

OPTIMIZERS = ['CTR', 'CTRNFC', 'MinMaxLD']

parser = argparse.ArgumentParser(description='Plots link occupancy')
parser.add_argument('--file', type=str, help='Metric file')
parser.add_argument('--sofile', type=str, help='Library file for parser', default='libmetrics_parser.dylib')
args = parser.parse_args()

def PlotCDF(x, label):
    x = np.sort(x)
    y = np.arange(len(x))/float(len(x))
    plt.plot(x, y, label=label)            

def GetFlowDistribution(path_stretches, flow_counts_per_path, p):
    assert(len(path_stretches) == len(flow_counts_per_path))
    total_count = np.sum(flow_counts_per_path)
    limit = total_count * p

    values = sorted(zip(path_stretches, flow_counts_per_path))
    path_stretches, flow_counts_per_path = zip(*values)
    
    total = 0
    for i in range(len(path_stretches)):
        total += flow_counts_per_path[i]
        if total >= limit:
            return path_stretches[i]
    assert(False)
    
p = parser_wrapper.MetricsParser(args.file, args.sofile)
for opt in OPTIMIZERS:
    data = []
    path_stretch_data = p.Parse('opt_path_stretch_ms', '.*{}$'.format(opt))
    flow_count_data = p.Parse('opt_path_flow_count', '.*{}$'.format(opt))

    for k, path_stretch_v in sorted(path_stretch_data.items()):
        flow_count_v = flow_count_data[('opt_path_flow_count', k[1])]
        tm_id = k[1].split(':')[1]

        v = np.sum(path_stretch_v[1] * flow_count_v[1])
        #v = GetFlowDistribution(path_stretch_v[1], flow_count_v[1], 1.0)
        data.append((tm_id, v))

    data = sorted(data)
    print data
    plt.plot([i[1] for i in data], label=opt)
plt.legend()
plt.show()
    



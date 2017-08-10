from collections import defaultdict
from scipy import interpolate
import numpy as np
import matplotlib.pylab as plt
import parser_wrapper
import glob
import itertools
import matplotlib.patches as mpatches
import argparse

PERCENTILES = [100, 99, 95]
OPTIMIZERS = ['MinMax', 'B4', 'B4FC', 'CTR', 'CTRNFC', 'MinMaxLD']

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
for metric, flow_count_scale, log_scale in [('opt_path_stretch_rel', True, True), ('opt_path_k_index', False, True)]:
    data_list = []
    for opt in OPTIMIZERS:
        data = defaultdict(dict)
        path_stretch_data = p.Parse(metric, '.*{}$'.format(opt))
        flow_count_data = p.Parse('opt_path_flow_count', '.*{}$'.format(opt))

        for k, path_stretch_v in sorted(path_stretch_data.items()):
            flow_count_v = flow_count_data[('opt_path_flow_count', k[1])]
            tm_id = k[1].split(':')[1]

            for percentile in PERCENTILES:
                if flow_count_scale:
                    v = GetFlowDistribution(path_stretch_v[1], flow_count_v[1], float(percentile/100.0))
                else:
                    v = np.percentile(path_stretch_v[1], percentile)
                data[percentile][tm_id] = (v, tm_id)
        data_list.append(data)

    by_percentile = defaultdict(list)
    for data in data_list:
        for percentile, data_for_percentile in data.items():
            by_percentile[percentile].append([i[1] for i in sorted(data_for_percentile.items())])

    for percentile, data in by_percentile.items():
        by_percentile[percentile] = sorted(zip(*data))

    for percentile, data in by_percentile.items():
        plt.figure()
        plt.title('P_{}_{}'.format(percentile, metric))
        for i in range(len(OPTIMIZERS)):
            y = [k[i][0] for k in data]
            opt = OPTIMIZERS[i]
            PlotCDF(y, opt)
        plt.legend()
        if log_scale:
            plt.xscale('symlog')

plt.show()
    



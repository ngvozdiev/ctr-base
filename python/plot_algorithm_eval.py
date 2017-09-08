from collections import defaultdict
from scipy import interpolate
import numpy as np
import matplotlib.pylab as plt
import parser_wrapper
import glob
import itertools
import matplotlib.patches as mpatches
import argparse
import matplotlib

matplotlib.rcParams.update({'font.size': 12})

PERCENTILES = [100, 99, 95]
OPTIMIZERS = ['MinMax', 'MinMaxLD', 'CTRNFC', 'B4', 'CTR']
OPT_LABELS = ['MinMax (max delay)', 'MinMax (min delay)', 'CTR (no flow counts)', 'B4', 'CTR' ]
FOCUS_ON = 'Iij\.0003'

parser = argparse.ArgumentParser(description='Plots link occupancy')
parser.add_argument('--file', type=str, help='Metric file')
parser.add_argument('--sofile', type=str, help='Library file for parser', default='libmetrics_parser.dylib')
args = parser.parse_args()

def PlotCDF(x, label):
    x = np.sort(x)
    y = np.arange(len(x))/float(len(x))
    plt.plot(x, y, label=label, linewidth=2)            

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
for metric, flow_count_scale, log_scale, x_label, x_min, x_max, y_min in [('opt_path_stretch_rel', True, True, 'per-flow shortest path delay multiplier', 0.8, 100, 0.3), ('opt_path_k_index', False, True, 'k index per path', -1, 11000, 0.0), ('opt_path_count', False, False, 'per-aggregate path count', 0, 10, 0.0), ('opt_path_stretch_ms', True, False, 'per-flow stretch from shortest path (ms) ', 0.1, 50, 0.5)]:
    data_list = []
    for opt in OPTIMIZERS:
        data = defaultdict(dict)
        print metric, opt, 'A'
        path_stretch_data = p.Parse(metric, '.*{}$'.format(opt))
        flow_count_data = p.Parse('opt_path_flow_count', '.*{}$'.format(opt))
        print metric, opt

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
            opt = OPT_LABELS[i]
            PlotCDF(y, opt)
        plt.legend()
        plt.xlabel(x_label)
        plt.ylabel('CDF')
        plt.xlim([x_min, x_max])
        plt.ylim([y_min, 1.05])
        if log_scale:
            plt.xscale('symlog')

plt.figure()
plt.ylabel('CDF')
plt.xlabel('CTR runtime (ms)')
plt.ylim([-0.05, 1.05])

data = []
runtime_data = p.Parse('ctr_runtime_ms', '.*'.format(opt))
for _, v in runtime_data.items():
    data.append(v[1][0])
PlotCDF(data, 'CTR (cold run)')

data = []
runtime_data = p.Parse('ctr_runtime_cached_ms', '.*'.format(opt))
for _, v in runtime_data.items():
    data.append(v[1][0])
PlotCDF(data, 'CTR')

plt.legend()
plt.xscale('log')

# plt.figure()
# for i in range(len(OPTIMIZERS)):
#     opt = OPTIMIZERS[i]
#     data = defaultdict(dict)
#     path_stretch_data = p.Parse('opt_path_stretch_ms', '.*{}.*{}$'.format(FOCUS_ON, opt))
#     flow_count_data = p.Parse('opt_path_flow_count', '.*{}.*{}$'.format(FOCUS_ON, opt))
#     print path_stretch_data.keys()
#     assert(len(path_stretch_data) == 1)

#     for k, path_stretch_v in sorted(path_stretch_data.items()):
#         flow_count_v = flow_count_data[('opt_path_flow_count', k[1])]
#         tm_id = k[1].split(':')[1]

#         xs = []
#         ys = []
#         for percentile in range(101):
#             v = GetFlowDistribution(path_stretch_v[1], flow_count_v[1], float(percentile/100.0))
#             xs.append(v)
#             ys.append(percentile)
#         plt.plot(xs, ys, label=OPT_LABELS[i], linewidth=2)

# plt.ylabel('CDF')
# plt.xlabel('per-flow path stretch (ms)')
# plt.legend()
plt.show()

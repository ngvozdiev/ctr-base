import os 
import glob
import subprocess
import argparse

parser = argparse.ArgumentParser(description='Generates TMs based on other TMs.')
parser.add_argument('--tm_util_location', type=str, default='./tm_util', help='Location of tm_util')
parser.add_argument('--csf', type=float, default=1.3, help='Commodity scale factor for generated matrices')
parser.add_argument('--graph_location', type=str, required=True, help='Where to find the graphs')
parser.add_argument('--base_tm_location', type=str, required=True, help='Where to find the base TMs')
args = parser.parse_args()

for f in os.listdir(args.graph_location):
    if f.endswith(".graph"):
        graph_base = f.split('.graph')[0]
        demand_files = glob.glob("{}/{}.*.demands".format(args.base_tm_location, graph_base))
        i = 0
        for demand_file in demand_files:
            if 'tm_util' in demand_file:
                continue

            out = '{}.tm_util_{}.demands'.format(graph_base, i)
            i += 1

            cmd = './tm_util --topology {}/{} --traffic_matrix {} --target_commodity_scale_factor {} --output {}'.format(args.graph_location, f, demand_file, args.csf, out)
            print cmd
            subprocess.check_call(cmd, shell=True)



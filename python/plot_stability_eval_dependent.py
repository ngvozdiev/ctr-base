import subprocess
import numpy as np
import timeit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.widgets import CheckButtons

import matplotlib
matplotlib.rcParams.update({'font.size': 14})

fig, axs = plt.subplots(2, sharex=True)

#plt.plot([0], [0])
#plt.figure('2')
plt.subplots_adjust(bottom=0.35)

cd_capacity_init = 1
cd_delay_init = 10
cd_demand_init = 1
cd_flows_init = 1000
ab_capacity_init = 1
ab_delay_init = 10
ab_demand_init = 1
ab_flows_init = 1000
decay_factor_init = 0

axcolor = 'lightgoldenrodyellow'
ax_cd_capacity = plt.axes([0.25, 0.1, 0.65, 0.01], facecolor=axcolor)
ax_cd_delay = plt.axes([0.25, 0.12, 0.65, 0.01], facecolor=axcolor)
ax_cd_demand = plt.axes([0.25, 0.14, 0.65, 0.01], facecolor=axcolor)
ax_cd_flows = plt.axes([0.25, 0.16, 0.65, 0.01], facecolor=axcolor)
ax_ab_capacity = plt.axes([0.25, 0.18, 0.65, 0.01], facecolor=axcolor)
ax_ab_delay = plt.axes([0.25, 0.2, 0.65, 0.01], facecolor=axcolor)
ax_ab_demand = plt.axes([0.25, 0.22, 0.65, 0.01], facecolor=axcolor)
ax_ab_flows = plt.axes([0.25, 0.24, 0.65, 0.01], facecolor=axcolor)
ax_decay_factor = plt.axes([0.25, 0.26, 0.65, 0.01], facecolor=axcolor)

s_cd_capacity = Slider(ax_cd_capacity, 'C->D capacity', 0.1, 2, valinit=cd_capacity_init)
s_cd_delay = Slider(ax_cd_delay, 'C->D delay', 1, 200, valinit=cd_delay_init)
s_cd_demand = Slider(ax_cd_demand, 'C->D demand', 0.1, 2, valinit=cd_demand_init)
s_cd_flows = Slider(ax_cd_flows, 'C->D flows', 100, 2000, valinit=cd_flows_init)
s_ab_capacity = Slider(ax_ab_capacity, 'A->B capacity', 0.1, 2, valinit=ab_capacity_init)
s_ab_delay = Slider(ax_ab_delay, 'A->B delay', 1, 200, valinit=ab_delay_init)
s_ab_demand = Slider(ax_ab_demand, 'A->B demand', 0.1, 2, valinit=ab_demand_init)
s_ab_flows = Slider(ax_ab_flows, 'A->B flows', 100, 2000, valinit=ab_flows_init)
s_decay_factor = Slider(ax_decay_factor, 'Decay factor', 0, 1, valinit=decay_factor_init)

def Update(val):
    cd_capacity = s_cd_capacity.val
    cd_delay = s_cd_delay.val
    cd_demand = s_cd_demand.val
    cd_flows = s_cd_flows.val
    ab_capacity = s_ab_capacity.val
    ab_delay = s_ab_delay.val
    ab_demand = s_ab_demand.val
    ab_flows = s_ab_flows.val
    decay_factor = s_decay_factor.val
    opt = radio.value_selected
    dampen_ratio = False
    if check.lines[0][0].get_visible():
        dampen_ratio = True

    cmd = '../build/stability_eval_dependent --cd_link_gbps {} --cd_link_ms {} --cd_aggregate_gbps {} --cd_aggregate_flows {} --ab_link_gbps {} --ab_link_ms {} --ab_aggregate_gbps {} --ab_aggregate_flows {} --decay_factor {} --steps {} --opt="{}" --dampen_ratio={}'.format(cd_capacity, cd_delay, cd_demand, int(cd_flows), ab_capacity, ab_delay, ab_demand, int(ab_flows), decay_factor, 10, opt, dampen_ratio)
    #print cmd
    start_time = timeit.default_timer()
    output = subprocess.check_output(cmd, shell=True)
    elapsed = timeit.default_timer() - start_time
    #print elapsed
    split = output.strip().split('\n')
    print output

    axs[0].clear()
    axs[0].set_ylabel('fraction')
    axs[1].clear()
    axs[1].set_ylabel('Mbps')

    num_aggregates = int(split[0])
    aggregate_data = []
    for aggregate_str in split[1:1 + num_aggregates]:
        label, vs = aggregate_str.split(':')
        aggregate_data.append((label, vs))

    num_paths = int(split[1 + num_aggregates])
    path_data = []
    for path_str in split[2 + num_aggregates : 2 + num_aggregates + num_paths]:
        label, xy = path_str.split(':')
        path_data.append((label, xy))

    aggregate_data.sort()
    path_data.sort()

    for label, xy in path_data:
        x, y = zip(*eval(xy))
        axs[0].plot(x, y, '--', label=label, linewidth=2)
    axs[0].legend()
    axs[0].set_ylim([0, 1])
    axs[0].set_xlabel('run #')

    for label, vs in aggregate_data:
        y = eval(vs)
        axs[1].plot(y, label=label)
    axs[1].set_ylim([0, axs[1].get_ylim()[1]])
    axs[1].legend()

s_cd_capacity.on_changed(Update)
s_cd_delay.on_changed(Update)
s_cd_demand.on_changed(Update)
s_cd_flows.on_changed(Update)
s_ab_capacity.on_changed(Update)
s_ab_delay.on_changed(Update)
s_ab_demand.on_changed(Update)
s_ab_flows.on_changed(Update)
s_decay_factor.on_changed(Update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def Reset(event):
    s_cd_capacity.reset()
    s_cd_delay.reset()
    s_cd_demand.reset()
    s_cd_flows.reset()
    s_ab_capacity.reset()
    s_ab_delay.reset()
    s_ab_demand.reset()
    s_ab_flows.reset()
    s_decay_factor.reset()
button.on_clicked(Reset)

rax = plt.axes([0.025, 0.1, 0.1, 0.1], facecolor=axcolor)
radio = RadioButtons(rax, ('CTR', 'B4', 'B4(P)', 'MinMax'), active=0)
radio.on_clicked(Update)

rax = plt.axes([0.025, 0.2, 0.1, 0.05])
check = CheckButtons(rax, ['Preserve ratio'], [False])
check.on_clicked(Update)

#Reset(None)
#plt.tight_layout()
#plt.figure('1')
plt.show()

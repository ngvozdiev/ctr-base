import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.65)
plt.plot([0], [0])

cd_capacity_init = 1
cd_delay_init = 10
cd_demand_init = 1
cd_flows_init = 1000
ab_capacity_init = 1
ab_delay_init = 10
ab_demand_init = 1
ab_flows_init = 1000

axcolor = 'lightgoldenrodyellow'
ax_cd_capacity = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_cd_delay = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_cd_demand = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
ax_cd_flows = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
ax_ab_capacity = plt.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
ax_ab_delay = plt.axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)
ax_ab_demand = plt.axes([0.25, 0.4, 0.65, 0.03], facecolor=axcolor)
ax_ab_flows = plt.axes([0.25, 0.45, 0.65, 0.03], facecolor=axcolor)

s_cd_capacity = Slider(ax_cd_capacity, 'C->D capacity', 0.1, 2, valinit=cd_capacity_init)
s_cd_delay = Slider(ax_cd_delay, 'C->D delay', 1, 20, valinit=cd_delay_init)
s_cd_demand = Slider(ax_cd_demand, 'C->D demand', 0.1, 2, valinit=cd_demand_init)
s_cd_flows = Slider(ax_cd_flows, 'C->D flows', 100, 2000, valinit=cd_flows_init)
s_ab_capacity = Slider(ax_ab_capacity, 'A->B capacity', 0.1, 2, valinit=ab_capacity_init)
s_ab_delay = Slider(ax_ab_delay, 'A->B delay', 1, 20, valinit=ab_delay_init)
s_ab_demand = Slider(ax_ab_demand, 'A->B demand', 0.1, 2, valinit=ab_demand_init)
s_ab_flows = Slider(ax_ab_flows, 'A->B flows', 100, 2000, valinit=ab_flows_init)

def Update(val):
    cd_capacity = s_cd_capacity.val
    cd_delay = s_cd_delay.val
    cd_demand = s_cd_demand.val
    cd_flows = s_cd_flows.val
    ab_capacity = s_ab_capacity.val
    ab_delay = s_ab_delay.val
    ab_demand = s_ab_demand.val
    ab_flows = s_ab_flows.val

    cmd = '../build/stability_eval_dependent --cd_link_gbps {} --cd_link_ms {} --cd_aggregate_gbps {} --cd_aggregate_flows {} --ab_link_gbps {} --ab_link_ms {} --ab_aggregate_gbps {} --ab_aggregate_flows {}'.format(cd_capacity, cd_delay, cd_demand, int(cd_flows), ab_capacity, ab_delay, ab_demand, int(ab_flows))
    output = subprocess.check_output(cmd, shell=True)
    per_path = output.strip().split('\n')
    ax.clear()
    data = []
    for path_str in per_path:
        label, xy = path_str.split(':')
        data.append((label, xy))
    data.sort()

    for label, xy in data:
        x, y = zip(*eval(xy))
        ax.plot(x, y, label=label)
    ax.legend()

s_cd_capacity.on_changed(Update)
s_cd_delay.on_changed(Update)
s_cd_demand.on_changed(Update)
s_cd_flows.on_changed(Update)
s_ab_capacity.on_changed(Update)
s_ab_delay.on_changed(Update)
s_ab_demand.on_changed(Update)
s_ab_flows.on_changed(Update)

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
button.on_clicked(Reset)

Reset(None)
plt.show()

import matplotlib.pylab as plt
import numpy as np

DTYPE = [('id', 'S256'), ('demand_delta', 'f8'), ('demand_count', 'u8'), ('flow_delta', 'f8'), ('fd_p50', 'f8'), ('fd_p90', 'f8'), ('fd_p95', 'f8'), ('fd_p100', 'f8'), ('add_p50', 'f8'), ('add_p90', 'f8'), ('add_p95', 'f8'), ('add_p100', 'f8'), ('up_p50', 'f8'), ('up_p90', 'f8'), ('up_p95', 'f8'), ('up_p100', 'f8'), ('rm_p50', 'f8'), ('rm_p90', 'f8'), ('rm_p95', 'f8'), ('rm_p100', 'f8')]

ctr_data = np.loadtxt('se.data.ctr', dtype=DTYPE)
b4_data = np.loadtxt('se.data.b4', dtype=DTYPE)

ctr_data.sort(order='demand_delta')
print ctr_data[-10:]
b4_data.sort(order='demand_delta')

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(ctr_data['demand_delta'], label='CTR')
axarr[0].plot(b4_data['demand_delta'], label='B4')
axarr[0].legend()
axarr[0].set_ylabel('Demand delta')
axarr[1].plot(ctr_data['flow_delta'], label='CTR')
axarr[1].plot(b4_data['flow_delta'], label='B4')
axarr[1].set_ylabel('Flow delta')
axarr[1].legend()

f, axarr = plt.subplots(2, sharex=True, sharey=True)

def PlotPerAggregateDeltas(v, ax):
    ax.plot(v['fd_p50'], label='p50')
    ax.plot(v['fd_p90'], label='p90')
    ax.plot(v['fd_p95'], label='p95')
    ax.plot(v['fd_p100'], label='p100')
    ax.legend()
    ax.set_xlabel('TM rank')
    ax.set_ylabel('Per-aggregate delta')

ctr_data.sort(order='fd_p100')
b4_data.sort(order='fd_p100')
    
PlotPerAggregateDeltas(ctr_data, axarr[0])
PlotPerAggregateDeltas(b4_data, axarr[1])

plt.show()

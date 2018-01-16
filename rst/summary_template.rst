Load {{tm_load}}, locality {{tm_locality}}
-----------------

.. image:: {{path_ratios_location}}

This CDF shows how much large the total delay experienced by flows is compared 
to the lowest possible total propagation delay. The lowest possible total 
propagation delay is achieved when all flows are routed on their respective 
shortest paths, but may result in queueing delay, if there is not enough 
capacity to fit all flows on their shortst paths. The plot is CDF with one point 
per traffic matrix for each of the different schemes 
(:download:`tgz <{{path_ratios_location_tgz}}>`).  

Interesting points (traffic matrices) from the plot above:

{{path_ratios_interesting_table}}

.. image:: {{max_path_stretch_rel_location}}

This plot again has one point per traffic matrix. For each traffic matrix we 
pick the maximum path stretch out of the stretches of all flows in the network. 
We plot a CDF of these values for each of the different schemes 
(:download:`tgz <{{max_path_stretch_rel_location_tgz}}>`). 

Interesting points (traffic matrices) from the plot above:

{{max_path_stretch_rel_interesting_table}}

.. image:: {{path_count_location}}

This plot shows the number of path that each aggregate has. If this number is 
one the aggregate will be routed over only one path, if this number is more than 
one the aggregate will have to be split. There is one point per aggregate from 
all traffic matrices (:download:`tgz <{{path_count_location_tgz}}>`).

.. image:: {{link_utilization_location}}

Link utilization (:download:`tgz <{{link_utilization_location_tgz}}>`).

.. image:: {{link_scales_location}}

Minimal link scale when keeping delay within X\% of optimal 
(:download:`tgz <{{link_scales_location_tgz}}>`).

Interesting points (traffic matrices) from the plot above:

{{link_scales_interesting_table}}

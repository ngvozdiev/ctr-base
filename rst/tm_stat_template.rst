.. _interesting_{{interesting_index}}:

{{topology_name}} (seed {{tm_seed}}, locality {{locality}}, load {{load}})
--------------------------------------------------------------------------------

Here is a summary of the traffic matrix::

{{tm_summary}} 

The following graphs are not specific to any routing or traffic engineering 
scheme---they describe properties of the combination of topology and traffic 
matrix.

.. image:: {{demand_sizes_location}}

Here is a CDF of the sizes of the individual aggregates in Mbps. There 
is a point for each ingress-egress pair with non-zero demand in the traffic 
matrix (:download:`download plot data <{{demand_sizes_location_tgz}}>`).

.. image:: {{cumulative_demands_location}}

A plot of the length of aggregates' shortest paths versus their cumulative 
demand. On the x axis non-zero aggregates in the traffic matrix are sorted by 
the lengths of their respective shortest paths. On the y axis the ordered 
demands' volumes are summed, so that the y value of the rightmost point is the 
total demand in the entire TM and the leftmost point is the demand of the 
aggregate with the shortest path of smallest length. This plot is useful in 
judging how local (geographically) aggregates are
(:download:`download plot data <{{cumulative_demands_location_tgz}}>`).

.. image:: {{cumulative_demands_hop_location}}

Same as the previous plot, but each aggregate's shortest path is represented as 
hop count instead of as delay 
(:download:`download plot data <{{cumulative_demands_hop_location_tgz}}>`). 

.. image:: {{sp_utilization_location}}

The distrubition of link utilizations when all aggregates are routed on the 
shortest pahts. This is a general measure of how loaded the network is. If the 
maximum link utilization (rightmost point) is above 1 single shortest-path 
routing will not fit the traffic matrix and traffic engineering is needed 
(:download:`download plot data <{{sp_utilization_location_tgz}}>`).

.. image:: {{total_delay_at_link_scale_location}}

This plot shows how adding headroom to links affects the delay of the network. 
If a scheme that minimizes the total delay experienced by all flows is used and 
headroom is added to all links (e.g., to cope wih variability), then total delay 
will increase, as less flows will fit on low-delay paths. The plot below 
captures this dependency. On the x axis there is link scale and on the y axis 
is total delay, normalized so that delay experienced at no headroom is 1 
(:download:`download plot data <{{total_delay_at_link_scale_location_tgz}}>`).

Routing Schemes
^^^^^^^^^^^^^^^

The rest of the graphs compare how differnt routing schemes perform on this 
traffic matrix. There is a separate curve for each different scheme.  

.. image:: {{link_delay_vs_utilization_location}}

A scatter plot with one point for each link in the topology. On the x axis the 
links are ranked based on propagation delay. Each point's y axis value is the 
link's utilization under the routing scheme being tested. If all points' y 
values are below 1, the routing scheme was able to fit the traffic matrix and 
there will be no congestion (barring variability) 
(:download:`download plot data <{{link_delay_vs_utilization_location_tgz}}>`).

.. image:: {{absolute_stretch_location}}

A CDF with each flow's absolute path stretch---the difference between the flow's 
path and its shortest path 
(:download:`download plot data <{{absolute_stretch_location_tgz}}>`).

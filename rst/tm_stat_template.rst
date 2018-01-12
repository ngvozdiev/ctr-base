.. _interesting_{{interesting_index}}:

TM on {{topology_name}} (seed {{tm_seed}}, locality {{locality}}, load {{load}})
------------------------

Here is a summary of the traffic matrix:

{{tm_summary}} 

This is a CDF of the sizes of the individual aggregates in Mbps. There is a 
point for each ingress-egress pair with non-zero demand in the traffic matrix 
(:download:`tgz <{{demand_sizes_location_tgz}}>`):

.. image:: {{demand_sizes_location}}

Interesting points (aggregates) from the CDF above:

{{demand_sizes_interesting_table}}

A plot of the length of aggregates' shortest paths versus their cumulative 
demand. On the x axis non-zero aggregates in the traffic matrix are sorted by 
the lengths of their respective shortest paths. On the y axis the ordered 
demands' volumes are summed, so that the y value of the rightmost point is the 
total demand in the entire TM and the leftmost point is the demand of the 
aggregate with the shortest path of smallest length. This plot is useful in 
judging how local (geographically) aggregates are. 
(:download:`tgz <{{cumulative_demands_location_tgz}}>`):

.. image:: {{cumulative_demands_location}}

The distrubition of link utilizations when all aggregates are routed on the 
shortest pahts. This is a general measure of how loaded the network is. If the 
maximum link utilization (rightmost point) is above 1 single shortest-path 
routing will not fit the traffic matrix and traffic engineering is needed 
(:download:`tgz <{{sp_utilization_location_tgz}}>`):

.. image:: {{sp_utilization_location}}

Interesting points (links) from the CDF above:

{{sp_utilization_interesting_table}}

This plot shows how adding headroom to links affects the delay of the network. 
If a scheme that minimizes the total delay experienced by all flows is used and 
headroom is added to all links (e.g., to cope wih variability), then total delay 
will increase, as less flows will fit on low-delay paths. The plot below 
captures this dependency. On the x axis there is link scale and on the y axis 
is total delay, normalized so that delay experienced at no headroom is 1 
(:download:`tgz <{{total_delay_at_link_scale_location_tgz}}>`):

.. image:: {{total_delay_at_link_scale_location}}

.. _interesting_{{interesting_index}}:

Seed {{tm_seed}} on topology {{topology_name}}
------------------------

.. csv-table:: Summary
   :header: "load", "locality", "demand count", "fraction of all pairs"
   :widths: 10, 10, 15, 15

{{#summary_table_rows}}
   "{{load_value}}", "{{locality_value}}", "{{demand_count_value}}", "{{demand_fraction_value}}"{{/summary_table_rows}}

The distribution of the sizes of the individual demands (:download:`tgz <{{demand_sizes_location_tgz}}>`):

.. image:: {{demand_sizes_location}}

A plot of each aggregate's cumulative demand size versus the length of 
the aggregate's shortest path (:download:`tgz <{{cumulative_demands_location_tgz}}>`):

.. image:: {{cumulative_demands_location}}

The distrubition of link utilizations when all aggregates are routed on the 
shortest pahts (:download:`tgz <{{sp_utilization_location_tgz}}>`):

.. image:: {{sp_utilization_location}}

Total delay vs link scale (:download:`tgz <{{total_delay_at_link_scale_location_tgz}}>`):

.. image:: {{total_delay_at_link_scale_location}}

Summary for load {{tm_load}}, locality {{tm_locality}}
-----------------

Path ratios (:download:`tgz <{{path_ratios_location_tgz}}>`):

.. image:: {{path_ratios_location}}

.. csv-table:: Maximum ratio
   :header: "optimizer", "value", "link"
   :widths: 10, 10, 15
   
{{#max_ratio_table_rows}}
   "{{optimizer}}", "{{value}}", "{{link}}"{{/max_ratio_table_rows}}
   
.. csv-table:: Median ratio
   :header: "optimizer", "value", "link"
   :widths: 10, 10, 15
   
{{#med_ratio_table_rows}}
   "{{optimizer}}", "{{value}}", "{{link}}"{{/med_ratio_table_rows}}

Relative path stretch (:download:`tgz <{{path_stretch_rel_location_tgz}}>`):

.. image:: {{path_stretch_rel_location}}

Max relative path stretch (:download:`tgz <{{max_path_stretch_rel_location_tgz}}>`):

.. image:: {{max_path_stretch_rel_location}}

Path count (:download:`tgz <{{path_count_location_tgz}}>`):

.. image:: {{path_count_location}}

Link utilization (:download:`tgz <{{link_utilization_location_tgz}}>`):

.. image:: {{link_utilization_location}}

Minimal link scale when keeping delay within X\% of optimal (:download:`tgz <{{link_scales_location_tgz}}>`):

.. image:: {{link_scales_location}}


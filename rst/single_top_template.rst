{{topology_name}}
-----------------

Topology summary::

{{topology_summary}}

Distribution of propagation delays of all shortest paths (:download:`tgz <{{node_distances_location_tgz}}>`):

.. image:: {{node_distances_location}}

.. toctree::
   :maxdepth: 1
{{#traffic_matrices}}
   tm_seed_{{tm_seed}}/index {{/traffic_matrices}}

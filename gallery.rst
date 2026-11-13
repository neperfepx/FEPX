.. _gallery:

Gallery
=======

Images of FEPX simulations, visualized with `Neper <https://neper.info>`_.  Click on an image to enlarge it, and follow the link below it to see how it is made.

.. _gallery_crystal_behavior:

Crystal Behavior
----------------

.. grid:: 1 2 3 3
   :gutter: 3
   :class-row: fepx-gallery-row

   .. grid-item::

      .. image:: imgs/gallery/crystal-strength.png
         :class: fepx-gallery-img
         :target: imgs/gallery/crystal-strength.png

      .. container:: fepx-gallery-caption

         Slip-system strength after a multiaxial load cycle, from :doc:`doc/tutorials/multiaxial_clr`

   .. grid-item::

      .. image:: imgs/gallery/crystal-slip.png
         :class: fepx-gallery-img
         :target: imgs/gallery/crystal-slip.png

      .. container:: fepx-gallery-caption

         Accumulated slip on one slip system after two cycles, from :doc:`doc/tutorials/restart`

   .. grid-item::

      .. image:: imgs/gallery/crystal-plastic-strain.png
         :class: fepx-gallery-img
         :target: imgs/gallery/crystal-plastic-strain.png

      .. container:: fepx-gallery-caption

         Equivalent plastic strain under multiaxial loading, from :doc:`doc/tutorials/multiaxial_csr`

   .. grid-item::

      .. image:: imgs/gallery/crystal-hcp.png
         :class: fepx-gallery-img
         :target: imgs/gallery/crystal-hcp.png

      .. container:: fepx-gallery-caption

         Stress-strain curves and slip rates of an hcp single crystal in three orientations, from :doc:`doc/tutorials/hcp_single_crystal`

.. _gallery_microstructures:

Microstructures
---------------

.. grid:: 1 2 3 3
   :gutter: 3
   :class-row: fepx-gallery-row

   .. grid-item::

      .. image:: imgs/gallery/micro-cruciform.png
         :class: fepx-gallery-img
         :target: imgs/gallery/micro-cruciform.png

      .. container:: fepx-gallery-caption

         Polycrystalline cruciform specimen, colored by grain, from :doc:`doc/tutorials/cruciform`

   .. grid-item::

      .. image:: imgs/gallery/micro-periodic.png
         :class: fepx-gallery-img
         :target: imgs/gallery/micro-periodic.png

      .. container:: fepx-gallery-caption

         Periodic polycrystal, colored by grain, from :doc:`doc/tutorials/periodic_rve`

   .. grid-item::

      .. image:: imgs/gallery/micro-phases.png
         :class: fepx-gallery-img
         :target: imgs/gallery/micro-phases.png

      .. container:: fepx-gallery-caption

         Dual-phase polycrystal: austenite (fcc) and ferrite (bcc), from :doc:`doc/tutorials/multiaxial_clr`

   .. grid-item::

      .. image:: imgs/gallery/micro-orientations.png
         :class: fepx-gallery-img
         :target: imgs/gallery/micro-orientations.png

      .. container:: fepx-gallery-caption

         Orientations varying within the grains, set per element, from :doc:`doc/tutorials/overwriting`

.. _gallery_loading:

Loading
-------

.. grid:: 1 2 3 3
   :gutter: 3
   :class-row: fepx-gallery-row

   .. grid-item::

      .. image:: imgs/gallery/load-cruciform.png
         :class: fepx-gallery-img
         :target: imgs/gallery/load-cruciform.png

      .. container:: fepx-gallery-caption

         Stress field of a cruciform specimen under biaxial loading, from :doc:`doc/tutorials/cruciform`

   .. grid-item::

      .. image:: imgs/gallery/load-periodic.png
         :class: fepx-gallery-img
         :target: imgs/gallery/load-periodic.png

      .. container:: fepx-gallery-caption

         Strain field under periodic boundary conditions, from :doc:`doc/tutorials/periodic_rve`

   .. grid-item::

      .. image:: imgs/gallery/load-unloading.png
         :class: fepx-gallery-img
         :target: imgs/gallery/load-unloading.png

      .. container:: fepx-gallery-caption

         Residual stress field after unloading to zero stress, from :doc:`doc/tutorials/mixed_targets`

   .. grid-item::

      .. image:: imgs/gallery/load-cyclic.png
         :class: fepx-gallery-img
         :target: imgs/gallery/load-cyclic.png

      .. container:: fepx-gallery-caption

         Cyclic multiaxial loading, run in two parts with a restart, from :doc:`doc/tutorials/restart`

.. _gallery_output:

Output
------

.. grid:: 1 2 3 3
   :gutter: 3
   :class-row: fepx-gallery-row

   .. grid-item::

      .. image:: imgs/gallery/output-stress.png
         :class: fepx-gallery-img
         :target: imgs/gallery/output-stress.png

      .. container:: fepx-gallery-caption

         Stress field under grip boundary conditions, from :doc:`doc/tutorials/simple_rve`

   .. grid-item::

      .. image:: imgs/gallery/output-work.png
         :class: fepx-gallery-img
         :target: imgs/gallery/output-work.png

      .. container:: fepx-gallery-caption

         Plastic work field, from :doc:`doc/tutorials/multiaxial_csr`

   .. grid-item::

      .. image:: imgs/gallery/output-defrate.png
         :class: fepx-gallery-img
         :target: imgs/gallery/output-defrate.png

      .. container:: fepx-gallery-caption

         Equivalent deformation rate field, from :doc:`doc/tutorials/overwriting`

   .. grid-item::

      .. image:: imgs/gallery/output-python.png
         :class: fepx-gallery-img
         :target: imgs/gallery/output-python.png

      .. container:: fepx-gallery-caption

         Mean stress of each grain, computed in Python from the results, from :doc:`doc/tutorials/postprocessing`

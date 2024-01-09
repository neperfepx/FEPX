.. _file_and_directory_formats_link:

File and Directory Formats
==========================


.. _simulation_directory:

Simulation Directory (:file:`.sim`)
-----------------------------------

Here are details on the :file:`.sim` simulation directory (the :file:`.sim` extension is entirely optional).  The directory contains *inputs* and *results* on the mesh *nodes* and *elements*, over a certain number of *simulation steps*. It is structured as follows:

.. code-block:: console

  simulation.sim
  |-- inputs
  |   |-- job.sh
  |   |-- simulation.cfg
  |   |-- simulation.msh
  |   `-- simulation.tess
  |-- results
  |   |-- elts
  |   |   |-- ori
  |   |   |   |-- ori.step0
  |   |   |   |-- ori.step1
  |   |   |   `-- ...
  |   |   `-- ...
  |   `-- nodes
  |       |-- coo
  |       |   |-- coo.step0
  |       |   |-- coo.step1
  |       |   `-- ...
  |       `-- ...
  `-- [restart]

where

- :file:`inputs` is an input file directory containing the Neper tessellation file (:file:`simulation.tess`, if found in the input directory), the mesh file (:file:`simulation.msh`), the FEPX configuration file (:file:`simulation.cfg`), and all script files (:file:`*.sh`, likely including a job submission file).

- :file:`results` is the result directory.

- :file:`results/nodes` is the node result directory.

- :file:`results/elts` is the element result directory.

- :file:`results/*/<res>` is a directory for result :data:`<res>` of an entity (nodes or elements). The directory contains one file per simulation step, named :file:`<res>.step<nb>`, where :data:`nb` is the step number, ranging from 0 (for the initial state) to the total number of steps.

- :file:`restart` is the restart directory.  It is present only if :data:`restart on` was used.

Results can have integer values, real values, vectorial values or tensorial values. In the result files, values for the different entities (nodes or elements) are written on successive lines, with components written on successive columns (space delimited).  The components are written as follows:

  - vector, :math:`v`: :math:`v_1`, :math:`v_2`, :math:`v_3`;
  - symmetrical tensor, :math:`t`: :math:`t_{11}`, :math:`t_{22}`, :math:`t_{33}`, :math:`t_{23}`, :math:`t_{31}`, :math:`t_{12}` (Voigt notation);
  - skew-symmetrical tensor, :math:`t`: :math:`t_{12}`, :math:`t_{13}`, :math:`t_{23}`;
  - non-symmetrical tensor, :math:`t`:  :math:`t_{11}`, :math:`t_{12}`, :math:`t_{13}`, :math:`t_{21}`, :math:`t_{22}`, :math:`t_{23}`, :math:`t_{31}`, :math:`t_{32}`, :math:`t_{33}`.

The directory also contains a hidden file, :file:`.sim`, containing information on the simulation and the content of the simulation directory.  This file is only for internal use and is formatted as follows:

.. code-block:: console

  ***sim
   **format
     <format>
   **input
    *tess
     <tess_file>
    *tesr
     <tesr_file>
    *msh
     <msh_file>
    *ori
     <ori_file>
    *bcs
     <bcs_file>
    *phase
     <phase_file>
    *config
     <config_file>
   **general
     <cell_nb> <node_nb> <elt_nb> <elset_nb> <part_nb>
    *orides
     <orientation_descriptor>
   **entity <entity>                                     \
    *member                                              |
     <member_nb>                                         |
     <member1> <member2> ...                             | section repeated for each entity
    *result                                              |
     <result_nb>                                         |
     <result1> <result2> ...                             /
   **orispace
    *rodrigues <space_file>
   **step
     <step_nb>
  ***end

.. _mesh_file:

Mesh File (:file:`simulation.msh`)
----------------------------------

Here are details on the native :file:`.msh` (adapted from Gmsh's msh format version :data:`2.2`).  Developers should note that read and write functions are available as `neut_msh_fscanf` and `neut_msh_fprintf`, defined in directories :file:`neut/neut_msh/neut_msh_fscanf` and :file:`neut/neut_msh/neut_msh_fprintf`.

.. code-block:: console

  $MeshFormat
  2.2 <file_type> <data_size>
  $EndMeshFormat
  $MeshVersion
  <mesh_version>
  $EndMeshVersion
  $Domain
  <domain>
  $EndDomain
  $Topology
  <reconstruct_topology>
  $EndTopology
  $Nodes
  <number_of_nodes>
  <node_id> <node_x> <node_y> <node_z>
  ...
  $EndNodes
  $Elements
  <number_of_elements>
  <elt_id> <elt_type> <number_of_tags> <tag1> ... <elt_id_node1> ...
  ...
  $EndElements
  $Periodicity
  <number_of_periodicities>
  <secondary_node_id> <primary_node_id> <per_vect_x> <per_vect_y> <per_vect_z>
  ...
  $EndPeriodicity
  $NSets
  <number_of_nsets>
  <nset1_label>
  <nset_node_nb>
  <nset_node1>
  <nset_node2>
  ...
  <nset2_label>
  ...
  $EndNSets
  $Fasets
  <number_of_fasets>
  <faset1_label>
  <faset_elt_nb>
  <faset_elt_id> <faset_elt_id_node1> ...
  ...
  <faset2_label>
  ...
  $EndFasets
  $NodePartitions
  <number_of_nodes>
  <node_id> <node_partition>
  ...
  $EndNodePartitions
  $PhysicalNames
  <number_of_physical_names>
  <physical_dimension> <physical_id> <physical_name>
  ...
  $EndPhysicalNames
  $ElsetOrientations
  <number_of_elsets> <orientation_descriptor>
  <elset_id> <ori_des1> ...
  ...
  $EndOrientations
  $ElsetCrySym
  <crysym>
  $EndElsetCrySym
  $ElementOrientations
  <number_of_elements> <orientation_descriptor>
  <element_id> ori_des1> ...
  ...
  $EndElementOrientations
  $Groups
  <group_entity>
  <number_of_group_entities>
  <entity_id group>
  ...
  $EndGroups

where

- :data:`$MeshFormat` denotes the beginning of a mesh format field.

- :data:`<file_type>` is equal to :data:`0` for an ASCII file and :data:`1` for a binary file.

- :data:`<data_size>` is an integer equal to the size of the floating point numbers used in the file (= :data:`sizeof (double)`).

- :data:`$EndMeshFormat` denotes the end of a mesh format field.

- :data:`$MeshVersion` denotes the beginning of a mesh version field.

- :data:`<mesh_version>` is the mesh file version (currently :data:`2.2.3`).

- :data:`$EndMeshVersion` denotes the end of a mesh version field.

- :data:`$Domain` denotes the beginning of an optional domain field.

- :data:`<domain>` is the domain.

- :data:`$EndDomain` denotes the end of an optional domain field.

- :data:`$Topology` denotes the beginning of an optional topology field.

- :data:`<reconstruct_topology>` is a boolean indicating whether the topology is to be reconstructed upon parsing or not (use :data:`0` to solve parsing issues).

- :data:`$EndTopology` denotes the end of an optional topology field.

- :data:`$Nodes` denotes the beginning of a node field.

- :data:`<number_of_nodes>` is the number of nodes.

- :data:`<node_id>` is the identifier of a node and ranges from :data:`1` to :data:`<number_of_nodes>`.

- :data:`<node_x>`, :data:`<node_y>` and :data:`<node_z>` are the three coordinates of a node (real numbers).

- :data:`$EndNodes` denotes the end of a node field.

- :data:`$Elements` denotes the beginning of an element field.

- :data:`<number_of_elements>` is the number of elements.

- :data:`<elt_type>` is an integer specifying the type of elements: :data:`15` for a 0D element, :data:`1` for a 1st-order 1D element (2 nodes), :data:`8` for a 2nd-order 1D element (3 nodes), :data:`2` for a 1st-order triangular element (3 nodes), :data:`3` for a 1st-order quadrangular element (4 nodes), :data:`9` for a 2nd-order triangular element (6 nodes), :data:`16` for a 2nd-order quadrangular element (8 nodes), :data:`10` for a 2nd-order quadrangular element (9 nodes), :data:`4` for a 1st-order tetrahedral element (4 nodes), :data:`5` for a 1st-order hexahedral element (8 nodes), :data:`11` for a 2nd-order tetrahedral element (10 nodes), :data:`17` for a 2nd-order hexahedral element (20 nodes), :data:`6` for a 1st-order prismatic element (6 nodes), :data:`18` for a 2nd-order prismatic element (15 nodes).

- :data:`<number_of_tags>` is the number of tags, and :data:`<tag#>` are the tags.  In the general case, the number of tags is equal to 3, the first and second tags are the elset and the third tag is the element partition.  The mesh partition is non-zero only for the higher-dimension elements of a mesh which was previously partitioned.

- :data:`<elt_id_node#>` are the nodes associated to an element.  The number of nodes depends on the element type (`<elt_type>`).

- :data:`$EndElements` denotes the end of an element field.

- :data:`$Periodicity` denotes the beginning of a periodicity field.

- :data:`<number_of_periodicities>` is the number of periodicities.

- :data:`<primary_node_id>` is the identifier of the primary node.

- :data:`<secondary_node_id>` is the identifier of the secondary node.

- :data:`<per_vect_x>` :data:`<per_vect_y>` :data:`<per_vect_z>` are the scaled components of the vector going from the primary node to the secondary node (-1, 0 or 1).

- :data:`$EndPeriodicity` denotes the end of a periodicity field.

- :data:`$NSets` denotes the beginning of an nset field.

- :data:`<number_of_nsets>` is the number of nsets.

- :data:`<nset#_label>` are the labels of the nsets.

- :data:`<nset_node_nb>` is the number of nodes of an nset.

- :data:`<nset_node_id#>` are the identifiers of the nodes of an nset.

- :data:`$EndNSets` denotes the end of an nset field.

- :data:`$Fasets` denotes the beginning of a faset field.

- :data:`<number_of_fasets>` is the number of fasets.

- :data:`<faset#_label>` are the labels of the fasets.

- :data:`<faset_elt_nb>` is the number of elements of a faset.

- :data:`<faset_elt_id>` are the identifiers of the elements of a faset (3D elements adjacent to the boundary).

- :data:`<faset_elt_id_node#>` are the nodes of an element of a faset.

- :data:`$EndFasets` denotes the end of a faset field.

- :data:`$NodePartitions` denotes the beginning of a node partition field.

- :data:`<nodeid_partition>` is the partition of node :data:`<id>` (ranging from 1 to the total number of partitions).

- :data:`$EndNodePartitions` denotes the end of a node partition field.

- :data:`$PhysicalNames` denotes the beginning of a physical name field.

- :data:`<number_of_physical_names>` is the number of physical names.  There are as many names as physical entities, and the physical entities correspond to all tessellation vertices, edges, faces and polyhedra (i.e., mesh 0D, 1D, 2D and 3D elsets).

- :data:`<physical_dimension>` is the dimension of a physical entity and can be equal to 0, 1, 2 or 3.

- :data:`<physical_id>` is the id of a physical entity.  It ranges from 1 to the number of 0D elsets (tessellation vertices) for the 0D entities, 1 to the number of 1D elsets (tessellation edges) for the 1D entities, 1 to the number of 2D elsets (tessellation faces) for the 2D entities and 1 to the number of 3D elsets (tessellation polyhedra) for the 3D entities.

- :data:`<physical_name>` is the name of a physical entity, under the form :data:`<verid>` for 0D elsets (tessellation vertices), :data:`<edgeid>` for 1D elsets (tessellation edges), :data:`<faceid>` for 2D elsets (tessellation faces) and :data:`<polyid>` for 3D elsets (tessellation polyhedra), where :data:`<id>` ranges from 1 to the number of elsets.

- :data:`$EndPhysicalNames` denotes the end of a physical name field.

- :data:`$ElsetOrientations` denotes the beginning of an elset orientation field.

- :data:`$EndElsetOrientations` denotes the end of an elset orientation field.

- :data:`<number_of_elsets>` is the number of elsets.

- :data:`<orientation_descriptor>` is the orientation descriptor.

- :data:`<elset_id>` is the elset id.

- :data:`<ori_des1>`, ... is the orientation, following :data:`<orientation_descriptor>`.

- :data:`$EndElsetOrientations` denotes the end of an elset orientation field.

- :data:`$ElsetCrySym` denotes the beginning of an elset crystal symmetry field.

- :data:`<crysym>` is the crystal symmetry (:data:`triclinic`, :data:`cubic` or :data:`hexagonal`).

- :data:`$EndElsetCrySym` denotes the end of an elset crystal symmetry field.

- :data:`$ElementOrientations` denotes the beginning of an element orientation field.

- :data:`<number_of_elements>` is the number of elements.

- :data:`<element_id>` is the element id.

- :data:`$EndElementOrientations` denotes the end of an element orientation field.

- :data:`$Groups` denotes the beginning of a group field.

- :data:`<group_entity>` is the entity for which groups are defined, which must be :data:`elset`.

- :data:`<number_of_group_entities>` is the number of group entities (number of elsets).

- :data:`<entity_id>` is the id of an entity.

- :data:`<group>` is the group of the entity.

- :data:`$EndGroups` denotes the end of a group field.

.. _ori_file:

Orientation File (:file:`simulation.ori`)
-----------------------------------------

The embedded orientation assignments within the :file:`simulation.msh` may be overridden via an external :file:`simulation.ori` file. This file contains formatting identical to the associated fields in the mesh file and is defined as:

For per-grain (or :data:`Elset`) orientations::

    $ElsetOrientations
    <number_of_ori_entities> <orientation_descriptor>:<orientation_convention>
    <entity_id> <ori_des1> ...
    ...
    $EndElsetOrientations

For per-element orientations::

    $ElementOrientations
    <number_of_ori_entities> <orientation_descriptor>:<orientation_convention>
    <entity_id> <ori_des1> ...
    ...
    $EndElementOrientations

where :data:`<number_of_ori_entities>` is the number of unique orientations defined in the section, :data:`<orientation_descriptor>` is the parameterization for the orientations (see options below), :data:`<orientation_convention>` describes the basis transformation route for the orientations provided, :data:`<entity_id>` is a unique 1-indexed identification number, and :data:`<ori_des*>` are the components of the unique orientation. Available options for :data:`<orientation_convention>` are: :data:`active` or :data:`passive`. Following the usual terminology, an active orientation assumes that which describes a basis transformation from the sample basis to the crystal basis (sample-to-crystal), while a passive orientation convention assumes that which describes a basis transformation from the crystal basis to the sample basis (“crystal-to-sample”).

The following :data:`<orientation_descriptor>` types are available (associated per-line formats are also described):

- For :data:`rodrigues`, each orientation is described by :math:`r_1, r_2, r_3`, where :math:`{\bf r} = {\bf t} \tan{ (\omega / 2)}` (see: :data:`axis-angle` for definitions of :math:`{\bf t}` and :math:`\omega`). The per-line format is::

    <entity_id> <r_1> <r_2> <r_3>

- For :data:`euler-bunge`, each orientation is described by :math:`\phi_1, \theta, \phi_2`, where :math:`\phi_1` is the rotation about the :math:`z` axis, :math:`\theta` is the rotation about the :math:`z^{\prime}` axis, and :math:`\phi_2` is the rotation about the :math:`z^{\prime \prime}` axis, all in degrees). The per-line format is::

    <entity_id> <phi_1> <Phi> <phi_2>

- For :data:`euler-kocks`, each orientation is described by :math:`\Psi, \Theta, \phi`, where :math:`\Psi` is the rotation about the :math:`z` axis, :math:`\Theta` is the rotation about the :math:`y^{\prime}` axis, and :math:`\phi` is the rotation about the :math:`z^{\prime \prime}` axis, all in degrees). The per-line format is::

    <entity_id> <Psi> <Theta> <phi>

- For :data:`axis-angle`, each orientation is described by :math:`t_1, t_2, t_3, \omega`, where :math:`\bf{t}` is the normalized axis of rotation and :math:`\omega` is the angle of rotation about said axis, in degrees. The per-line format is::

    <entity_id> <t_1> <t_2> <t_3> <omega>

- For :data:`quaternion`, each orientation is described by :math:`q_0, q_1, q_2, q_3`, where :math:`q_0 = \cos{(\omega / 2)}` and :math:`q_i = t_i \sin{(\omega / 2)}` for :math:`i = 1, 2, 3` (see: :data:`axis-angle` for definitions of :math:`{\bf t}` and :math:`\omega`). The per-line format is::

    <entity_id> <q_0> <q_1> <q_2> <q_3>

.. _phase_file:

Phase File (:file:`simulation.phase`)
-------------------------------------

The embedded grain/phase assignments within the :file:`simulation.msh` may be overridden via an external :file:`simulation.phase` file. This file contains formatting identical to the associated fields in the mesh file and is defined as::

    $Groups
    <group_entity>
    <number_of_group_entities>
    <entity_id> <group>
    ...
    $EndGroups

where :data:`<group_entity>` defines the phase assignment method and must always be defined as :data:`elset`, :data:`<number_of_group_entities>` is the number of grain/phase pairs defined in the field, :data:`<entity_id>` is a unique 1-indexed identification number, and :data:`<group>` is an 1-indexed value that defines the phase for a given grain.

.. _restart_output:

Restart Files (:file:`rst<N>.*`)
--------------------------------

If the :data:`print restart` command is present in the :file:`simulation.cfg` file, a set of additional restart files will be generated from the simulation. These files are written at the end of each prescribed step and contain necessary information to restart a given simulation (:ref:`restart_tutorial` for information on how to restart a simulation). Two types of restart files are generated, a control file, :file:`rst<N>.control`, and per-core field files, :file:`rst<N>.field.core*` (where :data:`<N>` indicates which simulation the files describe, 0 indexing). Both file types are unformatted (or binary) files and are generally unmodifiable. The structures of the data stored within both files for the various deformation modes follow.

Uniaxial Restart Control
~~~~~~~~~~~~~~~~~~~~~~~~

The :file:`rst<N>.control` file for uniaxial loading modes contains the following data in the given order::

    current_step <step>
    previous_load_array <load_x> <load_y> <load_z>
    step_complete_flag <logical>
    previous_timestep_value <time>
    current_incr <increment>
    current_time <time>
    surface_1_loads <load_x> <load_y> <load_z>
    ...
    surface_6_loads <load_x> <load_y> <load_z>
    previous_prescribed_load <load>
    current_surface_areas <area_surf_1> ... <area_surf_6>
    initial_surface_areas <area_surf_1> ... <area_surf_6>

Multiaxial CSR Restart Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :file:`rst<N>.control` file for multiaxial constant strain rate loading modes contains the following data in the given order::

    current_step <step>
    current_load_array <load_x> <load_y> <load_z>
    previous_load_array <load_x> <load_y> <load_z>
    step_complete_flag <logical>
    previous_timestep_value <time>
    current_incr <increment>
    current_time <time>
    surface_1_loads <load_x> <load_y> <load_z>
    ...
    surface_6_loads <load_x> <load_y> <load_z>
    current_surface_areas <area_surf_1> ... <area_surf_6>
    initial_surface_areas <area_surf_1> ... <area_surf_6>
    current_mesh_lengths <length_x> <length_y> <length_z>
    initial_mesh_lengths <length_x> <length_y> <length_z>
    current_control_velocity <vel_x> <vel_y> <vel_z>
    s_pert_mag <vel>
    t_pert_mag <vel>

Multiaxial CLR Restart Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :file:`rst<N>.control` file for multiaxial constant load rate loading modes contains the following data in the given order::

    current_step <step>
    current_load_array <load_x> <load_y> <load_z>
    previous_load_array <load_x> <load_y> <load_z>
    first_incr_in_step <logical>
    current_incr <increment>
    current_time <time>
    surface_1_loads <load_x> <load_y> <load_z>
    ...
    surface_6_loads <load_x> <load_y> <load_z>
    current_surface_areas <area_surf_1> ... <area_surf_6>
    initial_surface_areas <area_surf_1> ... <area_surf_6>
    current_mesh_lengths <length_x> <length_y> <length_z>
    initial_mesh_lengths <length_x> <length_y> <length_z>
    current_control_velocity <vel_x> <vel_y> <vel_z>
    previous_control_action <integer>
    current_control_action <integer>
    initial_load_dwell_velocity <vel_x> <vel_y> <vel_z>
    initial_unload_dwell_velocity <vel_x> <vel_y> <vel_z>

Restart Field Data
~~~~~~~~~~~~~~~~~~

All loading modes also write field data on a per-core basis to :file:`rst<N>.field.core*` files. These files contain the necessary field variable information in order to spatially define the total state of the virtual sample at the time of printing. The following field data arrays are written to the files in the given order::

    coords <coords>
    velocity <velocity>

    c0_angs <orientation>
    c_angs <orientation>
    rstar <rotation>
    rstar_n <rotation>
    wts <weight>
    crss <crss>
    crss_n <crss>

    gela_kk_bar <strain>
    gsig_vec_n <stress>
    pela_kk_bar <strain>
    psig_vec_n <stress>
    e_elas_kk_bar <strain>
    sig_vec_n <stress>

    eqstrain <strain>
    eqplstrain <strain>
    gamma <shear>

    el_work_n <work>
    el_workp_n <work>
    el_work_rate_n <work_rate>
    el_workp_rate_n <work_rate>

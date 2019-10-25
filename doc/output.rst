.. _output:

Understanding and Analysing |tf| Output
#######################################

|tf| produces several different types of trees. A descendant tree, a progenitor tree,
and a cross catalog. The output is set by the runtime mode of operation, whether
the code is building a descendant tree, progenitor tree or simple comparing two input catalogs.

When run with MPI, each thread will write the snapshots that have been processed by that mpi thread
when writing HDF output. For ASCII output, either a single continuous ascii file is written or
each thread writes a file, adding the rank of the mpi thread writing the file.

.. topic:: Standard files

    * ``.snapshot_%03d.VELOCIraptor.tree``: a tree file


TreeFile
========

The exact format of a tree file depends on whether the code produces a descendant tree
or progenitor tree. The information contain also depends on the output format. As the suggetsion
is to produce HDF format unless otherwise required we only list the HDF output in detail.

+------------------------------+---------------------------------------------------------------------------------------------+
| Name                         | Comments                                                                                    |
+==============================+=============================================================================================+
| `Header Attributes`                                                                                                        |
+------------------------------+---------------------------------------------------------------------------------------------+
| Number_of_snapshots          |  Number of snapshots in the tree                                                            |
+------------------------------+---------------------------------------------------------------------------------------------+
| Total_number_of_halos        |  Total number of halos across all snapshots                                                 |
+------------------------------+---------------------------------------------------------------------------------------------+
| Merit_limit                  |  Merit limit used to determine whether a connection is viable to be a primary connection    |
+------------------------------+---------------------------------------------------------------------------------------------+
| Number_of_steps              |  Number of snapshots searched for primary connections                                       |
+------------------------------+---------------------------------------------------------------------------------------------+
| Search_next_step_criterion   |  Integer indicating type of criterion used to keep searching for primary connection         |
+------------------------------+---------------------------------------------------------------------------------------------+
| Merit_limit_for_next_step    |  Merit limit below which more snaphots are search for viable primary connection             |
+------------------------------+---------------------------------------------------------------------------------------------+
| Core_fraction                |  Fraction of most bound particles used to calculate merits                                  |
+------------------------------+---------------------------------------------------------------------------------------------+
| Core_min_number_of_particles |  Minumum number of most bound particles used to calculate merits                            |
+------------------------------+---------------------------------------------------------------------------------------------+
| Description                  |  String describing how tree was produced                                                    |
+------------------------------+---------------------------------------------------------------------------------------------+
| `Tree Data with arrays typically the size of number of halos in given snapshot`                                            |
+------------------------------+---------------------------------------------------------------------------------------------+
| ID                           |  Tree Halo IDs (index of halo + 1 + TEMPORALHALOIDVAL * Snapshot_value)                     |
+------------------------------+---------------------------------------------------------------------------------------------+
| OrigID                       |  Original ID in halo catalog (IF compiled with HALOIDNOTINDEX. Otherwise not present)       |
+------------------------------+---------------------------------------------------------------------------------------------+
| Npart                        |  Number of particles in a halo. Only produced if desired.                                   |
+------------------------------+---------------------------------------------------------------------------------------------+
| NumDescen/NumProgen          |  Number of descendants/progenitors                                                          |
+------------------------------+---------------------------------------------------------------------------------------------+
| DescenOffsets/ProgenOffsets  |  An offset array indicating where a halo's connections begin associated descen/progen array |
+------------------------------+---------------------------------------------------------------------------------------------+
| `Tree Data with arrays the size of the number of viable connections found. This can be larger than the number of halos`    |
+------------------------------+---------------------------------------------------------------------------------------------+
| Descendants/Progenitors      |  Array of Tree ID connections. Halos can have 0,1,>1 connections. Array read using          |
|                              |  the NumDescen and DescenOffsets arrays                                                     |
+------------------------------+---------------------------------------------------------------------------------------------+
| Merit                        |  Merit of the connection                                                                    |
+------------------------------+---------------------------------------------------------------------------------------------+
| DescenNpart/ProgenNpart      |  Number of particles in descendant/progenitor. Only produced if desired.                    |
+------------------------------+---------------------------------------------------------------------------------------------+

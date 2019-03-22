.. _usage:

Using **TreeFrog**
######################

.. _running:

Running the code
================

**TreeFrog** is a stand alone executable (with the executable named stf (or STructure Finder for historical reasons)).
It can be run in serial, with OpenMP, or MPI APIs. A typical command to start the code looks like:
::

 ./treefrog < args >

When compiled with OpenMP, setting the environment variable ``OMP_NUM_THREADS`` will set the number of threads in the openmp sections.

With MPI using 8 MPI threads:
::

 mpirun -np 8 ./treefrog < args >

where here we assume that the parallel
environment uses the ``mpirun`` command to start MPI
applications. Depending on the operating system, other commands may be
required for this task, e.g. ``srun`` on some Cray machines. Note that
the code can in principle be started using an arbitrary number of
mpi threads, but the mpi decomposition is most efficient for powers of 2.

The output produced by TreeFrog will typically consist of several files containing:
bulk properties of structures found; particles belonging to these structures; and several
additional files containing configuration information.

When running in MPI, currently each mpi thread writes its own output.

.. _cmdargs:

Arguments
---------

The code has several command line arguements. To list the arguments, type
::

    ./treefrog -?

.. topic:: The main arguments that can be passed are:

    **-i** ``< file name of input file that contains list of files to be processed >``

    **-s** ``< number of files (snaphots) >``

    **-I** ``< input format [2 VELOCIraptor, 1 SUSSING, 3 nIFTy, 4 Void] >``

    **-o** ``< output base name >``

    **-C** ``< configuration file name (see`` :ref:`configoptions` ``) >``

Of these arguements, only an input file, number of files (snapshots) in input
and an output name must be provided.
In such a case, default values for all other confirguration options are used.
We suggest you do NOT run the code in this fashion.
Instead we suggest the code be run with at least a configuration file passed.
::

    ./treefrog -i input -s numsnaps -o output -C configfile.txt

This configuration file is an ascii file that lists keywords and values.
A list of keywords, along with a description is presented below in :ref:`configoptions`.
A more typical command for a large cosmological simulation might be something like
::

    export OMP_NUM_THREADS=4
    mpirun -np 2 ./treefrog -i listofsnaphots -s 128 -o halotree -C configfile.txt > treefrog.log

.. _configoptions:

Configuration File
------------------

An example configuration file can be found the examples directory within the repository
(see for instance :download:`sample <../examples/treefrog_sample.configuration>`). This sample file lists
all the options. *Only the keywords listed here will be used, all other words/characters
are ignored*.

.. warning:: Note that if misspell a keyword it will not be used.


There are numerous key words that can be passed. Here we list them, grouped into several categories:
:ref:`Outputs <config_output>`,
:ref:`Inputs <config_input>`,
:ref:`Parameters related to type of search <config_search_type>`,
:ref:`Miscellaneous <config_misc>`,
:ref:`MPI <config_mpi>`.

.. _config_output:

.. topic:: Output related
    ``Output format = 2/0``
        * Flag indicating whether output is 2 HDF, 0 ASCII
    ``Output data included = 1/0/2``
        * Flag indicating whether to produce standard output, minial output, or extensive output. Standard includes merits. Extensive includes merits and nubmer of particles

.. _config_input:

.. topic:: Input related

    ``Input_tree_format = 2/1/3/4``
        * Type of input halo catalog. 2 is VELOCIraptor input, 1 SUSSING, 3 nIFTY, 4 Void.
    ``VELOCIraptor_input_format = 2/0/1``
        * Input format of a VELOCIraptor catalog, 2 HDF, 0 ASCII, 1 binary.
    ``VELOCIraptor_input_field_sep_files = 0/1``
        * Whether VELOCIraptor catalog has separate files for field halos and subhalos.
    ``VELOCIraptor_input_num_files_per_snap = 0/1``
        * Whether there is more than one file per VELOCIraptor catalog (if it was run in MPI mode)

.. _particle_id_options:

.. topic:: Parameters related Particle IDs (which are used to cross correlate catalogs)
    ``Max_ID_Value = ``
        * TreeFrog allocates array of size Max_ID_Value to cross correlate particles thus specify maximim ID and code will allocate an array of size max ID of either ints or long ints (depending on compilation options) to cross correlate. If value not set must set an id to index mapping.
    ``Mapping= 0/1/-1 ``
        * Can construct a memory efficient ID to index map (computationally expensive but reduces) memory by prodiving a map or by having code produce a map. No map 0, memory efficient treefrog built map -1, user defined (must alter code) map 1

.. _config_search_type:

.. topic:: Parameters related to type of search

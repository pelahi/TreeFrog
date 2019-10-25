.. _usage:

Using |tf|
##########

.. _running:

Running the code
================

|tf| is a stand alone executable. It can be run in serial, with OpenMP, or MPI APIs.
A typical command to start the code looks like:
::

 ./treefrog < args >

When compiled with OpenMP, setting the environment variable ``OMP_NUM_THREADS`` will set the number of threads in the openmp sections.

With MPI using 8 MPI threads:
::

 mpirun -np 8 ./treefrog < args >

where here we assume that the parallel
environment uses the ``mpirun`` command to start MPI
applications. Depending on the operating system, other commands may be
required for this task, e.g. ``srun`` on some Cray machines.

The output produced by |tf| will consist of either a single ASCII file or multiple
HDF5 files (on per snapshot). The information contained in it this output is a
raw halo merger tree listing the connections halos have at one snapshot in time
to later (or earlier) snapshots. More details of the output can be found in :ref:`output`.

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
and an output name must be provided. In such a case, default values for all other confirguration options are used.
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


By default the code is works natively with output from the **VELOCIraptor** halo finder
(which can be obtained via `<https://github.com/pelahi/VELOCIraptor-STF>`_, and which has associated
`online documetation <https://velociraptor-stf.readthedocs.io/en/latest/>`_).

The code is able to process other types of halo finder input (see :ref:`AHF`
for a description of running |tf| on output from AHF).

Finally, |tf| can be used to produce halo merger trees used as input for semi-analytic
models of galaxy formation (see :ref:`sammergertree`).


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
:ref:`Particle ID options <config_particle_id_options>`,
:ref:`Tree construction options <config_tree_construction>`,
:ref:`Merit options <config_merit>`,
:ref:`Temporal linking options <config_temporal_linking>`,
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

.. _config_particle_id_options:

.. topic:: Particle IDs options

    How to handle particle ids which are used to cross correlate catalogs.

    ``Max_ID_Value = 1073741824``
        * TreeFrog allocates array of size Max_ID_Value to cross correlate particles thus specify maximim ID and code will allocate an array of size max ID of either ints or long ints (depending on compilation options) to cross correlate. If value not set must set an id to index mapping.
    ``Mapping= 0/1/-1``
        * Can construct a memory efficient ID to index map (computationally expensive but reduces) memory by prodiving a map or by having code produce a map. No map 0, memory efficient treefrog built map -1, user defined (must alter code) map 1

.. _config_tree_construction:

.. topic:: Tree construction options

    ``Tree_direction = 1``
        * Integer indicating direction in which to process snapshots and build the tree. Descendant [1], Progenitor [0], or Both [-1].
    ``Particle_type_to_use = -1``
        * Integer describing particle types to use when calculating merits. All [-1], Gas [0], Dark Matter [1], Star [4].
    ``Default_values = 1``
        * Whether to use default cross matching & merit options when building the tree. 1/0 for True/False.

.. _config_merit:

.. topic:: Merit options

    Related to what merit function to use to define matches

    ``Merit_type = 6``
        * Integer specifying merit function to use. Several options available using variations of two specific merits:
            * the shared number of particles :math:`\mathcal{N}_{A_{i}B_{j}} = N_{A_{i}\bigcap B_{j}}^2/(N_{A_{i}}N_{B_{j}})`
            * the rank ordering of particles :math:`\mathcal{S}_{A_{i}B_{j},A_{i}} = \sum_{l}^{N_{A_{i}\bigcap B_{j}}} 1/\mathcal{R}_{l,A_{i}}`
        * Optimal descendant tree merit is combination of both rank ordered in both directions [6], common (progenitor tree) merit in is using the shared merit [1].
    ``Core_match_type = 2``
        * Integer flag indicating the type of core matching used. Off [0], core-to-all [1], core-to-all followed by core-to-core [2], core-to-core only [3].
    ``Particle_core_fraction = 0.4``
        * Fraction of particles to use when calculating merits. Assumes some meaningful rank ordering to input particle lists and uses the first :math:`f_{\rm TF}` fraction.
    ``Particle_core\_min_numpart = 5``
        * Minimum number of particles to use when calculating merit if core fraction matching enabled.

.. _config_temporal_linking:

.. topic:: Temporal linking options

    Related to how code searches for candidate links across multiple snapshots.

    ``Nsteps_search_new_links = 1``
        * Number of snapshots to search for links.
    ``Multistep_linking_criterion = 3``
        * Integer specifying the criteria used when deciding whether more snapshots should be searched for candidate links. Criteria depend on tree direction.
            * **Descendant Tree**: continue searching if halo is: missing descendant [0]; missing descendant or descendant merit is low [1]; missing descendant or missing primary descendant [2]; missing a descendant, a primary descendant or primary descendant has poor merit [3].
            * **Progenitor tree**: continue searchign if halo is: missing progenitor[0]; missing progenitor or progenitor merit is low [1].
    ``Merit_limit_continuing_search = 0.025``
        * Float specifying the merit limit a match must meed if using ``Multistep_linking_criterion = 1 (progenitor) or 3 (descendant)``.

.. _config_mpi:

.. topic:: MPI related options

    Related to MPI options

.. _config_misc:

.. topic:: Miscellaneous options

    Miscellaneous options

    ``Verbose = 0/1/2``
        * Indicates how verbose the code is while running. 0 is minimial, 1 verbose, 2 chatterbox.

.. _othercats:

Running on Other Catalogs
=========================

.. _AHF:

**AHF** catalogue
-----------------

To process the ASCII output produced by the `AHF halo finder <http://popia.ft.uam.es/AHF/Download.html>`_, the mpi flag needs to be switched off and the flag that tells TreeFrog the ID's do not correspond to an index (as with AHF halo ID's) needs to be switched on. These flags are
::

	TF_MPI:BOOL=ON
	TF_HALOIDNOTINDEX:BOOL=OFF

There are several configuration options that must be set. The input format must be set appropriately.
::

	Input_tree_format = 3

The next one is the maximum particle id value that should be set to the total amount of all particles in the simulation.
::

	Max_ID_Value


The code can run with
::

	./treefrog -C ../treefrog.cfg -i <filelist> -s <numsnapshots> -o <baseoutputfilename>

Where `<filelist>` is a file containing the _particles files for each snapshot from the **AHF** output.

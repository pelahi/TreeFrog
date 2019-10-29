.. _getting:

Getting |tf|
############

|tf| is currently hosted in `GitHub <https://github.com/pelahi/TreeFrog>`_.
To get a copy you can clone the repository::

    git clone https://github.com/pelahi/TreeFrog

|tf|'s compilation system is based on `cmake <https://www.cmake.org/>`_. ``cmake`` will
check that you have a proper compiler (anything supporting C++14 or later should do),
and scan the system for all required dependencies.

To compile |tf| run (assuming you are inside the ``TreeFrog`` directory already)::

 $> mkdir build
 $> cd build
 $> cmake ..
 $> make

With ``cmake`` you can also specify additional compilation flags.
For example, if you want to generate the fastest possible code
you can try this::

 $> cmake .. -DCMAKE_CXX_FLAGS="-O3 -march=native"

You can also specify a different installation directory like this::

 $> cmake .. -DCMAKE_INSTALL_PREFIX=~/my/installation/directory

Other ``cmake`` options that can be given in the command-line include:

A list of compile time options is found below in :ref:`compileoptions`.

.. _requirements:

Requirments
===========

|tf| depends on:

* `GSL <https://www.gnu.org/software/gsl/>`_ - the GNU Scientific Library
* **NBodylib** - a internal scientific library included with |tf| as a submodule.

Optional requirements
---------------------

For parallel use may need the following libraries are required for compilation
depending on the compilation flags used:

* MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g. `Open MPI <https://www.open-mpi.org/>`_, `MPICH <http://www-unix.mcs.anl.gov/mpi/mpich/>`_ or
  `LAM <http://www.lam-mpi.org/>`_.

* `OpenMP <http://www.openmp.org/>`_ - API, generally included with many compilers

|tf| also can output in a variety of formats: ASCII, and HDF.
HDF can be enabled and disabled, and requires libraries.

* `Hiearchical Data Format (HDF) <https://www.hdfgroup.org/>`_ - self describing data format.

.. _compileoptions:

Compilation Options
===================

These can be passed to ``cmake``

.. topic:: External library flags

    * Parallel APIs can be enabled by setting
        * For MPI
            | ``TF_MPI``: boolean to compile with MPI support
            | ``MPI_LIBRARY``: specify library path to MPI
            | ``MPI_EXTRA_LIBRARY``: Extra MPI libraries to link against
        * For OpenMP
            | ``TF_OPENMP``: boolean to compile with OpenMP support
            | ``OpenMP_CXX_FLAGS``: string, compiler flag that enables OpenMP


    * Enable input/output formats
        * For HDF
            | ``TF_HDF5``: boolean on whether to include HDF support
            | ``HDF5_ROOT``: specify a local directory containing HDF library.

    * To set directories of required libraries
        * Set the directories of the following libraries
            | ``GSL_DIR =``

.. topic:: Internal precision and data structure flags

    * Adjust precision in stored variables and calculations
        * all integers are 64 bit integers. Enable this if dealing with more than MAXINT total number of particles
            ``TF_LONG_INT``: boolean on whether to use long ints. Default is ON
        * Use unsigned ids.
            ``TF_UNSIGNED_IDS``: boolean on whether to use unsigned (long) ints

.. topic:: Operation flags

    * Adjust **TreeFrog** operation
        * Assume halo indices and snapshot do not map to a halo ID and halo ID must be read from input.
            ``TF_HALOIDNOTINDEX``: boolean on whether ids and indices map. Default is OFF

.. topic:: Executable flags

    * Enable debugging
        ``DEBUG``: boolean on whether to run with debug flags on and no optimisation.

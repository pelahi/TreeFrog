.. _sammergertree:

Producing SAM digestible merger trees
#####################################

One of the key uses of halo merger trees is to produce synthetic galaxies
using semi-analytic models such as `shark <https://github.com/ICRAR/shark>`_,
`sage <https://github.com/darrencroton/sage>`_ or
`meraxes <https://www.ph.unimelb.edu.au/~smutch/papers/meraxes/meraxes.html>`_.

We provide some background to the production of halo merger trees. Readers simply
interested in producing an output can skip the background and simply follow the
steps described in :ref:`samhmtwalkabletree` and subsequent sections
or just follow the steps listed in :ref:`samhmtsummarysteps`.

We note that the default mode of |tf| is to process snapshots walking forward in
time, identifying descendants but it can operate in a mode walking backward in
time identifying progenitors. Halo merger trees for SAMs are better built identifying
links by walking forward in time. Thus we focus here on using `Descendant Trees` but
it is possible to use `Progenitor Trees`.

.. _samhmtbackground:

Background
==========

The raw output from |tf| is more akin to a graph than a simple halo merger tree.
The file contains for all halos the optimal connection and also secondary connections.
This extra information is useful for exploring a variety of topics, from the performance
of halo finders to how particles can be exchanged between mergers. However, this
extra complexity is not necesary for SAMs.

The raw tree also does not list for a given halo, the first or progenitor or last
descendant but instead just the immediate connection. Most SAM models also make
use of this information.

Currently, |tf| makes use of python tools available in the tools directory
to produce input that is digestible by such SAMs, converting the raw tree to
one which is easily naviable by such codes. There are also scripts availbe in the
examples directory that can be used.

.. note::

   This process may change so that |tf| natively produces a raw tree along with
   the halo merger trees used by a SAM.

.. _samhmtwalkabletree:

Generating Walkable Tree
========================

To produce information that stores for each halo, its progenitor (Tail),
descendant (Head), first progenitor (Root Tail), and final descendant (Root Head)
we make use of the function

.. code-block:: python

    #for processing a standard descendant based tree frog tree
    velociraptor_python_tools.BuildTemporalHeadTailDescendant(numsnaps: int,
        tree: list, numhalos: np.array, halodata: list,
        TEMPORALHALOIDVAL : long = 1000000000000, ireverseorder : bool = False,
        iverbose : int = 1)
    #for processing a standard progenitor based tree
    velociraptor_python_tools.BuildTemporalHeadTail(numsnaps: int,
        tree: list, numhalos: np.array, halodata: list,
        TEMPORALHALOIDVAL : long = 1000000000000, iverbose : int = 1)

This code uses the raw |tf| data and the information from halo catalogues to
produce a simple wakable tree. An example script is located
:download:`here <../examples/example_produce_walkabletree.py>`.

The process is as follows:

.. code-block:: python

    #load the velociraptor python tools
    import sys
    sys.path.append(TreeFrog/tools/)
    import velociraptor_python_tools as vpt

    # to determine how many snapshots there are, how many snapshots
    # were used to identify links and the temporal halo ID,
    # we can use the firt file from TreeFrog
    numsnaps = 100
    # how do ids map to snapshot, ie. the temporal halo id.
    TEMPORALHALOIDVAL = 1000000000000L

    # load tree data. Here we assume data is in HDF5 format and follows
    # a specific naming convention
    basetreename = 'treefrog'
    # produce a file listing the treefrog output
    snaptreelist=open('snaptreelist.txt','w')
    for i in range(numsnaps):
        snaptreelist.write(basetreefname+'.snapshot_%03d.VELOCIraptor\n'%i)
    snaptreelist.close()

    # read a treefrog descendant tree
    reverseorderflag = False
    formatflag = 2 #HDF5
    iverbose = 0 #not verbose
    meritinfoflag = True #merit information present
    treefrogdata = vpt.ReadHaloMergerTreeDescendant('snaptreelist.txt',
        reverseorderflag, formatflag, iverbose, meritinfoflag)

    # read halo catalog. Here we assume that VELOCIraptor has been used
    # we also assume a specific naming convention
    # allocate data structures to store the information
    numhalos=np.zeros(numsnaps, dtype=np.int64)
    halodata = [None for i in range(numsnaps)]
    scalefactors = np.zeros(numsnaps)
    # as we do not need all the information in the halo catalogs
    # only request a subset of the fields
    requestedfields = ['ID', 'hostHaloID']

    # load halo properties file (this also assumes HDF input)
    iverbose = 0
    separatefilesforhaloandsubhalos = 0
    for i in range(numsnaps):
        fname='snapshot_%03d.VELOCIraptor'%i
        halodata[i],numhalos[i] = vpt.ReadPropertyFile(fname, formatflag,
            separatefilesforhaloandsubhalos, iverbose, requestedfields)
        scalefactors[i] = halodata[i]['SimulationInfo']['ScaleFactor']

This loads all the data necessary to make a walkable tree

.. code-block:: python

    #build the walkable tree
    vpt.BuildTemporalHeadTailDescendant(numsnaps,
        treefrogdata, numhalos, halodata, )

We now save the data

.. code-block:: python

    # We store the information related to
    # how the tree was built in a dictonary.
    # here values are hard coded but can be taken
    # from the input.
    DescriptionInfo={
            'Title':'Walkable Tree',
            'TreeBuilder':'TreeFrog',
            'TreeBuilder_version':1.20,
            'Temporal_linking_length':NSNAPSEARCH,
            'Temporal_halo_id_value':TEMPORALHALOIDVAL,
            'HaloFinder':'VELOCIraptor',
            'HaloFinder_version':1.11,
            'Particle_num_threshold':20,
            }
    # write file
    outputfname = 'walkablehalomergertree.hdf5'
    vpt.WriteWalkableHDFTree(outputfname, numsnaps, treefrogdata,
        numhalos, halodata, scalefactors, DescriptionInfo)

Using the :download:`script <../examples/example_produce_walkabletree.py>`
simply requires altering it to the desired naming convention and running it.

.. code-block:: shell

    #we set the appropriate variables
    base_treefrog_filename='treedir/treefrog'
    #base halo catalog where we assume the names are in directory and follow
    #a specific naming convention
    halocatalog_dir='halos'
    output_filename='treedir/walkabletree.hdf5'
    python3 ${base_treefrog_filename} ${halocatalog_dir} ${output_filename}

.. _samhmtforshark:

Generating Input for **shark**
------------------------------

The semi-analytic code shark is designed to load the this walkable tree and
the halo catalogues. No further process of |tf| is required.

.. _samhmtforest:

Generating Forest
=================

Some SAMs require more information to process output. This can range from just
extra links to quickly navigate halo catalogs. We focus here on producing output
that contains not only a halo's progneitor and descendant ID but also a Forest ID.
The idea of a halo forest is a collection of halo merger trees that have interacted
with each other at some point in cosmic time. The interaction is typically taken to
be that a halo has become a subhalo of another halo at some point. However,
such a concept can be generalised to halos that entire some factor of the virial
radius of another halo. Here we limit the forest to objects that have become
subhalos of another halo as defined by the FOF envelop.

.. _samhmtforsage:

Generating Input for **sage**
-----------------------------

.. _samhmtformeraxes:

Generating Input for **meraxes**
--------------------------------

.. _samhmtsummarysteps:

Summary of steps
================

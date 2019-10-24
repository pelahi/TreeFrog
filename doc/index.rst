.. TreeFrog documentation master file, created by
   sphinx-quickstart on Mon Jul 31 10:13:40 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TreeFrog
########

.. image:: treefrog.png
  :scale:  75 %
  :align:  left

|tf| is a C++ halo merger tree builder using MPI and OpenMP APIs.
The repository also contains several associated analysis tools in python, example configuration files
and analysis python scripts (and sample jupyter notebooks).

If you are using |tf| please cite the following paper, which describe the code in full::

    @ARTICLE{10.1017/pasa.2019.18,
        author = {{Elahi}, Pascal J. and {Poulton}, Rhys J.~J. and {Tobar}, Rodrigo J. and
        {Ca{\~n}as}, Rodrigo and {Lagos}, Claudia del P. and {Power}, Chris and
        {Robotham}, Aaron S.~G.},
        title = "{Climbing halo merger trees with TreeFrog}",
        journal = {\pasa},
        keywords = {dark matter, methods: numerical, galaxies: evolution, galaxies: halos, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Astrophysics of Galaxies},
        year = "2019",
        month = "Jan",
        volume = {36},
        eid = {e028},
        pages = {e028},
        doi = {10.1017/pasa.2019.18},
        archivePrefix = {arXiv},
        eprint = {1902.01527},
        primaryClass = {astro-ph.IM},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2019PASA...36...28E},
    }

An online entry can also be found at `NASA's ADS service <https://ui.adsabs.harvard.edu/abs/2019PASA...36...28E/abstract>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting
   usage
   output
   sammergertree
   dev

.. _pascaljelahi@gmail.com: mailto:pascaljelahi@gmail.com

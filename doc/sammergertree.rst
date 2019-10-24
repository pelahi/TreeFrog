.. _sammergertree:

Producing SAM digestible merger trees
#####################################

One of the key uses of halo merger trees is to produce synthetic galaxies
using semi-analytic models such as `shark <https://github.com/ICRAR/shark>`_,
`sage <https://github.com/darrencroton/sage>`_ or
`meraxes <https://www.ph.unimelb.edu.au/~smutch/papers/meraxes/meraxes.html>`_.

The input to such codes is more than a raw tree with immediate connections between
halos at one snapshot to halos at a later snapshot. Currently, |tf| makes use of
python tools to produce input that is digestible by such SAMs, converting the raw
tree to one which is easily naviable by such codes. This process may change so that
|tf| natively produces a raw tree (which is actually more akin to a graph)
along with the paired down halo merger trees. We outline the process to convert
a raw tree input to a halo merger tree here.

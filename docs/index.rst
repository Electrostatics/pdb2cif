.. old_pdb documentation master file, created by
   sphinx-quickstart on Sun Feb 21 19:43:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

old_pdb
=======

This code reads and writes the `old PDB format <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_.

Development of this code was greatly helped by `the Chimera documentation <https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html>`_.

-----------------------
Origin of this software
-----------------------

This software is derived, in part, from the `PDB2PQR code <http://pdb2pqr.readthedocs.io>`_ which is supported by `NIH <htp://www.nih.gov>`_ grant GM069702.

------------
Installation
------------

This python package can be installed via `setuptools <https://pypi.org/project/setuptools/>`_, ``pip install .``

-------
Testing
-------

The software can be tested with `pytest <https://docs.pytest.org/en/stable/>`_ by running::

   python -m pytest

from the top-level directory.

------------------
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`

--------
Contents
--------

.. toctree::
   :maxdepth: 2

   api/index
   changelog

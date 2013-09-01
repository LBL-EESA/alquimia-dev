.. alquimia documentation master file, created by
   sphinx-quickstart on Wed Jan 23 20:43:25 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

..
   Alquimia Copyright (c) 2013, The Regents of the University of California, 
   through Lawrence Berkeley National Laboratory (subject to receipt of any 
   required approvals from the U.S. Dept. of Energy).  All rights reserved.
   
   Alquimia is available under a BSD license. See LICENSE.txt for more
   information.
   
   If you have questions about your rights to use or distribute this software, 
   please contact Berkeley Lab's Technology Transfer and Intellectual Property 
   Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
   
   NOTICE.  This software was developed under funding from the U.S. Department 
   of Energy.  As such, the U.S. Government has been granted for itself and 
   others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
   license in the Software to reproduce, prepare derivative works, and perform 
   publicly and display publicly.  Beginning five (5) years after the date 
   permission to assert copyright is obtained from the U.S. Department of Energy, 
   and subject to any subsequent five (5) year renewals, the U.S. Government is 
   granted for itself and others acting on its behalf a paid-up, nonexclusive, 
   irrevocable, worldwide license in the Software to reproduce, prepare derivative
   works, distribute copies to the public, perform publicly and display publicly, 
   and to permit others to do so.
   
   Authors: Benjamin Andre <bandre@lbl.gov>
    
    

Welcome to alquimia's documentation!
====================================

Alquimia will use `Semantic Versioning <http://semver.org/>`_ for its public API.

The Alquimia API and library will use the `Google C++ style guide
<http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml>`_
where appropriate. In particular, the naming conventions should be
followed consistently on both sides (driver and engine) of the
interface wrapper code. Engines are obviously allowed to use their own
internal coding styles in their own code.

Contents:

.. toctree::
   :maxdepth: 2

   api/APIv0
   api/APIv0_functions
   api/APIv0_structures
   api/APIv0_constants
   api/APIv0_c_utils

   batch-tests

   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Glossary
========

.. glossary::

    driver
        the transport simulator or other driver, e.g. Amanzi-u, Amanzi-s, ParFlow

    engine
        the backend geochemistry engine, e.g. PFloTran, CrunchFlow, toughreact, stomp, PHREEQC, ...


.. alquimia documentation master file, created by
   sphinx-quickstart on Wed Jan 23 20:43:25 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

..
   Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
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

Description
-----------

Alquimia is an biogeochemistry API and wrapper library. The aim is to
provide a unified interface to the biogeochemistry routines from high
quality reactive transport simulators such as `PFLOTRAN
<https://bitbucket.org/pflotran/pflotran-dev>`_ or `CrunchFlow
<http://www.csteefel.com/CrunchFlowIntroduction.html>`_ , allowing any
subsurface flow and transport simulator to access a range of
functionality.

It is not an implementation of a biogeochemistry reaction library, and
does not do any geochemical calculations.

Availability
------------

The official repository for Alquimia API documentation and reference
API is `<https://bitbucket.org/berkeleylab/alquimia>`_.

The alquimia API is intended to be open and implemented by anyone. The
reference implementation developed by `Lawrence Berkeley National
Laboratory <http://www.lbl.gov>`_ is open source, licensed under a BSD
:ref:`license`.

Versioning
----------

Alquimia API and reference implementation use `Semantic Versioning
<http://semver.org/>`_.

Style Guide
-----------

The Alquimia API and library will use the `Google C++ style guide
<http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml>`_
where appropriate. In particular, the naming conventions should be
followed consistently on both sides (driver and engine) of the
interface wrapper code. Engines are obviously allowed to use their own
internal coding styles in their own code.

Guidance for use
----------------

.. toctree::
   :maxdepth: 2

   guidance

Contents
--------

.. toctree::
   :maxdepth: 2

   api/APIv1
   api/APIv1_functions
   api/APIv1_structures
   api/APIv1_constants
   api/APIv1_c_utils

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


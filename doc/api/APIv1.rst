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


API
===

Alquimia will use `Semantic Versioning <http://semver.org/>`_ for its public API.

Required parts of the interface API include:

* :doc:`function <APIv1_functions>` call signatures
* :doc:`data structures <APIv1_structures>`
* :doc:`constants <APIv1_constants>`

The required parts of the interface will be compiled into
``libalquimia_c.a`` and ``libalquimia_fortran.a``.

There is also an optional C utility library containing reuseable code
for common tasks (allocating memory, printing data)

* :doc:`c utility library <APIv1_c_utils>`



Required Functionality
~~~~~~~~~~~~~~~~~~~~~~

Alquimia exposes the following functionality from the underlying
geochemistry :term:`engine`\ s in a generic interface:

* Reading of geochemical reaction data

* Basis management (includes reading database, swapping basis, etc.)

* Constraint processing for boundary/initial constraints. (Say you set
  a boundary concentration to be equilibrated with a mineral. The
  total component concentration for the species equilibrated with the
  mineral needs to be calculated and sent back to Amanzi transport.)

* Reaction step: operator split and global implicit

* Access to user selected geochemical data for output, i.e. pH, mineral SI, reaction rates

* Additional information about the library, e.g. is it thread safe, so
  amanzi can create multiple copies using openmp? Does it support
  temperature dependent chemistry, porosity updates, etc.

Alquimia has two parts.
    * An engine independent API consisting of :doc:`function
      <APIv1_functions>` call signatures and :doc:`data structures
      <APIv1_structures>`.
    * An optional :doc:`utility <APIv1_c_utils>` library to handle
      data memory allocation/freeing, printing structs, etc.


Division of labor / responsibilities for implementing alquimia
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Driver Process Kernel
---------------------

Implemented by the developers of the :term:`driver`.

* Manage global storage, including reading spatially and/or temporally
  varying material properties, geochemical conditions, etc.
* loop through space
* manage time stepping, substepping, etc
* unpack and move data from the mesh dependent storage into alquimia
  data transfer containers
* I/O

Alquimia Interface
------------------

* Defines the engine independent API, including function call
  signatures and data transfer :doc:`containers <APIv1_structures>`.

* The API must be compatible with mixed language programing, C++
  calling fortran or fortran calling C, etc.

* Handles the engine dependent details for setup, shutdown, processing
  constraints, reaction stepping, etc.

* At the start of on operation, it unpacks data from the alquimia data
  transfer containers and packages them into the correct format for
  the engine. At the end of the operation, in packages the results
  back into the alquimia containers.

* **Alquimia does NOT do any geochemical calculations.**

* Implementation details : Each engine interface will handle the
  alquimia operations in it's native implementation language.


Engine
------

* The :term:`engine` is responsibly for all geochemistry calculations, including:
    * managing the reaction network (database reading, basis swapping).
    * constraint processing.
    * reaction stepping (OS).
    * returning reaction step jacobian and rhs evaluations (GI).
    * provideding auxiliary output, e.g. pH, mineral saturation index, reaction rates, etc

* The maintainers of each engine are responsible for providing a
  wrapper library that **EXACTLY** conforms to the alquimia API.

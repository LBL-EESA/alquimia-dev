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


Alquimia Constants
==================

These are the named constants that alquimia uses and their required values. 

These should be implemented as C constants or Fortran parameters, not
preprocessor macros. Because they are constant, they can be
global. Implementations should be described in header files
``alquimia_constants.h`` or module files ``alquimia_constants.F90``
and compiled into ``libalquimia_c.a`` or ``libalquimia_fortran.a``
respectively.


Alquimia Error Codes
~~~~~~~~~~~~~~~~~~~~

Every alquimia function call should be followed by checking the
AlquimiaEngineStatus.error against the following codes:

+---------------------------------------+-----------+-----------------------------------+
| **name**                              | **value** | **meaning**                       |
+---------------------------------------+-----------+-----------------------------------+
| kAlquimiaNoError                      | 0         |no error                           |
+---------------------------------------+-----------+-----------------------------------+
| kAlquimiaErrorInvalidEngine           | 1         |unknown engine requested           |
+---------------------------------------+-----------+-----------------------------------+
| kAlquimiaErrorUnknownConstraintName   | 2         |the user requested a named         |
|                                       |           |constraint which could not be found|
|                                       |           |in the constraint list             |
+---------------------------------------+-----------+-----------------------------------+
| kAlquimiaErrorUnsupportedFunctionality| 3         |the user requested functionality   |
|                                       |           |that is not supported with the     |
|                                       |           |current engine                     |
+---------------------------------------+-----------+-----------------------------------+
| kAlquimiaErrorEngineIntegrity         | 4577      |pointer to the engine's internal   |
|                                       |           |state did not pass integrity check |
+---------------------------------------+-----------+-----------------------------------+

Alquimia String Lengths
~~~~~~~~~~~~~~~~~~~~~~~

+--------------------------+---------------+
| **name**                 | **value**     |
+--------------------------+---------------+
| kAlquimiaMaxStringLength | 512           |
+--------------------------+---------------+
| kAlquimiaMaxWordLength   | 32            |
+--------------------------+---------------+

.. _AlquimiaStrings:

Alquimia Strings
~~~~~~~~~~~~~~~~

String comparisons should be case insensitive. The only exception is
chemical symbols, which should follow standard conventions, i.e. "H+",
not "h+", "CO" = carbon monoxide, "Co" = cobalt, "co" = undefined.

+--------------------------------------+----------------------------+
| **name**                             | **value**                  |
+--------------------------------------+----------------------------+
| kAlquimiaStringPFloTran              | "PFloTran"                 |
+--------------------------------------+----------------------------+
|  kAlquimiaStringCrunchFlow           |  "CrunchFlow"              |
+--------------------------------------+----------------------------+
| kAlquimiaStringTotalAqueous          | "total_aqueous"            |
+--------------------------------------+----------------------------+
| kAlquimiaStringTotalSorbed           | "total_sorb"               |
+--------------------------------------+----------------------------+
| kAlquimiaStringFree                  | "free"                     |
+--------------------------------------+----------------------------+
| kAlquimiaStringpH                    | "pH"                       |
+--------------------------------------+----------------------------+
| kAlquimiaStringMineral               | "mineral"                  |
+--------------------------------------+----------------------------+
| kAlquimiaStringGas                   | "gas"                      |
+--------------------------------------+----------------------------+
| kAlquimiaStringCharge                | "charge"                   |
+--------------------------------------+----------------------------+



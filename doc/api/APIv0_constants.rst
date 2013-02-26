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
| kAlquimiaStringTotalSorbed           | "total_sorbed"             |
+--------------------------------------+----------------------------+
| kAlquimiaStringTotalAqueousPlusSorbed| "total_aqueous_plus_sorbed"|
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



Alquimia Constants
==================

These are the named constants that alquimia uses and their required values. 

Alquimia Error Codes
~~~~~~~~~~~~~~~~~~~~

+-------------------------------+-----------+-----------------------------------+
| **name**                      | **value** | **meaning**                       |
+-------------------------------+-----------+-----------------------------------+
| kAlquimiaNoError              | 0         |no error                           |
+-------------------------------+-----------+-----------------------------------+
| kAlquimiaErrorInvalidEngine   | 1         |unknown engine requested           |
+-------------------------------+-----------+-----------------------------------+
| kAlquimiaErrorEngineIntegrity | 4577      |pointer to the engine's internal   |
|                               |           |state did not pass integrity check |
+-------------------------------+-----------+-----------------------------------+

Alquimia String Lengths
~~~~~~~~~~~~~~~~~~~~~~~

+----------+---------------+---------------+
| **name**                 | **value**     |
+----------+---------------+---------------+
| kAlquimiaMaxStringLength | 512           |
+----------+---------------+---------------+
| kAlquimiaMaxWordLength   | 32            |
+----------+---------------+---------------+


Alquimia Strings
~~~~~~~~~~~~~~~~

+---------------------------+--------------+
| **name**                  | **value**    |
+---------------------------+--------------+
| kAlquimiaStringPFloTran   | "PFloTran"   |
+---------------------------+--------------+
| kAlquimiaStringCrunchFlow | "CrunchFlow" |
+---------------------------+--------------+
| kAlquimiaStringpH         | "pH"         |
+---------------------------+--------------+
| kAlquimiaStringMineral    | "mineral"    |
+---------------------------+--------------+
| kAlquimiaStringGas        | "gas"        |
+---------------------------+--------------+
| kAlquimiaStringCharge     | "charge"     |
+---------------------------+--------------+

NOTE(bja): Except for species names, strings should be case insensitive? The demo driver (and C++ chem) use lower case, but pflotran uses a mixture of StringToUpper and StringToLower internally.


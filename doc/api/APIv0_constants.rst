Alquimia Constants
==================

These are the named constants that alquimia uses and their required values. 

Alquimia Error Codes
~~~~~~~~~~~~~~~~~~~~

Every alquimia function call should be followed by checking the
AlquimiaEngineStatus.error against the following codes:

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

String comparisons should be case insensitive. The only exception is
chemical sysmols, which should follow standard conventions, i.e. "H+",
not "h+", "CO" = carbon monoxide, "Co" = cobalt, "co" = undefined.

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

NOTE(bja): 


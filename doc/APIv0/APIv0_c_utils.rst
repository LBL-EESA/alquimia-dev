Alquimia C Utilities Library
============================

Details of the optional alquimia C utility library. This is not a
required part of the alquimia API.

Alquimia Interface
~~~~~~~~~~~~~~~~~~

C struct with function pointers to the alquimia interface functions
(``alquimia_intefrace.h`` and ``alquimia_interface.c``) and engine
state for the desired engine.

+---------------------------+------------------+
| **variable**              | **type**         |
+---------------------------+------------------+
| Setup                     | function pointer |
+---------------------------+------------------+
| Shutdown                  | function pointer |
+---------------------------+------------------+
| ProcessCondition          | function pointer |
+---------------------------+------------------+
| ReactionStepOperatorSplit | function pointer |
+---------------------------+------------------+
| GetAuxiliaryOutput        | function pointer |
+---------------------------+------------------+
| GetEngineMetaData         | function pointer |
+---------------------------+------------------+
| engine_state              | void*            |
+---------------------------+------------------+

Creating an Alquimia Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function CreateAlquimiaInterface(), ``alquimia_intefrace.h``, will
create an alquimia interface struct with the correct function pointers for a
given engine.

.. code-block:: c

    void CreateAlquimiaInterface(char* engine_name,
        struct AlquimiaInterface* interface,
        struct AlquimiaEngineStatus* status)


Alquimia Data
~~~~~~~~~~~~~

Convenience data structure holding all the log lived Alquimia data in
a single struct, ``alquimia_intefrace.h``. Eases storage of multiple
copies of data (e.g. when threading with OpenMP), but not part of the
interface.

+----------------------------+---------------------+
| **variable**               | **type**            |
+----------------------------+---------------------+
| AlquimiaSizes              | sizes               |
+----------------------------+---------------------+
| AlquimiaState              | state               |
+----------------------------+---------------------+
| AlquimiaMaterialProperties | material_properties |
+----------------------------+---------------------+
| AlquimiaAuxiliaryData      | aux_data            |
+----------------------------+---------------------+
| AlquimiaMetaData           | meta_data           |
+----------------------------+---------------------+

Alquimia: Memory
~~~~~~~~~~~~~~~~

C library for allocating and freeing memory in the alquimia structs,
``alquimia_memory.h`` and ``alquimia_memory.c``.

Alquimia Print Utils
~~~~~~~~~~~~~~~~~~~~

C library of common utilities for working with the contents of the
alquimia structs, ``alquimia_utils.h`` and
``alquimia_utils.c``. Primarily functions for printing the contents of
alquimia structs.

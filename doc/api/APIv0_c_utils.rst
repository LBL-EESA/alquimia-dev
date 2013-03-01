Alquimia C Utilities Library
============================

Alquimia provides an optional C utility library,
``libalquimia_cutils.a``. This is not a required part of the alquimia
API, but a collection of reusable code for common operations when
working with alquimia from C or C++.

Alquimia Interface
------------------

Code to ease working the alquimia interface, ``alquimia_intefrace.h``
and ``alquimia_interface.c``.

Struct: Alquimia Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~

``AlquimiaInterface {...}`` is a C structure containing function pointers
 to the alquimia interface functions for the desired engine.
 
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
| GetProblemMetaData        | function pointer |
+---------------------------+------------------+

Creating an Alquimia Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function ``CreateAlquimiaInterface(...)``, will create an alquimia
interface structure with the correct function pointers for a given
engine.  All ``#ifdef HAVE_XXX`` code for a particular geochemistry
engine is contained in this function.

.. code-block:: c

    void CreateAlquimiaInterface(const char* const engine_name,
        struct AlquimiaInterface* interface,
        struct AlquimiaEngineStatus* status)


Struct: Alquimia Data
~~~~~~~~~~~~~~~~~~~~~

``AlquimiaData {...}`` is a convenience data structure holding all the long
lived Alquimia data in a single structure. Eases storage of multiple
copies of data (e.g. when threading with OpenMP), but not part of the
formal alquimia API.

.. code-block:: c

    struct AlquimiaData {
      void* engine_state;
      struct AlquimiaSizes sizes;
      struct AlquimiaEngineFunctionality functionality;
      struct AlquimiaState state;
      struct AlquimiaMaterialProperties material_properties;
      struct AlquimiaAuxiliaryData aux_data;
      struct AlquimiaProblemMetaData meta_data;
      struct AlquimiaAuxiliaryOutputData aux_output;
    };


Alquimia Memory
----------------

``alquimia_memory.h`` and ``alquimia_memory.c`` contain a C library
for allocating and freeing memory in the alquimia structures.

If an AlquimiaData struct is created and the 'sizes' member is
correctly filled out, then a call to

.. code-block:: c

   void AllocateAlquimiaData(struct AlquimiaData* data)

will allocate all the internal memory necessary.

Individual containers can be allocated by calling the appropriate
allocate function:

.. code-block:: c

    void AllocateAlquimiaXXX(const struct AlquimiaSizes* const sizes,
                             struct AlquimiaXXX* xxx)

where XXX is the name of an alquimia structure. Every allocate
function has a corresponding free function

.. code-block:: c

    void FreeAlquimiaXXX(struct AlquimiaXXX* xxx)

Alquimia Utils
--------------

``alquimia_utils.h`` and ``alquimia_utils.c`` contain common utilities
for working with the contents of the alquimia data.

Printing
~~~~~~~~

Calling ``PrintAlquimiaXXX`` will pretty-print the contents of alquimia data structure XXX to the screen.

.. code-block:: c

    void PrintAlquimiaXXX(const struct AlquimiaXXX* const xxx)

Strings
~~~~~~~

Compare two alquimia strings. True if they are equivalent, false otherwise.

.. code-block:: c

  bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                            const char* const str2);



Species Name-Index Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~

Determine the **C** index for a particular species name. Sets ``index = -1``
if the name is not in the provided list.

.. code-block:: c

  void AlquimiaFindIndexFromName(const char* const name,
                                 const struct AlquimiaVectorString* const names,
                                 int* index);


If the **engine** index is needed, then the driver needs to use the index_base offset from the AlquimiaEngineFunctionality structure, i.e.

.. code-block:: c

  engine_index = c_index + functionality.index_base;


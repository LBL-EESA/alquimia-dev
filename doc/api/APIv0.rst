API
===

Alquimia will use `Semantic Versioning <http://semver.org/>`_ for its public API.


* :doc:`function <APIv0_functions>` call signatures
* :doc:`data structures <APIv0_structures>`
* :doc:`constants <APIv0_constants>`
* :doc:`c utility library <APIv0_c_utils>`

Required Functionality
~~~~~~~~~~~~~~~~~~~~~~

Alquimia exposes the following generic functionality from the underlying libraries:

* Reading of geochemical reaction data

* Basis management (includes reading database, swapping basis, etc.)

* Constraint processing for boundary/initial constraints. (Say you set a boundary concentration to be equilibrated with a mineral. The total component concentration for the species equilibrated with the mineral needs to be calculated and sent back to Amanzi transport.)

* Reaction step(ping), operator split and global implicit

* Access to user selected geochemical data for output, i.e. pH, mineral SI, reaction rates

* Additional information about the library, e.g. is it thread safe, so amanzi can create multiple copies using openmp? Does it support temperature dependent chemistry, porosity updates, etc.

Notes
~~~~~

* "Driver" : the transport simulator or other driver, e.g. amanzi

* "Engine" : the backend geochemistry engine, e.g. pflotran, crunchflow, toughreact, stomp, phreeqc, ...

* Alquimia has two parts.
    * An engine independent API consisting of :doc:`function <APIv0_functions>` call signatures and :doc:`data structures <APIv0_structures>`.
    * An optional :doc:`utility <APIv0_c_utils>` library to handle data memory allocation/freeing, printing structs, etc.

* Implementation details:
    * Each client will have some sort of process kernel interface to handle the native memory --> alquimia memory.
    * Each engine interface will handle the alquimia --> engine memory mapping and engine specific function calls. Most of this side of the interface will take place in the language of the engine.


Division of labor / responsibilities for implementing alquimia
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PK
--

Amanzi-U, Amanzi-S, ParFlow etc. Implemented by the developers of the PK.

* Manage global storage, including reading spatially and/or temporally varying material properties, geochemical conditions, etc.
* loop through space
* manage time stepping, substepping, etc
* unpack and move data from the mesh dependent storage into alquimia data transfer containers
* I/O

Alquimia Interface
------------------

defines the engine independent API, including function call signatures and data transfer [[APIv0_Structures | containers]].

Handles the engine dependent details for setup, shutdown, processing constraints, reaction stepping, etc. 

**Alquimia does NOT do any geochemical calculations**

Engine
------

PFloTran, CrunchFlow, phreeqc, tough react, etc.

Responsibly for all geochemistry calculations, including database reading, basis swapping, setting up reaction networks, constraint processing, reaction stepping (OS), returning reaction step jacobian and rhs evaluations (GI).

The maintainers of each engine are responsible for providing a wrapper library that conforms to the alquimia API.

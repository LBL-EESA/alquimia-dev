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
   
   Authors: Benjamin Andre <bandre@lbl.gov>, Sergi Molins <smolins@lbl.gov, 
   Jeffrey Johnson <jnjohnson@lbl.gov>


Guidance for use of Alquimia
============================

Alquimia aims at providing a flexible interface to access a range of 
capabilities available in geochemical engines. There is no single way to use
Alquimia. The developer that wishes to use the API to add geochemical capabilities
to a driver can decide at what level the driver will control aspects of the 
geochemical problem. In the simplest approach, the transport driver will hand 
over complete control over the geochemical problem to the geochemical engine 
via Alquimia; in particular, in assigning initial concentrations as well as 
parameter values to the different geochemical processes. However, this may be 
not be desired in all cases, e.g. when the driver wants to perform sensitivity 
analyses or inverse modeling on geochemical input parameters. In this scenario,
Alquimia provides structures and functions that afford lower level control over 
parameter values used by the geochemical engine. In both scenarios, however, the
stochiometric relationships defining reactive processes and model formulations 
are always specified by the geochemical engine. In fact, these can be different 
for different engines, which would make it difficult to develop a generic 
interface. 

In most applications, the concept of geochemical condition is central in providing
the connection between the driver and the engine. A geochemical condition is a 
consistent set of concentrations of aqueous, sorbed and mineral species 
based on a number of constraints, e.g. mineral equilibrium, gas equilibrium, 
solution charge balance, pH value, etc. Speciation or *processing* of this
set of constraint fully defined the geochemical state in a porous medium cell,  
both in the aqueous and the solid phases.

In the following, guidance is given as to how to use Alquimia in a number of 
possible modes of operation, ranging from hands-off to hands-on control of 
geochemical parameters by the driver.

Hands-off mode
---------------------

NOTE: This mode is in development.

In this mode of operation, complete control of the geochemical problem
is handed over to Alquimia and thus to the geochemical engine. This mode of 
operation is the simplest to code. Because the geochemical engine
has full control of the geocehmical problem, a user who is already familiar with
the structure of the input file of a particular engine can benefit 
from this mode of operation. Thus, the hands-off coupling
is recommended to developers interested in ease of implementation and
current users of particular geochemical engine that are nevertheless interested
in accessing capabilities not available said geochemical code (e.g. a specific
discretization method for transport, a mesh type, etc.).

The connection between the input of the driver's input and the geochemical 
engine input file is done via named geochemical conditions. Specifically, 
the names of geochemical conditions are given in the driver's input. These 
names are then passed to Alquimia. The condition constraints are provided in the 
geochemical engine's input file. Once, Alquimia hands the name of the condition
to the engine, this processes the condition and assigns the 
concentration values and material property values to the cell of interest.

The suggested call sequence in this mode of operation is given in the following 
pseudocode: 

.. code-block:: c

   //alquimia initialization (alquimia_interface.h)     
   AllocateAlquimiaEngineStatus(AlquimiaStatus);
   CreateAlquimiaInterface(Alquimia, 
              AlquimiaStatus);
   Alquimia.Setup(AlquimiaData.sizes);

   //alquimia data allocation (alquimia_memory.h)
   AllocateAlquimiaData(AlquimiaData);

   //driver initialization (driver-specific calls)
   readDriverInput(domainSize,
              conditionNames);
   allocateDriverState(domainSize,
              AlquimiaData.sizes, 
              DriverState);

   //processing conditions - assigning initial values to cells & boundaries
   for(int cell=0;cell<domainSize;cell++) {
     // send in the name of the condition to apply to this cell, 
     // get out the state (i.e. set of concentrations)
     Alquimia.ProcessCondition(conditionNames[cell],
                                                  AlquimiaData.state);
     // copy to driver state
     CopyAlquimiaToDriver(AlquimiaData.state,
                DriverState[cell]);
   }
   //time stepping (operator splitting)
   for(int itime=0;itime<endTime;itime++) {
     // driver (e.g. flow, transport) advance
     DriverTimeStep(DriverState);
     // loop over domain (Alquimia is cell-by-cell)
     for(int cell=0;cell<domainSize;cell++) { 
          // copy to alquimia state
          CopyDriverToAlquimia(DriverState[cell],
                      AlquimiaData.state);
          // Alquimia chemistry advance      
          Alquimia.ReactionStepOperatorSplit(AlquimiaData.state);
          // copy to driver state
          CopyAlquimiaToDriver(AlquimiaData.state,
                      DriverState[cell]);
     }
   }

*Note:* Some :doc:`Alquimia functions <api/APIv1_functions>` have been used 
here but the signature is not that of the API. The signature used here is only
for the purpose of illustrating the hands-off coupling, with an emphasis on the 
arguments that are central to this approach. The exact signature of the 
functions :doc:`as defined by API <api/APIv1_functions>` is required for actual
coding.
   
Hands-on mode
-------------
   
In this mode of operation, control of the geochemical problem is shared between
the driver and the geochemical engine. This mode of 
operation is harder to code as it requires a lower-level interaction between
the driver and the engine. For example, the driver needs to check whether the
size and order of the geochemical parameters in the driver are compatible with
those in the geochemical engine. Further, because the geochemical engine does 
not have full control of the geocehmical problem, this approach may be prone to 
error if the division of responsibilites is not made clear to the user. Thus, 
the hands-on coupling is recommended to developers familiar with geochemical 
modeling, to users who may not be familiar with a particular geochemical 
engine's input format, and for applications that require changing the values of 
input parameters for geochemistry without interacting with the geochemical 
engine's input file. In any case, however, the stochiometric relationships 
defining reactive processes and model formulations (e.g. reaction rate 
expressions) are always specified by the geochemical engine. 

The connection between the input of the driver's input and the geochemical 
engine input file is again done via named geochemical conditions. However, in 
this case, this not done by passing the names of the names of geochemical 
conditions but actually assembling the geochemical conditions constraints. 
These conditions are then passed to Alquimia (a data structure is provided for 
that). Alquimia takes care of translating the condition constraint to the 
geochemical engine's internal format, which processes the condition as it were 
provided in its own input file. The concentration values and material property 
values are assigned to the cell of interest.

The suggested call sequence in this mode of operation is given in the following 
pseudocode: 

.. code-block:: c

   //alquimia initialization (alquimia_interface.h)
   AllocateAlquimiaEngineStatus(AlquimiaStatus);
   CreateAlquimiaInterface(Alquimia, 
              AlquimiaStatus);
   Alquimia.Setup(AlquimiaData.sizes);

   //alquimia data allocation (alquimia_memory.h)
   AllocateAlquimiaData(AlquimiaData);

   // obtain alquimia problem meta data
   Alquimia.GetProblemMetaData(AlquimiaData.meta_data);

   //driver initialization (driver-specific calls)
   readInput(domainSize,
             chemistrySizes,
             chemistryProblemMetaData,
             initialConcentrations,
             initialChemicalProperties,
             initialConstraints);
   allocateDriverState(domainSize,
              AlquimiaData.sizes, 
              DriverState,
              DriverChemicalProperties);

   // check compatibility of sizes, order and names of parameters
   if (AlquimiaData.sizes != chemistrySizes) {
     error;
   }
   if (AlquimiaData.meta_data != chemistryProblemMetaData) {
     createMap(AlquimiaData.meta_data,
               chemistryProblemMetaData);
   }
   // assemblage of conditions w/ constraints
   for(int cell=0;cell<domainSize;cell++) {
     AlquimiaGeochemicalCondition[cell] = 
       AssembleGeochemicalCondition(initialConcentrations[cell],
                                    initialMaterialProperties[cell],
                                    initialConstraints[cell]);
   }                                                                                          

   //processing conditions - assigning initial values to cells & boundaries
   for(int cell=0;cell<domainSize;cell++) {
     // copy to Drive ChemicalProperties to Alquimia Properties (can be cell-by-cell)
     CopyDriverToAlquimia(DriverChemicalProperties,
                              AlquimiaData.properties);
     // send in condition assembled by Driver to apply to this cell, 
     // get out the state (i.e. set of concentrations)    
     AlquimiaProcessCondition(AlquimiaGeochemicalCondition[cell],
                              AlquimiaData.state);
     // copy to driver state
     CopyAlquimiaToDriver(AlquimiaData.state,
                              DriverState[cell]);
   }
   //time stepping (operator splitting)
   for(int itime=0;itime<endTime;itime++) {
     // driver (e.g. flow, transport) advance
     DriverTimeStep(DriverState)
     // loop over domain (Alquimia is cell-by-cell)
     for(int cell=0;cell<domainSize;cell++) { 
          // copy to alquimia state
          CopyDriverToAlquimia(DriverState[cell],
                      AlquimiaData.state);
          // Alquimia chemistry advance      
          Alquimia.ReactionStepOperatorSplit(AlquimiaData.state);
          // copy to driver state
          CopyAlquimiaToDriver(AlquimiaData.state,
                      DriverState[cell]);
     }
   }

*Note:* Some :doc:`Alquimia functions <api/APIv1_functions>` and 
:doc:`Alquimia structures <api/APIv1_structures>` have been used 
here but the signature is not that of the API. The signature used here is only
for the purpose of illustrating the hands-on coupling, with an emphasis on the 
arguments that are central to this approach. The exact signature of the 
functions :doc:`as defined by the API <api/APIv1_functions>` is required for 
actual coding. The same applies to 
:doc:`Alquimia structures <api/APIv1_structures>`.

   
Fine-grained control of chemistry feedback processes    
----------------------------------------------------

The discussion in this section is only intended to give a high level view of the 
interaction between the driver and Alquimia, especially regarding problem setup
and initialization. Naturally, the interaction between the driver and Alquimia 
can go beyond what is described here. In general, feedback processes between 
flow or transport and geochemical processes can be considered. As part of 
Alquimia.Setup, the AlquimiaData.functionality data structure is returned.
This contains information about the functionality supported by the engine. 
It is up to the driver to decide (at run time if desired) whether this 
functionality is to be used. For example, the geochemical driver may update 
porosity based on changes in mineral volume changes. This updated porosity is 
returned after the call to Alquimia.ReactionStepOperatorSplit. The driver should 
use AlquimiaData.functionality to be aware whether this is the case in any 
particular simulation and decide to discard those changes or to use them. For
some engines, this functionality can depend on the options provided in its
input file, thus it is considered good practice check (at run time) what 
behavior to expect from an engine in any particular simulation.

As it is apparent from this discussion, the driver can alter the data that comes
out of and goes into Alquimia at will. This provides the lowest level of control
of the driver on the geochemical problem. However, this needs to be used with
care. For example, arbitrarily altering concentrations may result in convergence
issues next time the nonlinear geochemical problem is solved.

Alquimia state and properties
---------------------------------------

In the previous section, porosity is given as an example of a variable that 
can change in an Alquimia time step. In Alquimia, a distinction is made 
between variables defined in the AlquimiaData.state and AlquimiaData.properties
structures. AlquimiaData.state contains variables that can change in time. For
example, aqueous concentrations are part of the AlquimiaData.state but also mineral 
volume fractions. In contrast, AlquimiaData.properties contains variables 
that are constant over time (within Alquimia) but that may be different in 
different parts of the domain considered by the driver. This include aqueous
saturation but also Kd coefficients.

Often in the development for flow and transport, a different definition is given
for State variables and material properties. For example, material properties 
can be properties that are associated the solid phase (thus, immobile), while 
state variables are those that are subject to transport (thus, mobile). In this 
definition, concentrations are state variables, but mineral volume fractions are
material properties. The developer must be careful in handling this correctly.

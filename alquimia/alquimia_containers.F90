
!
! Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
! through Lawrence Berkeley National Laboratory (subject to receipt of any 
! required approvals from the U.S. Dept. of Energy).  All rights reserved.
! 
! Alquimia is available under a BSD license. See LICENSE.txt for more
! information.
!
! If you have questions about your rights to use or distribute this software, 
! please contact Berkeley Lab's Technology Transfer and Intellectual Property 
! Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
! 
! NOTICE.  This software was developed under funding from the U.S. Department 
! of Energy.  As such, the U.S. Government has been granted for itself and 
! others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
! license in the Software to reproduce, prepare derivative works, and perform 
! publicly and display publicly.  Beginning five (5) years after the date 
! permission to assert copyright is obtained from the U.S. Department of Energy, 
! and subject to any subsequent five (5) year renewals, the U.S. Government is 
! granted for itself and others acting on its behalf a paid-up, nonexclusive, 
! irrevocable, worldwide license in the Software to reproduce, prepare derivative
! works, distribute copies to the public, perform publicly and display publicly, 
! and to permit others to do so.
! 
! Authors: Benjamin Andre <bandre@lbl.gov>
!

! **************************************************************************** !
!
! Alquimia Containers module
!
! Author: Benjamin Andre
!
! WARNINGS:
!
!   * The alquimia data structures defined in the this are dictated by
!     the alquimia API! Do NOT change them unless you make
!     corresponding changes to the c containers (and doc).
!
!   * The order of the data members matters! If num_primary is the
!     first member and num_minerals is the third, then they must be in
!     those positions on both the c and fortran side of the interface!
!
!   * The names of the containers and their members should be the same
!     for both c and fortran. The language interface doesn't require
!     this, but it makes it easier for the reader to understand what is
!     going on.
!
! **************************************************************************** !

module AlquimiaContainers_module

  use, intrinsic :: iso_c_binding

  implicit none

  integer (c_int), parameter :: kAlquimiaMaxStringLength = 512
  integer (c_int), parameter :: kAlquimiaMaxWordLength = 32

  integer (c_int), parameter :: kAlquimiaNoError = 0
  integer (c_int), parameter :: kAlquimiaErrorInvalidEngine = 1
  integer (c_int), parameter :: kAlquimiaErrorUnknownConstraintName = 2
  integer (c_int), parameter :: kAlquimiaErrorUnsupportedFunctionality = 3
  integer (c_int), parameter :: kAlquimiaErrorEngineIntegrity = 4577

  character (13), parameter :: kAlquimiaStringTotalAqueous = 'total_aqueous'
  character (12), parameter :: kAlquimiaStringTotalSorbed = 'total_sorbed'
  character (25), parameter :: kAlquimiaStringTotalAqueousPlusSorbed = 'total_aqueous_plus_sorbed'
  character (4), parameter :: kAlquimiaStringFree = 'free'
  character (2), parameter :: kAlquimiaStringPH = 'pH'
  character (7), parameter :: kAlquimiaStringMineral = 'mineral'
  character (3), parameter :: kAlquimiaStringGas = 'gas'
  character (6), parameter :: kAlquimiaStringCharge = 'charge'

  type, public, bind(c) :: AlquimiaVectorDouble
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type AlquimiaVectorDouble

  type, public, bind(c) :: AlquimiaVectorInt
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type AlquimiaVectorInt

  type, public, bind(c) :: AlquimiaVectorString
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type AlquimiaVectorString

  type, public, bind(c) :: AlquimiaSizes
     integer (c_int) :: num_primary
     integer (c_int) :: num_sorbed
     integer (c_int) :: num_minerals
     integer (c_int) :: num_aqueous_complexes
     integer (c_int) :: num_aqueous_kinetics
     integer (c_int) :: num_surface_sites
     integer (c_int) :: num_ion_exchange_sites
     integer (c_int) :: num_isotherm_species
     integer (c_int) :: num_aux_integers
     integer (c_int) :: num_aux_doubles
  end type AlquimiaSizes

  type, public, bind(c) :: AlquimiaState
     real (c_double) :: water_density
     real (c_double) :: porosity
     real (c_double) :: temperature
     real (c_double) :: aqueous_pressure
     type (AlquimiaVectorDouble) :: total_mobile
     type (AlquimiaVectorDouble) :: total_immobile
     type (AlquimiaVectorDouble) :: mineral_volume_fraction
     type (AlquimiaVectorDouble) :: mineral_specific_surface_area
     type (AlquimiaVectorDouble) :: surface_site_density
     type (AlquimiaVectorDouble) :: cation_exchange_capacity
  end type AlquimiaState

  type, public, bind(c) :: AlquimiaProperties
     real (c_double) :: volume
     real (c_double) :: saturation
     type (AlquimiaVectorDouble) :: isotherm_kd
     type (AlquimiaVectorDouble) :: freundlich_n
     type (AlquimiaVectorDouble) :: langmuir_b
     type (AlquimiaVectorDouble) :: mineral_rate_cnst
!!     real (c_double) :: solid_density
     type (AlquimiaVectorDouble) :: aqueous_kinetic_rate_cnst
  end type AlquimiaProperties

  type, public, bind(c) :: AlquimiaAuxiliaryData 
     type (AlquimiaVectorInt) :: aux_ints
     type (AlquimiaVectorDouble) :: aux_doubles
  end type AlquimiaAuxiliaryData

  type, public, bind(c) :: AlquimiaEngineStatus
     integer (c_int) :: error
     type (c_ptr) :: message
     logical (c_bool) :: converged
     integer (c_int) :: num_rhs_evaluations
     integer (c_int) :: num_jacobian_evaluations
     integer (c_int) :: num_newton_iterations
  end type AlquimiaEngineStatus

  type, public, bind(c) :: AlquimiaEngineFunctionality
     logical (c_bool) :: thread_safe
     logical (c_bool) :: temperature_dependent
     logical (c_bool) :: pressure_dependent
     logical (c_bool) :: porosity_update
     logical (c_bool) :: operator_splitting
     logical (c_bool) :: global_implicit
     integer (c_int) :: index_base
  end type AlquimiaEngineFunctionality

  type, public, bind(c) :: AlquimiaProblemMetaData
     type (AlquimiaVectorString) :: primary_names
     type (AlquimiaVectorInt)    :: positivity
     type (AlquimiaVectorString) :: mineral_names
     type (AlquimiaVectorString) :: surface_site_names
     type (AlquimiaVectorString) :: ion_exchange_names
     type (AlquimiaVectorString) :: isotherm_species_names
     type (AlquimiaVectorString) :: aqueous_kinetic_names
  end type AlquimiaProblemMetaData

  type, public, bind(c) :: AlquimiaAuxiliaryOutputData
     real (c_double) :: pH
     type (AlquimiaVectorDouble) :: aqueous_kinetic_rate
     type (AlquimiaVectorDouble) :: mineral_saturation_index
     type (AlquimiaVectorDouble) :: mineral_reaction_rate
     type (AlquimiaVectorDouble) :: primary_free_ion_concentration
     type (AlquimiaVectorDouble) :: primary_activity_coeff
     type (AlquimiaVectorDouble) :: secondary_free_ion_concentration
     type (AlquimiaVectorDouble) :: secondary_activity_coeff
  end type AlquimiaAuxiliaryOutputData

  type, public, bind(c) :: AlquimiaAqueousConstraint
     type (c_ptr) :: primary_species_name
     type (c_ptr) :: constraint_type
     type (c_ptr) :: associated_species
     real (c_double) :: value
  end type AlquimiaAqueousConstraint

  type, public, bind(c) :: AlquimiaAqueousConstraintVector
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type AlquimiaAqueousConstraintVector

  type, public, bind(c) :: AlquimiaMineralConstraint
     type (c_ptr) :: mineral_name
     real (c_double) :: volume_fraction
     real (c_double) :: specific_surface_area
  end type AlquimiaMineralConstraint

  type, public, bind(c) :: AlquimiaMineralConstraintVector
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type AlquimiaMineralConstraintVector

  type, public, bind(c) :: AlquimiaGeochemicalCondition
     type (c_ptr) :: name
     type (AlquimiaAqueousConstraintVector) :: aqueous_constraints
     type (AlquimiaMineralConstraintVector) :: mineral_constraints
  end type AlquimiaGeochemicalCondition

end module AlquimiaContainers_module

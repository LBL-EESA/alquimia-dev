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
  integer (c_int), parameter :: kAlquimiaErrorEngineIntegrity = 4577

  type, public, bind(c) :: AlquimiaVectorDouble
     integer (c_int) :: size
     type (c_ptr) :: data
  end type AlquimiaVectorDouble

  type, public, bind(c) :: AlquimiaVectorInt
     integer (c_int) :: size
     type (c_ptr) :: data
  end type AlquimiaVectorInt

  type, public, bind(c) :: AlquimiaVectorString
     integer (c_int) :: size
     type (c_ptr) :: data
  end type AlquimiaVectorString

  type, public, bind(c) :: AlquimiaSizes
     integer (c_int) :: num_primary
     integer (c_int) :: num_sorbed
     integer (c_int) :: num_kinetic_minerals
     integer (c_int) :: num_aqueous_complexes
     integer (c_int) :: num_surface_sites
     integer (c_int) :: num_ion_exchange_sites
  end type AlquimiaSizes

  type, public, bind(c) :: AlquimiaState
     real (c_double) :: water_density
     real (c_double) :: saturation
     real (c_double) :: porosity
     real (c_double) :: temperature
     real (c_double) :: aqueous_pressure
     type (AlquimiaVectorDouble) :: total_primary
     type (AlquimiaVectorDouble) :: total_sorbed
     type (AlquimiaVectorDouble) :: free_ion
     type (AlquimiaVectorDouble) :: mineral_volume_fraction
     type (AlquimiaVectorDouble) :: mineral_specific_surface_area
     type (AlquimiaVectorDouble) :: cation_exchange_capacity
     type (AlquimiaVectorDouble) :: surface_site_density
  end type AlquimiaState

  type, public, bind(c) :: AlquimiaMaterialProperties
     real (c_double) :: volume
     type (AlquimiaVectorDouble) :: isotherm_kd
     type (AlquimiaVectorDouble) :: freundlich_n
     type (AlquimiaVectorDouble) :: langmuir_b
  end type AlquimiaMaterialProperties

  type, public, bind(c) :: AlquimiaAuxiliaryData 
     type (AlquimiaVectorDouble) :: primary_activity_coeff
     type (AlquimiaVectorDouble) :: secondary_activity_coeff
     type (AlquimiaVectorDouble) :: ion_exchange_ref_cation_conc
     type (AlquimiaVectorDouble) :: surface_complex_free_site_conc
  end type AlquimiaAuxiliaryData

  type, public, bind(c) :: AlquimiaEngineStatus
     integer (c_int) :: error
     type (c_ptr) :: message
     logical (c_bool) :: converged
     integer (c_int) :: num_rhs_evaluations
     integer (c_int) :: num_jacobian_evaluations
     integer (c_int) :: num_newton_iterations
  end type AlquimiaEngineStatus

  type, public, bind(c) :: AlquimiaMetaData
     logical (c_bool) :: thread_safe
     logical (c_bool) :: temperature_dependent
     logical (c_bool) :: pressure_dependent
     logical (c_bool) :: porosity_update
     logical (c_bool) :: operator_splitting
     logical (c_bool) :: global_implicit
     integer (c_int) :: index_base
     type (AlquimiaVectorInt) :: primary_indices
     type (AlquimiaVectorString) :: primary_names
     type (AlquimiaVectorInt) :: mineral_indices
     type (AlquimiaVectorString) :: mineral_names
  end type AlquimiaMetaData

  type, public, bind(c) :: AlquimiaAuxiliaryOutputData
     real (c_double) pH
     type (AlquimiaVectorDouble) mineral_saturation_index
     type (AlquimiaVectorDouble) mineral_reaction_rate
  end type AlquimiaAuxiliaryOutputData

  type, public, bind(c) :: AlquimiaConstraint
     type (c_ptr) :: primary_species
     type (c_ptr) :: constraint_type
     type (c_ptr) :: associated_species
     real (c_double) :: value
  end type AlquimiaConstraint

  type, public, bind(c) :: AlquimiaCondition
     type (c_ptr) :: name
     integer (c_int) :: num_constraints
     type (c_ptr) :: constraints
  end type AlquimiaCondition

end module AlquimiaContainers_module

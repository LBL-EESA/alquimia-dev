
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
! PFLOTRAN Alquimia Inteface module
!
! Author: Benjamin Andre
!
! Notes:
!
!  * Public function call signatures, including intent, are dictated
!    by the alquimia API.
!
!  * alquimia data structures defined in the AlquimiaContainers_module
!    (alquimia_containers.F90) are dictated by the alquimia API.
!
!  * All function calls involving pflotran native data structures must
!    be private!
!
!  * (bja) 2012-12 - different fortran compilers use different name
!    mangling conventions for fortran modules:
!
!    gfortran : ___modulename_MOD_procedurename
!
!    intel : _modulename_mp_procedurename_
!
!    as a consequence we can't call the alquimia interface directly
!    from C, but need to use the simple wrappers in
!    pflotran_alquimia_wrappers.F90!
!
! **************************************************************************** !

module PFLOTRANAlquimiaInterface_module

  ! pflotran modules
  use Option_module, only : option_type
  use Reaction_Aux_module, only : reaction_rt_type
  use Reaction_Base_module, only : reaction_base_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Material_Aux_class, only : material_auxvar_type
  use Transport_Constraint_module, only : tran_constraint_list_type
  use Transport_Constraint_RT_module, only : tran_constraint_coupler_rt_type

  use PFLOTRAN_Constants_module
#include "finclude/petscsys.h"
#include "petsc/finclude/petscsys.h"
  implicit none

  public :: Setup, &
       Shutdown, &
       ProcessCondition, &
       ReactionStepOperatorSplit, &
       GetAuxiliaryOutput, &
       GetProblemMetaData

  private :: &
       SetupPFLOTRANOptions, &
       SetEngineFunctionality, &
       SetAlquimiaSizes, &
       InitializeTemperatureDependence, &
       InitializePFLOTRANReactions, &
       ReadPFLOTRANConstraints, &
       ProcessPFLOTRANConstraint, &
       ConvertAlquimiaConditionToPflotran, &
       CopyAlquimiaToAuxVars, &
       CopyAuxVarsToAlquimia, &
       GetAuxiliaryDataSizes, &
       PackAlquimiaAuxiliaryData, &
       UnpackAlquimiaAuxiliaryData, &
       PrintTranConstraint, &
       PrintAqueousSpeciesConstraint

  integer(kind=8), private, parameter :: integrity_check_value = 5784429996817932654_8
  !integer(kind=int64), private, parameter :: integrity_check_value = &
  !     b"0101000001000110011011000110111101010100011100100110000101101110"

  type, private :: PFLOTRANEngineState
     ! This is the data structure that stores the persistent data for
     ! pflotran, (e.g. reaction network).
     !
     ! We pass it back and forth with c as a void pointer so we don't
     ! have to a global variable.
     !
     ! It is NOT part of the alquimia interface, and the driver code
     ! should not use or depend on it!
     !
     ! NOTE(bja): these are fortran pointers, so this struct can not
     ! be unpacked on the c side!
     integer(kind=8) :: integrity_check
     logical :: hands_off
     type(option_type), pointer :: option
     class(reaction_rt_type), pointer :: reaction
     type(reactive_transport_auxvar_type), pointer :: rt_auxvar
     type(global_auxvar_type), pointer :: global_auxvar
     class(material_auxvar_type), pointer :: material_auxvar
     type(tran_constraint_list_type), pointer :: transport_constraints
     class(tran_constraint_coupler_rt_type), pointer :: constraint_coupler 
  end type PFLOTRANEngineState

contains


! **************************************************************************** !
!
! Public interface routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine Setup(input_filename, hands_off, pft_engine_state, sizes, &
                 functionality, status)
!  NOTE: Function signature is dictated by the alquimia API.
!
!  NOTE: Assumes that MPI_Init() and / or PetscInitialize() have already
!    been called by the driver

  use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_bool

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_Aux_module, only : reaction_rt_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type, &
       RTAuxVarInit
  use Global_Aux_module, only : global_auxvar_type, GlobalAuxVarInit
  use Material_Aux_class, only : material_auxvar_type, MaterialAuxVarInit
  use Option_module, only : option_type, OptionCreate
  use Input_Aux_module, only : input_type, InputCreate, InputDestroy
  use Transport_Constraint_module, only : tran_constraint_list_type
  use Transport_Constraint_RT_module, only : tran_constraint_coupler_rt_type, &
                                             TranConstraintCouplerRTCreate, &
                                             TranConstraintRTCreate

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  logical (c_bool), value, intent(in) :: hands_off
  type (c_ptr), intent(out) :: pft_engine_state
  type (AlquimiaSizes), intent(out) :: sizes
  type (AlquimiaEngineFunctionality), intent(out) :: functionality
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFLOTRANEngineState), pointer :: engine_state
  PetscBool :: option_found
  character(len=kAlquimiaMaxStringLength) :: string
  class(reaction_rt_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(global_auxvar_type), pointer :: global_auxvar
  class(material_auxvar_type), pointer :: material_auxvar
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(tran_constraint_list_type), pointer :: transport_constraints
  class(tran_constraint_coupler_rt_type), pointer :: constraint_coupler 

  integer :: len
  integer :: i

  write (*, '(a)') "PFLOTRAN_Alquimia_Setup() : "

  ! setup pflotran's option object, including mpi
  option => OptionCreate()
  call SetupPFLOTRANOptions(input_filename, option)

  write (*, '(a, a)') "  Reading : ", trim(option%input_filename)
  input => InputCreate(IN_UNIT, option%input_filename, option)

  call InitializeScreenOutput(option, input)

  !
  ! initialize chemistry
  !
  call InitializeTemperatureDependence(option, input)
  call InitializePFLOTRANReactions(option, input, reaction)

  !
  ! create the various storage containers pflotran will need later
  !
  ! NOTE(bja) : batch chem --> one cell

  ! global_auxvar --> cell by cell temperature, pressure, saturation, density
  allocate(global_auxvar)
  call GlobalAuxVarInit(global_auxvar, option)

  ! material_auxvar --> cell by cell volume, porosity, etc.
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar, option)

  ! assign default state values, not really needed?
  global_auxvar%pres = option%flow%reference_pressure
  global_auxvar%temp = option%flow%reference_temperature
  global_auxvar%den_kg = option%flow%reference_density(option%liquid_phase)
  global_auxvar%sat = option%flow%reference_saturation  

  material_auxvar%porosity = option%flow%reference_porosity

  ! rt_auxvar --> cell by cell chemistry data
  allocate(rt_auxvar)
  call RTAuxVarInit(rt_auxvar, reaction, option)

  ! NOTE(bja) : constraint coupler is only needed for printing the
  ! processed constraints to pflotran.out. But destroying the coupler
  ! destroys the auxvar as well, so we need to keep it around long
  ! term.
  constraint_coupler => TranConstraintCouplerRTCreate(option)

  !
  ! Read the constraints so we can finish using the input file.
  !
  call ReadPFLOTRANConstraints(option, input, reaction, transport_constraints)

  ! close the input file because we don't need it any more
  call InputDestroy(input)

  call SetAlquimiaSizes(reaction, sizes)
  !FIXME(bja): need some way to identify pflotran ion exchange sites
  !when there are more than one. Use the mineral name?
  if (sizes%num_ion_exchange_sites > 1) then
     status%error = kAlquimiaErrorUnsupportedFunctionality
     call f_c_string_ptr("ERROR: pflotran interface only supports a single ion exchange site!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  call SetEngineFunctionality(reaction, option, functionality)

  !
  ! save pflotran's persistent data to a struct so the driver can
  ! store it for us
  !
  allocate(engine_state)
  engine_state%integrity_check = integrity_check_value
  engine_state%hands_off = hands_off
  engine_state%option => option
  engine_state%reaction => reaction
  engine_state%rt_auxvar => rt_auxvar
  engine_state%global_auxvar => global_auxvar
  engine_state%material_auxvar => material_auxvar
  engine_state%constraint_coupler => constraint_coupler
  engine_state%transport_constraints => transport_constraints

  pft_engine_state = c_loc(engine_state)

  status%error = kAlquimiaNoError
  call f_c_string_ptr("Alquimia::PFLOTRAN::Setup() : successful.", &
       status%message, kAlquimiaMaxStringLength)
  
end subroutine Setup


! **************************************************************************** !
subroutine Shutdown(pft_engine_state, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! pflotran
  use Option_module, only : OptionDestroy
  use Driver_module, only : DriverDestroy
  use Reaction_aux_module, only : ReactionDestroy
  use Transport_Constraint_Base_module, only : tran_constraint_coupler_base_type
  use Transport_Constraint_module, only : TranConstraintCouplerDestroy, &
                                          TranConstraintListDestroy

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFLOTRANEngineState), pointer :: engine_state
  class(tran_constraint_coupler_base_type), pointer :: constraint_coupler_base

  !write (*, '(a)') "PFLOTRANAlquimiaInterface::Shutdown() : "

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !TODO(geh): remove when TranConstraintAddToList() has been refactored
  !           with target instead of pointer
  constraint_coupler_base => engine_state%constraint_coupler
  call TranConstraintCouplerDestroy(constraint_coupler_base)
  call TranConstraintListDestroy(engine_state%transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(engine_state%rt_auxvar)
  !call GlobalAuxVarDestroy(engine_state%global_auxvar)
  call ReactionDestroy(engine_state%reaction, engine_state%option)
  call DriverDestroy(engine_state%option%driver)
  call OptionDestroy(engine_state%option)

  deallocate(engine_state)
  nullify(engine_state)

  status%error = kAlquimiaNoError

end subroutine Shutdown


! **************************************************************************** !
subroutine ProcessCondition(pft_engine_state, condition, properties, &
     state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  use c_f_interface_module, only : c_f_string_ptr, f_c_string_ptr

  use AlquimiaContainers_module

  ! pflotran
  use String_module, only : StringCompareIgnoreCase
  use Transport_Constraint_module, only : TranConstraintAddToList
  use Transport_Constraint_RT_module, only : tran_constraint_rt_type
  use Transport_Constraint_Base_module, only : tran_constraint_base_type
  use Option_module, only : PrintMsg

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaGeochemicalCondition), intent(in) :: condition
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent (inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  real (c_double), pointer :: data(:)
  character (kAlquimiaMaxStringLength) :: name
  class(tran_constraint_rt_type), pointer :: tran_constraint
  class(tran_constraint_base_type), pointer :: tran_constraint_base
  real (c_double) :: constraint_value
  PetscReal :: porosity, volume
  integer :: i
  type(PFLOTRANEngineState), pointer :: engine_state
  PetscInt, parameter :: phase_index = 1
  logical, parameter :: copy_auxdata = .false.

  !write (*, '(a)') "PFLOTRANAlquimiaInterface::ProcessCondition() : "

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  ! NOTE(bja): the data stored in alquimia's aux_data is uninitialized
  ! at this point, so don't want to copy it! (copy_auxdata = false)
  call CopyAlquimiaToAuxVars(copy_auxdata, engine_state%hands_off, &
       state, aux_data, properties, &
       engine_state%reaction, engine_state%global_auxvar, &
       engine_state%material_auxvar, engine_state%rt_auxvar)

  !
  ! process the condition
  !
  call  c_f_string_ptr(condition%name, name)
  engine_state%option%io_buffer = "processing : " // trim(name)
  call PrintMsg(engine_state%option)

  if (condition%aqueous_constraints%size > 0) then
     ! the driver is supplying the constraint data, so we need to
     ! construct an object in pflotran's internal format.
     tran_constraint => &
       ConvertAlquimiaConditionToPflotran(engine_state%option, &
                                          engine_state%reaction, condition)
     ! add to the list to ensure memory gets cleaned up.
     !TODO(geh): remove when TranConstraintAddToList() has been refactored
     !           with target instead of pointer
     tran_constraint_base => tran_constraint
     call TranConstraintAddToList(tran_constraint_base, &
                                  engine_state%transport_constraints)
     engine_state%constraint_coupler%constraint => tran_constraint_base
  else
    ! the driver just supplied a name, so we check for a constraint
    ! with that name in the pflotran input file and use that.
    engine_state%option%io_buffer = "Looking for pflotran constraint : " // &
                                     trim(name)
     call PrintMsg(engine_state%option)
     tran_constraint_base => engine_state%transport_constraints%first
     do
        if (associated(tran_constraint_base)) then
           ! check the name of this constraint
           if (StringCompareIgnoreCase(tran_constraint_base%name, name)) then
              ! found the constraint we are looking for, bail from the
              ! loop with the current pointer
              engine_state%constraint_coupler%constraint => tran_constraint_base
              select type(c => tran_constraint_base)
                class is(tran_constraint_rt_type)
                  tran_constraint => c
              end select
              exit
           else
              ! check the next constraint
              tran_constraint_base => tran_constraint_base%next
           end if
        else
           ! end of the list (or empty list) without out finding a match.
           call f_c_string_ptr("INPUT ERROR: Could not find pflotran constraint : "//trim(name), &
                status%message, kAlquimiaMaxStringLength)
           status%error = kAlquimiaErrorUnknownConstraintName
           return
        end if
     end do
  end if

  if (associated(tran_constraint)) then
     !call PrintTranConstraint(tran_constraint)
     ! tran_constraint should be valid. Now we can ask pflotran to
     ! process it...
     call ProcessPFLOTRANConstraint( &
          engine_state%option, &
          engine_state%reaction, &
          engine_state%global_auxvar, &
          engine_state%material_auxvar, &
          engine_state%rt_auxvar, &
          engine_state%constraint_coupler)
     ! now repack the processed constraint data into the alquimia
     ! struct for the driver.
     call CopyAuxVarsToAlquimia( &
          engine_state%reaction, &
          engine_state%global_auxvar, &
          engine_state%rt_auxvar, &
          porosity, &
          state, aux_data)
     status%error = kAlquimiaNoError
  end if
end subroutine ProcessCondition


! **************************************************************************** !
subroutine ReactionStepOperatorSplit(pft_engine_state, &
     delta_t, properties, state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

#include "petsc/finclude/petscsys.h"
  use petscsys
  use, intrinsic :: iso_c_binding, only : c_ptr, c_double, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_module, only : RReact, RUpdateKineticState, RTAuxVarCompute

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  real (c_double), value, intent(in) :: delta_t
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFLOTRANEngineState), pointer :: engine_state
  PetscReal :: porosity, volume, vol_frac_prim
  PetscReal, allocatable :: guess(:)
  PetscInt :: i, num_newton_iterations, ierror
  PetscInt, parameter :: natural_id = -999
  PetscInt, parameter :: phase_index = 1
  logical, parameter :: copy_auxdata = .true.
  class(reaction_rt_type), pointer :: reaction

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  reaction => engine_state%reaction
  !write (*, '(a)') "F_PFLOTRANAlquimiaInterface::ReactionStepOperatorSplit() :"

  !call PrintState(state)
  
  call CopyAlquimiaToAuxVars(copy_auxdata, engine_state%hands_off, &
       state, aux_data, properties, &
       reaction, engine_state%global_auxvar, &
       engine_state%material_auxvar, engine_state%rt_auxvar)

  ! copy free ion primaries into initial guess array
  allocate(guess(reaction%ncomp))
  do i = 1, reaction%naqcomp
     guess(i) = engine_state%rt_auxvar%pri_molal(i)
  enddo
  do i = 1, reaction%immobile%nimmobile
     guess(i+reaction%offset_immobile) = engine_state%rt_auxvar%immobile(i)
  enddo

  vol_frac_prim = 1.0

  engine_state%option%tran_dt = delta_t

  !TODO(bja): if (.not.option%use_isothermal) then call RUpdateTempDependentCoefs...?

!!$  call RTAuxVarCompute(engine_state%rt_auxvar, &
!!$                       engine_state%global_auxvar, &
!!$                       engine_state%reaction, engine_state%option)

  call RReact(guess, engine_state%rt_auxvar, engine_state%global_auxvar, &
       engine_state%material_auxvar, num_newton_iterations, &
       reaction, natural_id, engine_state%option, &
       PETSC_FALSE, PETSC_FALSE, ierror)
  deallocate(guess)


  if (ierror /= 1) then

     call RUpdateKineticState(engine_state%rt_auxvar, engine_state%global_auxvar, &
          engine_state%material_auxvar, engine_state%reaction, engine_state%option)

     call CopyAuxVarsToAlquimia( &
          engine_state%reaction, &
          engine_state%global_auxvar, &
          engine_state%rt_auxvar, &
          porosity, &
          state, aux_data)

  ! Copy the diagnostic information into the status object.
  ! PFlotran doesn't do anything really fancy in its Newton step, so
  ! the numbers of RHS evaluations, Jacobian evaluations, and Newton 
  ! iterations are all the same.
     status%error = kAlquimiaNoError
     status%converged = .true.
     status%num_rhs_evaluations = num_newton_iterations
     status%num_jacobian_evaluations = num_newton_iterations
     status%num_newton_iterations = num_newton_iterations

  else
  ! not converged (ierror == 1)

     status%error = kAlquimiaNoError
     status%converged = .false.
     status%num_rhs_evaluations = num_newton_iterations
     status%num_jacobian_evaluations = num_newton_iterations
     status%num_newton_iterations = num_newton_iterations

  endif

end subroutine ReactionStepOperatorSplit



! **************************************************************************** !
subroutine GetAuxiliaryOutput( &
     pft_engine_state, &
     properties, &
     state, &
     aux_data, &
     aux_output, &
     status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_double, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_Mineral_module, only : RMineralSaturationIndex

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaAuxiliaryOutputData), intent(inout) :: aux_output
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFLOTRANEngineState), pointer :: engine_state
  integer :: i, ph_index
  real (c_double), pointer :: local_array(:)
  PetscReal :: porosity, volume
  logical, parameter :: copy_auxdata = .true.

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !write (*, '(a)') "PFLOTRANAlquimiaInterface::GetAuxiliaryOutput() :"

  ! FIXME(bja): right now, all info from the previous reaction step is
  ! still in the auxvars. We are assuming that the driver has called
  ! GetAuxiliaryOutput immediately after reaction step...!

  !call CopyAlquimiaToAuxVars(copy_auxdata, state, aux_data, properties, &
  !     engine_state%reaction, engine_state%global_auxvar, engine_state%rt_auxvar, &
  !     porosity, volume)

  ph_index = engine_state%reaction%species_idx%h_ion_id
  ! FIXME(bja): this violates the no geochemistry calculations in alquimia rule.
  if (ph_index > 0) then
     aux_output%pH = -log10(engine_state%rt_auxvar%pri_act_coef(ph_index) * &
          engine_state%rt_auxvar%pri_molal(ph_index))
  else
     ! FIXME(bja, 2013-07) need a meaningful N/A value
     aux_output%pH = -100.d0
  end if

  ! Aqueous kinetic rate FIXME: not yet supported!
!  call c_f_pointer(aux_output%aqueous_kinetic_rate%data, local_array, &
!       (/aux_output%aqueous_kinetic_rate%size/))
!  do i = 1, aux_output%aqueous_kinetic_rate%size
!    local_array(i) = engine_state%global_auxvar%reaction_rate(i)
!  end do

  ! mineral data
  call c_f_pointer(aux_output%mineral_reaction_rate%data, local_array, &
       (/aux_output%mineral_reaction_rate%size/))
  do i = 1, aux_output%mineral_reaction_rate%size
    local_array(i) = engine_state%rt_auxvar%mnrl_rate(i)
  end do

  call c_f_pointer(aux_output%mineral_saturation_index%data, local_array, &
       (/aux_output%mineral_saturation_index%size/))
  do i = 1, aux_output%mineral_saturation_index%size
     local_array(i) = RMineralSaturationIndex(i, engine_state%rt_auxvar, &
          engine_state%global_auxvar, &
          engine_state%reaction, engine_state%option)
  end do

  !
  ! primary data
  !
  call c_f_pointer(aux_output%primary_free_ion_concentration%data, local_array, &
       (/aux_output%primary_free_ion_concentration%size/))
  do i = 1, engine_state%reaction%naqcomp
     local_array(i) = engine_state%rt_auxvar%pri_molal(i)
  end do
  ! Immobile species have zero free ion conc by definition
  do i = 1, engine_state%reaction%immobile%nimmobile
    local_array(i+engine_state%reaction%offset_immobile) = 0.0
 end do

  call c_f_pointer(aux_output%primary_activity_coeff%data, local_array, &
       (/aux_output%primary_activity_coeff%size/))
  do i = 1, aux_output%primary_activity_coeff%size
     local_array(i) = engine_state%rt_auxvar%pri_act_coef(i)
  end do
  ! Immobile species have zero primary activity
  do i = 1, engine_state%reaction%immobile%nimmobile
    local_array(i+engine_state%reaction%offset_immobile) = 0.0
   end do

  !
  ! secondary aqueous complex data
  !
  call c_f_pointer(aux_output%secondary_free_ion_concentration%data, local_array, &
       (/aux_output%secondary_free_ion_concentration%size/))
  do i = 1, aux_output%secondary_free_ion_concentration%size
     local_array(i) = engine_state%rt_auxvar%sec_molal(i)
  end do


  call c_f_pointer(aux_output%secondary_activity_coeff%data, local_array, &
       (/aux_output%secondary_activity_coeff%size/))
  do i = 1, aux_output%secondary_activity_coeff%size
     local_array(i) = engine_state%rt_auxvar%sec_act_coef(i)
  end do

  status%error = kAlquimiaNoError
end subroutine GetAuxiliaryOutput


! **************************************************************************** !
subroutine GetProblemMetaData(pft_engine_state, meta_data, status)
  !  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_int, c_char, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr, f_c_string_chars

  use AlquimiaContainers_module

  use Reaction_Aux_module, only : general_rxn_type
  use Reaction_Microbial_Aux_module, only : microbial_rxn_type

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaProblemMetaData), intent(out) :: meta_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type (c_ptr), pointer :: name_list(:)
  character (len=kAlquimiaMaxWordLength), pointer :: pflotran_names(:)
  character (len=kAlquimiaMaxWordLength) :: dummy_names(1)
  character (c_char), pointer :: name
  integer :: i, list_size, id
  integer(c_int), pointer :: idata(:)
  type(PFLOTRANEngineState), pointer :: engine_state
  type(general_rxn_type), pointer :: cur_gen_rxn
  type(microbial_rxn_type), pointer :: cur_mic_rxn

  !write (*, '(a)') "PFLOTRAN_Alquimia_GetEngineMetaData() :"

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !
  ! copy primary indices and names
  !

! In Pflotran, reaction%ncomp includes primary species, colloids, and immobile species
! so if a simulation includes colloids or immobile species, this may need to change to reaction%naqcomp
  if (meta_data%primary_names%size /= engine_state%reaction%ncomp) then
     write (*, '(a, i3, a, i3, a)') "meta_data%primary_names%size (", &
          meta_data%primary_names%size, ") != pflotran%reaction%ncomp(", &
          engine_state%reaction%ncomp, ")"
  end if
  list_size = engine_state%reaction%naqcomp

  pflotran_names => engine_state%reaction%primary_species_names

  call c_f_pointer(meta_data%primary_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(pflotran_names(i)), &
          name, kAlquimiaMaxStringLength)     
  end do

  ! Immobile species names
  list_size = engine_state%reaction%ncomp
  
  pflotran_names => engine_state%reaction%immobile%names
  
  call c_f_pointer(meta_data%primary_names%data, name_list, (/list_size/))
  do i = 1, engine_state%reaction%immobile%nimmobile
     call c_f_pointer(name_list(i+engine_state%reaction%offset_immobile), name)
     call f_c_string_chars(trim(pflotran_names(i)), &
          name, kAlquimiaMaxStringLength)     
  end do


!
! positivity constraints
!

  if (meta_data%positivity%size /= engine_state%reaction%ncomp) then
     write (*, '(a, i3, a, i3, a)') "meta_data%positivity%size (", &
          meta_data%positivity%size, ") != pflotran%reaction%ncomp(", &
          engine_state%reaction%ncomp, ")"
      stop
  end if
  list_size = meta_data%positivity%size

  call c_f_pointer(meta_data%positivity%data, idata, (/list_size/))
  do i = 1, list_size
      if (i == engine_state%reaction%species_idx%h_ion_id) then
!       H+ component can be negative
        idata(i) = 0
      else
        idata(i) = 1 
      end if
  end do

  !
  ! copy mineral indices and names
  !

  if (meta_data%mineral_names%size /= engine_state%reaction%mineral%nkinmnrl) then
     write (*, '(a, i3, a, i3, a)') "meta_data%mineral_names%size (", &
          meta_data%mineral_names%size, ") != pflotran%reaction%mineral%nkinmnrl(", &
          engine_state%reaction%mineral%nkinmnrl, ")"
  end if
  list_size = meta_data%mineral_names%size

  pflotran_names => engine_state%reaction%mineral%mineral_names

  call c_f_pointer(meta_data%mineral_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(pflotran_names(i)), &
          name, kAlquimiaMaxStringLength)     
  end do

  !
  ! surface sites
  !
  if (meta_data%surface_site_names%size /= &
       engine_state%reaction%surface_complexation%nsrfcplxrxn) then
     write (*, '(a, i3, a, i3, a)') "meta_data%surface_site_names%size (", &
          meta_data%surface_site_names%size, ") != pflotran%reaction%surface_complexation%nsrfcplxrxn(", &
          engine_state%reaction%surface_complexation%nsrfcplxrxn, ")"
  end if

  list_size = meta_data%surface_site_names%size

  pflotran_names => engine_state%reaction%surface_complexation%srfcplxrxn_site_names

  call c_f_pointer(meta_data%surface_site_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(pflotran_names(i)), &
          name, kAlquimiaMaxStringLength)     
  end do


  !
  ! ion exchange
  !
  ! FIXME(bja): only support a single ion exchange site
  if (meta_data%ion_exchange_names%size /= &
       engine_state%reaction%neqionxrxn) then
     write (*, '(a, i3, a, i3, a)') "meta_data%ion_exchange_names%size (", &
          meta_data%ion_exchange_names%size, ") != pflotran%reaction%neqionxrxn(", &
          engine_state%reaction%neqionxrxn, ")"
  end if

  list_size = meta_data%ion_exchange_names%size

  dummy_names(1) = "X-"

  call c_f_pointer(meta_data%ion_exchange_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(dummy_names(i)), &
          name, kAlquimiaMaxStringLength)     
  end do


  !
  ! isotherm indices
  !
  list_size = meta_data%isotherm_species_names%size
  call c_f_pointer(meta_data%isotherm_species_names%data, name_list, &
       (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     id = engine_state%reaction%isotherm%eqkdspecid(i)
     call f_c_string_chars( &
          trim(engine_state%reaction%primary_species_names(id)), &
          name, kAlquimiaMaxStringLength)
  end do
  
  !
  ! aqueous kinetic names
  !
  list_size = meta_data%aqueous_kinetic_names%size
  call c_f_pointer(meta_data%aqueous_kinetic_names%data, name_list, &
       (/list_size/))
   ! General reactions
   cur_gen_rxn => engine_state%reaction%general_rxn_list
   i=1
   do
     if (.not.associated(cur_gen_rxn)) exit
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars( &
          trim(cur_gen_rxn%reaction), &
          name, kAlquimiaMaxStringLength)
    cur_gen_rxn => cur_gen_rxn%next
    i = i+1
   enddo
   ! Microbial reactions
   cur_mic_rxn => engine_state%reaction%microbial%microbial_rxn_list
   do
     if (.not.associated(cur_mic_rxn)) exit
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars( &
          trim(cur_mic_rxn%reaction), &
          name, kAlquimiaMaxStringLength)
    cur_mic_rxn => cur_mic_rxn%next
    i = i+1
   enddo

  status%error = 0
end subroutine GetProblemMetaData


! **************************************************************************** !
!
! Private work routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine SetupPFLOTRANOptions(input_filename, option)

  use, intrinsic :: iso_c_binding, only : c_char

  use c_f_interface_module, only : c_f_string_chars

  use Option_module, only : option_type
  use Driver_module, only : DriverCreate
  use Communicator_Aux_module, only : CommCreate
  use PFLOTRAN_Constants_module
  use petscsys
  implicit none



  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (option_type), intent(inout) :: option

  ! local variables
  character(len=MAXSTRINGLENGTH) :: filename_out
  PetscErrorCode :: ierr

  ! set the pflotran input file name
  call c_f_string_chars(input_filename, option%input_filename)

  ! setup the driver and comm
  option%driver => DriverCreate()
  option%driver%comm => CommCreate()
  option%comm => option%driver%comm

  option%global_prefix = option%input_filename

  !
  ! mpi
  !
  option%comm%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, option%comm%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, option%comm%global_commsize, ierr)
  call MPI_Comm_group(MPI_COMM_WORLD, option%comm%global_group, ierr)
  option%comm%mycomm = option%comm%global_comm
  option%comm%myrank = option%comm%global_rank
  option%comm%mycommsize = option%comm%global_commsize
  option%comm%mygroup = option%comm%global_group

  !
  ! output file
  !
  option%fid_out = DRIVER_OUT_UNIT

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out.alquimia'

  if (option%myrank == option%driver%io_rank .and. &
      option%driver%print_to_file) then
    open(option%driver%fid_out, file=filename_out, action="write", &
         status="unknown")
  endif

  !
  ! manually initialize some required values
  !
  option%nphase = 1
  option%liquid_phase = 1
  option%flow%reference_density(option%liquid_phase) = 997.16
  option%use_isothermal = PETSC_TRUE

end subroutine SetupPFLOTRANOptions

! **************************************************************************** !
subroutine SetEngineFunctionality(reaction, option, functionality)

  use AlquimiaContainers_module, only : AlquimiaEngineFunctionality

  use Reaction_aux_module, only : reaction_rt_type
  use Option_module, only : option_type

  implicit none

  ! function parameters
  class (reaction_rt_type), intent(in) :: reaction
  type (option_type), intent(in) :: option
  type (AlquimiaEngineFunctionality), intent(out) :: functionality

  functionality%thread_safe = .true.
  functionality%temperature_dependent = .not. option%use_isothermal
  functionality%pressure_dependent = .false.
  functionality%porosity_update = reaction%update_porosity
  functionality%operator_splitting = .true.
  functionality%global_implicit = .false.
  functionality%index_base = 1

end subroutine SetEngineFunctionality

! **************************************************************************** !
subroutine SetAlquimiaSizes(reaction, sizes)

  use AlquimiaContainers_module, only : AlquimiaSizes

  use Reaction_aux_module, only : reaction_rt_type

  implicit none

  ! function parameters
  class (reaction_rt_type), intent(in) :: reaction
  type (AlquimiaSizes), intent(out) :: sizes

! This may need to change to reaction%naqcomp if there are colloids and/or immobile species defined
  sizes%num_primary = reaction%ncomp
  if (reaction%nsorb > 0 .or. reaction%immobile%nimmobile > 0) then
     sizes%num_sorbed = reaction%ncomp
  else
     sizes%num_sorbed = 0
  end if
  sizes%num_minerals = reaction%mineral%nkinmnrl
  sizes%num_aqueous_complexes = reaction%neqcplx
  sizes%num_aqueous_kinetics = reaction%ngeneral_rxn + reaction%microbial%nrxn
  sizes%num_surface_sites = reaction%surface_complexation%nsrfcplxrxn
  sizes%num_ion_exchange_sites = reaction%neqionxrxn
  sizes%num_isotherm_species = reaction%isotherm%neqkdrxn
  
  call GetAuxiliaryDataSizes(reaction, &
       sizes%num_aux_integers, sizes%num_aux_doubles)

  !call PrintSizes(sizes)

end subroutine SetAlquimiaSizes

! **************************************************************************** !
subroutine InitializeScreenOutput(option, input)

  ! pflotran
  use Option_module, only : option_type, PrintMsg, PrintErrMsg
  use Input_Aux_module, only : input_type, InputReadPflotranString, InputReadWord, &
       InputCheckExit, InputError, InputErrorMsg, InputReadStringErrorMsg
  use String_module, only : StringToUpper
  use petscsys
  implicit none



  ! function parameters
  type(option_type), pointer, intent(inout) :: option
  type(input_type), pointer, intent(in) :: input

  ! local variables
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word, card

  ! default value: option%print_to_screen = PETSC_TRUE

  rewind(input%fid)

  do
    call InputReadPflotranString(input, option)
    if (InputError(input)) exit

    call InputReadWord(input, option, word, PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    select case(trim(card))
!....................
      case ('OUTPUT')
        do
          call InputReadPflotranString(input, option)
          call InputReadStringErrorMsg(input, option, card)
          if (InputCheckExit(input, option)) exit
          call InputReadWord(input, option, word, PETSC_TRUE)
          call InputErrorMsg(input, option, 'keyword', 'OUTPUT')
          call StringToUpper(word)
          select case(trim(word))
            case('SCREEN')
              call InputReadWord(input, option, word, PETSC_TRUE)
              call InputErrorMsg(input, option, 'time increment', 'OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%driver%print_to_screen = PETSC_FALSE
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,SCREEN for alquimia-pflotran interface.'
                  call PrintErrMsg(option)
               end select
            end select
         end do
      end select
   end do

end subroutine InitializeScreenOutput

! **************************************************************************** !
subroutine InitializeTemperatureDependence(option, input)

  ! pflotran
  use Option_module, only : option_type
  use Input_Aux_module, only : input_type, InputFindStringInFile, InputError
  use petscsys

  implicit none

  ! function parameters
  type(option_type), pointer, intent(inout) :: option
  type(input_type), pointer, intent(in) :: input

  ! local variables
  character(len=MAXSTRINGLENGTH) :: string

  ! NOTE(bja): we are assuming isothermal by default. then override
  ! from the input file if nonisothermal was specified

  option%use_isothermal = PETSC_TRUE

  string = 'NONISOTHERMAL'
  call InputFindStringInFile(input, option, string)
  if (.not. InputError(input)) then
     ! found a nonisothermal block
        option%use_isothermal = PETSC_FALSE
  end if

end subroutine InitializeTemperatureDependence

! **************************************************************************** !
subroutine InitializePFLOTRANReactions(option, input, reaction)

  ! pflotran
  use Reaction_module, only : ReactionInit, ReactionReadPass2
  use Reaction_Aux_module, only : reaction_rt_type, ACT_COEF_FREQUENCY_OFF
  use Reaction_Database_module, only : DatabaseRead, BasisInit
  use Option_module, only : option_type
  use Input_Aux_module, only : input_type, InputFindStringInFile, InputError
  use petscsys
  implicit none

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  type(input_type), pointer, intent(in) :: input
  class(reaction_rt_type), pointer, intent(out) :: reaction

  ! local variables
  character(len=MAXSTRINGLENGTH) :: string

  ! check for a chemistry block in the  input file
  string = "CHEMISTRY"
  call InputFindStringInFile(input, option, string)
  if (.not. InputError(input)) then
    ! found a chemistry block, initialize the chemistry.

    ! NOTE(bja): ReactionInit() only does a first pass through the
    ! input file to check for a select items
    call ReactionInit(reaction, input, option)
    ! rewind the input file to prepare for the second pass
    call InputFindStringInFile(input, option, string)
    ! the second pass through the input file to read the remaining blocks
    call ReactionReadPass2(reaction, input, option)
  else
     ! TODO(bja): no chemistry block --> fatal error
  endif
    
  if (associated(reaction)) then
    if (reaction%use_full_geochemistry) then
       call DatabaseRead(reaction, option)
       call BasisInit(reaction, option)    
    else
      ! NOTE(bja): do we need this for the batch chemistry driver?

      ! turn off activity coefficients since the database has not been read
      reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(reaction%primary_species_print(option%ntrandof))
      reaction%primary_species_print = PETSC_TRUE
    endif
  endif

end subroutine InitializePFLOTRANReactions


! **************************************************************************** !
subroutine ReadPFLOTRANConstraints(option, input, reaction, transport_constraints)
!  NOTE: We are just reading the data from the input file. No processing
!    is done here because we don't know yet if we are using these.

  use Reaction_Aux_module, only : reaction_rt_type
  use Option_module, only : option_type, PrintMsg, PrintErrMsg
  use Input_Aux_module, only : input_type, InputReadPflotranString, InputReadWord, &
       InputErrorMsg, InputError
  use String_module, only : StringToUpper
  use Transport_Constraint_module, only : tran_constraint_list_type, &
        TranConstraintInitList, TranConstraintAddToList
  use Transport_Constraint_RT_module, only : tran_constraint_rt_type, &
        TranConstraintRTRead, TranConstraintRTCreate
  use Transport_Constraint_Base_module, only : tran_constraint_base_type
  use petscsys
  implicit none

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  type(input_type), pointer, intent(inout) :: input
  class(reaction_rt_type), pointer, intent(inout) :: reaction
  type(tran_constraint_list_type), pointer, intent(inout) :: transport_constraints

  ! local variables
  PetscBool :: debug = PETSC_FALSE
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  class(tran_constraint_rt_type), pointer :: tran_constraint
  class(tran_constraint_base_type), pointer :: tran_constraint_base

  allocate(transport_constraints)
  call TranConstraintInitList(transport_constraints)

  ! look through the input file
  rewind(input%fid)        
  do
    call InputReadPflotranString(input, option)
    if (InputError(input)) exit

    call InputReadWord(input, option, word, PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    if (debug) then
       option%io_buffer = 'pflotran card:: ' // trim(card)
       call PrintMsg(option)
    end if

    select case(trim(card))
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call PrintErrMsg(option)
        endif
        tran_constraint => TranConstraintRTCreate(option)
        !TODO(geh): remove when TranConstraintAddToList() has been refactored
        !           with target instead of pointer
        tran_constraint_base => tran_constraint
        call InputReadWord(input, option, tran_constraint%name, PETSC_TRUE)
        call InputErrorMsg(input, option, 'constraint', 'name') 
        if (debug) then
           call PrintMsg(option, tran_constraint%name)
        end if
        call TranConstraintRTRead(tran_constraint, reaction, input, option)
        call TranConstraintAddToList(tran_constraint_base, &
                                     transport_constraints)
        nullify(tran_constraint)

      case default
         ! do nothing
    end select
  enddo

end subroutine ReadPFLOTRANConstraints

! **************************************************************************** !
subroutine ProcessPFLOTRANConstraint(option, reaction, global_auxvar, &
                                     material_auxvar, rt_auxvar, &
                                     constraint_coupler)

  use Reaction_module, only : ReactionProcessConstraint, &
        ReactionPrintConstraint, ReactionEquilibrateConstraint
  use Reaction_Aux_module, only : reaction_rt_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Material_Aux_class, only : material_auxvar_type
  use Transport_Constraint_RT_module, only : tran_constraint_rt_type, &
                                             tran_constraint_coupler_rt_type
  use Option_module, only : option_type, PrintMsg, PrintErrMsg
  use petscsys
  implicit none

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  class(reaction_rt_type), pointer, intent(inout) :: reaction
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvar
  class(material_auxvar_type), pointer, intent(inout) :: material_auxvar
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvar
  class(tran_constraint_coupler_rt_type), pointer, intent(inout) :: constraint_coupler 

  ! local variables
  class(tran_constraint_rt_type), pointer :: tran_constraint
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: use_prev_soln_as_guess
  PetscInt :: num_iterations
  integer :: i
  !
  ! process constraints
  !
  num_iterations = 0
  use_prev_soln_as_guess = PETSC_FALSE

  select type(c => constraint_coupler%constraint)
    class is(tran_constraint_rt_type)
      tran_constraint => c
  end select

  if (.not. associated(tran_constraint)) then
     ! TODO(bja) : report error
     option%io_buffer = "tran_constraint not associated"
     call PrintErrMsg(option)
  end if
  ! initialize constraints
  option%io_buffer = "initializing constraint : " // tran_constraint%name
  call PrintMsg(option)

  call ReactionProcessConstraint(reaction,tran_constraint,option)

  ! equilibrate
  option%io_buffer = "equilibrate constraint : " // &
    constraint_coupler%constraint%name
  call PrintMsg(option)
  call ReactionEquilibrateConstraint(rt_auxvar, global_auxvar, &
         material_auxvar, reaction, &
         tran_constraint, &
         num_iterations, &
         use_prev_soln_as_guess, &
         option)

  ! Pflotran does not seem to apply the constraint to immobile species in
  ! ReactionEquilibrateConstraint
  ! So we will apply it by hand here
  do i = 1, reaction%immobile%nimmobile
    rt_auxvar%immobile(i) = tran_constraint%immobile_species%constraint_conc(i)
  enddo

  ! link the constraint to the constraint coupler so we can print it
  constraint_coupler%global_auxvar => global_auxvar
  constraint_coupler%rt_auxvar => rt_auxvar
  constraint_coupler%num_iterations = num_iterations

  call ReactionPrintConstraint(constraint_coupler, reaction, option)

end subroutine ProcessPFLOTRANConstraint


! **************************************************************************** !
function ConvertAlquimiaConditionToPflotran(&
     option, reaction, alquimia_condition)
  use, intrinsic :: iso_c_binding
  use petscsys
  use c_f_interface_module, only : c_f_string_ptr

  use AlquimiaContainers_module

  ! pflotran
  use Option_module, only : option_type, PrintErrMsg, PrintMsg
  use Reaction_aux_module, only : reaction_rt_type, &
        aq_species_constraint_type, AqueousSpeciesConstraintCreate
  use Reaction_Mineral_Aux_module, only : mineral_constraint_type, &
        MineralConstraintCreate
  use String_module, only : StringCompareIgnoreCase
  use Transport_Constraint_RT_module, only : tran_constraint_rt_type, &
        TranConstraintRTCreate, &
        CONSTRAINT_FREE, CONSTRAINT_TOTAL, CONSTRAINT_TOTAL_SORB, &
        CONSTRAINT_PH, CONSTRAINT_MINERAL, &
        CONSTRAINT_GAS, CONSTRAINT_CHARGE_BAL, CONSTRAINT_TOTAL_AQ_PLUS_SORB
  use Reaction_Immobile_Aux_module, only : immobile_constraint_type, ImmobileConstraintCreate
  use petscsys

  implicit none

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  class(reaction_rt_type), pointer, intent(in) :: reaction
  type (AlquimiaGeochemicalCondition), intent(in) :: alquimia_condition

  ! Return value
  class (tran_constraint_rt_type), pointer :: ConvertAlquimiaConditionToPflotran

  ! local variables
  integer :: i,i_all,i_immobile,jimmobile
  character (kAlquimiaMaxStringLength) :: name, constraint_type
  character (kAlquimiaMaxStringLength) :: associated_species
  class (tran_constraint_rt_type), pointer :: tran_constraint
  type(aq_species_constraint_type), pointer :: pft_aq_species_constraint
  type(mineral_constraint_type), pointer :: pft_mineral_constraint
  type(immobile_constraint_type), pointer :: pft_immobile_species_constraint
  type (AlquimiaAqueousConstraint), pointer :: alq_aqueous_constraints(:)
  type (AlquimiaMineralConstraint), pointer :: alq_mineral_constraints(:)
  logical :: is_immobile

  call c_f_string_ptr(alquimia_condition%name, name)
  option%io_buffer = "building : " // trim(name)
  call PrintMsg(option)
  tran_constraint => TranConstraintRTCreate(option)
  tran_constraint%name = trim(name)
  ! NOTE(bja): requires_equilibration not used in pflotran?
  tran_constraint%equilibrate_at_each_cell = PETSC_FALSE

  !
  ! aqueous species
  !
  if (alquimia_condition%aqueous_constraints%size /= reaction%naqcomp+reaction%immobile%nimmobile) then
     option%io_buffer = 'Number of aqueous constraints ' // &
          'does not equal the number of primary chemical ' // &
          'components in constraint: ' // &
          trim(tran_constraint%name)
     call PrintErrMsg(option)
  end if

  ! NOTE(bja) : this is the container for ALL aqueous constraints
  pft_aq_species_constraint => &
       AqueousSpeciesConstraintCreate(reaction, option)

  pft_immobile_species_constraint => &
      ImmobileConstraintCreate(reaction%immobile, option)

  call c_f_pointer(alquimia_condition%aqueous_constraints%data, &
       alq_aqueous_constraints, (/alquimia_condition%aqueous_constraints%size/))

  i_immobile=0
  i=0
  do i_all = 1, alquimia_condition%aqueous_constraints%size
   call c_f_string_ptr(alq_aqueous_constraints(i_all)%primary_species_name, name)

   ! Check if it is an immobile species
   is_immobile = .FALSE.
   do jimmobile = 1, reaction%immobile%nimmobile
      if (StringCompareIgnoreCase(name, &
                        reaction%immobile%names(jimmobile), &
                        MAXWORDLENGTH)) then
         is_immobile = .TRUE.
         exit
      endif
   enddo

   if(is_immobile) then
      i_immobile=i_immobile+1
      pft_immobile_species_constraint%names(i_immobile) = trim(name)
      pft_immobile_species_constraint%constraint_conc(i_immobile) = alq_aqueous_constraints(i_all)%value
   else
      i=i+1
      pft_aq_species_constraint%names(i) = trim(name)

     pft_aq_species_constraint%constraint_conc(i) = alq_aqueous_constraints(i_all)%value

     call c_f_string_ptr(alq_aqueous_constraints(i_all)%constraint_type, constraint_type)

     call c_f_string_ptr(alq_aqueous_constraints(i_all)%associated_species, &
          associated_species)

     if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringFree)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_FREE
     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringTotalAqueousPlusSorbed)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_TOTAL_AQ_PLUS_SORB

     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringTotalAqueous)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_TOTAL

     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringTotalSorbed)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_TOTAL_SORB
     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringPH)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_PH

     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringMineral)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_MINERAL
        pft_aq_species_constraint%constraint_aux_string(i) = &
             trim(associated_species)

     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringGas)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_GAS
        pft_aq_species_constraint%constraint_aux_string(i) = &
             trim(associated_species)

     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringCharge)) then
        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_CHARGE_BAL

     else
        option%io_buffer = 'Constraint type: ' // trim(constraint_type) // &
             ' not recognized in constraint,concentration'
        call PrintErrMsg(option)
     end if
   end if
  end do
  tran_constraint%aqueous_species => pft_aq_species_constraint
  tran_constraint%immobile_species => pft_immobile_species_constraint

  !
  ! minerals
  !
  ! FIXME(bja): are these checks the correct thing to do in all cases...?
  if (alquimia_condition%mineral_constraints%size > 0 .and. &
       alquimia_condition%mineral_constraints%size /= &
       reaction%mineral%nkinmnrl) then
     option%io_buffer = &
          'Number of mineral constraints is not equal to ' // &
          'number of kinetic minerals in condition: ' // &
          trim(tran_constraint%name)
     call PrintErrMsg(option)
  end if

  pft_mineral_constraint => MineralConstraintCreate(reaction%mineral, option)

  call c_f_pointer(alquimia_condition%mineral_constraints%data, &
       alq_mineral_constraints, (/alquimia_condition%mineral_constraints%size/))
  do i = 1, alquimia_condition%mineral_constraints%size
     call c_f_string_ptr(alq_mineral_constraints(i)%mineral_name, name)
     pft_mineral_constraint%names(i) = trim(name)
     pft_mineral_constraint%constraint_vol_frac(i) = &
          alq_mineral_constraints(i)%volume_fraction
     pft_mineral_constraint%constraint_area(i) = &
          alq_mineral_constraints(i)%specific_surface_area
      pft_mineral_constraint%constraint_area_units(i) = 'm^2/m^3'
  end do
  tran_constraint%minerals => pft_mineral_constraint

  ConvertAlquimiaConditionToPflotran => tran_constraint

end function ConvertAlquimiaConditionToPflotran

! **************************************************************************** !
subroutine CopyAlquimiaToAuxVars(copy_auxdata, hands_off, &
  state, aux_data, prop, &
  reaction, global_auxvar, material_auxvar, rt_auxvar)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_aux_module, only : reaction_rt_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Material_Aux_class, only : material_auxvar_type

  implicit none

  ! function parameters
  logical, intent(in) :: copy_auxdata
  logical, intent(in) :: hands_off
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaProperties), intent(in) :: prop
  class(reaction_rt_type), pointer, intent(inout) :: reaction
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvar
  class(material_auxvar_type), pointer, intent(inout) :: material_auxvar
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvar

  ! local variables
  real (c_double), pointer :: data(:)
  integer :: i
  PetscInt, parameter :: phase_index = 1

  !write (*, '(a)') "PFLOTRAN_Alquimia_CopyAlquimiaToAuxVars() :"

  !
  ! state
  !
  global_auxvar%den_kg(1) = state%water_density
  global_auxvar%sat(1) = prop%saturation
  global_auxvar%temp = state%temperature
  global_auxvar%pres(1) = state%aqueous_pressure

  material_auxvar%porosity = state%porosity
  material_auxvar%volume = prop%volume

  !
  ! primary aqueous
  !
  call c_f_pointer(state%total_mobile%data, data, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvar%total(i, phase_index) = data(i)
  end do

  ! sorbed primary
  if (reaction%neqsorb > 0) then
     call c_f_pointer(state%total_immobile%data, data, (/reaction%naqcomp/))
     do i = 1, reaction%naqcomp
        rt_auxvar%total_sorb_eq(i) = data(i)
     end do
  end if
  
  ! immobile species
  if (reaction%immobile%nimmobile > 0) then
     call c_f_pointer(state%total_immobile%data, data, (/reaction%ncomp/))
     do i = 1, reaction%immobile%nimmobile
        rt_auxvar%immobile(i) = data(i + reaction%offset_immobile)
     end do
  end if

  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction%data, data, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     rt_auxvar%mnrl_volfrac(i) = data(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area%data, data, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     rt_auxvar%mnrl_area(i) = data(i)
  end do

  !
  ! in hands-off mode CEC and site density, as well as any property
  ! are not provided by the driver so copying them over would
  ! lose pflotran's input file values 
  !
  if_hands_off: if (.not. hands_off) then
  
  !
  ! ion exchange, CEC only present in reaction, not aux_vars?
  !
  call c_f_pointer(state%cation_exchange_capacity%data, data, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     reaction%eqionx_rxn_CEC(i) = data(i)
  end do

  !
  ! equilibrium surface complexation
  !
  call c_f_pointer(state%surface_site_density%data, data, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     reaction%surface_complexation%srfcplxrxn_site_density(i) = data(i)
  end do

  !
  ! isotherms
  !
  call c_f_pointer(prop%isotherm_kd%data, data, &
       (/prop%isotherm_kd%size/))
  do i = 1, prop%isotherm_kd%size
     reaction%isotherm%isotherm_rxn%eqisothermcoeff(i) = data(i)
  end do

  call c_f_pointer(prop%langmuir_b%data, data, &
       (/prop%langmuir_b%size/))
  do i = 1, prop%langmuir_b%size
     reaction%isotherm%isotherm_rxn%eqisothermlangmuirb(i) = data(i)
  end do

  call c_f_pointer(prop%freundlich_n%data, data, &
       (/prop%freundlich_n%size/))
  do i = 1, prop%freundlich_n%size
     reaction%isotherm%isotherm_rxn%eqisothermfreundlichn(i) = data(i)
  end do

  !
  ! mineral reaction rate constant
  !
  call c_f_pointer(prop%mineral_rate_cnst%data, data, &
       (/prop%mineral_rate_cnst%size/))
  do i = 1, prop%mineral_rate_cnst%size
     reaction%mineral%kinmnrl_rate_constant(i) = data(i)
  end do

  !
  ! aqueous kinetic reaction rate constant
  !
  call c_f_pointer(prop%aqueous_kinetic_rate_cnst%data, data, &
       (/prop%aqueous_kinetic_rate_cnst%size/))
  do i = 1, reaction%ngeneral_rxn
      reaction%general_kf(i) = data(i)
      reaction%general_kr(i) = 0.0d0
  end do
  do i=1, reaction%microbial%nrxn
    reaction%microbial%rate_constant(i) = data(i+reaction%ngeneral_rxn)
  end do

  end if if_hands_off
  
  if (copy_auxdata) then
     call UnpackAlquimiaAuxiliaryData(aux_data, reaction, rt_auxvar)
  end if

end subroutine CopyAlquimiaToAuxVars

! **************************************************************************** !
subroutine CopyAuxVarsToAlquimia(reaction, global_auxvar, rt_auxvar, &
     porosity, state, aux_data)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_aux_module, only : reaction_rt_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type

  implicit none

  ! function parameters
  class(reaction_rt_type), pointer, intent(in) :: reaction
  type(global_auxvar_type), pointer, intent(in) :: global_auxvar
  type(reactive_transport_auxvar_type), pointer, intent(in) :: rt_auxvar
  PetscReal, intent(in) :: porosity
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  real (c_double), pointer :: data(:)
  integer :: i
  PetscInt, parameter :: phase_index = 1

  !write (*, '(a)') "PFLOTRAN_Alquimia_CopyAuxVarsToAlquimia() :"

  !
  ! state
  !
  state%water_density = global_auxvar%den_kg(1)
  state%temperature = global_auxvar%temp
  state%aqueous_pressure = global_auxvar%pres(1)

  state%porosity = porosity

  !
  ! primary aqueous species
  !
  call c_f_pointer(state%total_mobile%data, data, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     data(i) = rt_auxvar%total(i, phase_index)
  end do

  !
  ! sorbed
  !
  if (reaction%neqsorb > 0) then
     call c_f_pointer(state%total_immobile%data, data, (/reaction%naqcomp/))
     do i = 1, reaction%naqcomp
        data(i) = rt_auxvar%total_sorb_eq(i)
     end do
  end if
  !
  ! immobile species
  !
  if (reaction%immobile%nimmobile > 0) then
     call c_f_pointer(state%total_immobile%data, data, (/reaction%ncomp/))
     do i = 1, reaction%immobile%nimmobile
        data(i + reaction%offset_immobile) = rt_auxvar%immobile(i)
     end do
  end if
  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction%data, data, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     data(i) = rt_auxvar%mnrl_volfrac(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area%data, data, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     data(i) = rt_auxvar%mnrl_area(i)
  end do

  !
  ! ion exchange, CEC only present in reaction, not aux_vars?
  !
  call c_f_pointer(state%cation_exchange_capacity%data, data, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     data(i) = reaction%eqionx_rxn_CEC(i)
  end do

  !
  ! equilibrium surface complexation, site density in reaction, not aux_vars
  !
  call c_f_pointer(state%surface_site_density%data, data, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     data(i) = reaction%surface_complexation%srfcplxrxn_site_density(i)
  end do

  ! NOTE(bja): isotherms are material properties, and can't be changed
  ! by chemistry. We don't need to copy theme here!

  ! NOTE(smr): reaction rates constants are properties, and can't be changed
  ! by chemistry. We don't need to copy theme here!

  call PackAlquimiaAuxiliaryData(reaction, rt_auxvar, aux_data)

end subroutine CopyAuxVarsToAlquimia

! **************************************************************************** !
!
! Auxiliary Data helper routines
!
! We are packaging data into the auxiliary data vectors in the following order:
!
! ints: none
!
! doubles:
!   free ion conc <N_primary>
!   primary_activity_coeff <N_primary>
!   secondary_activity_coeff <N_aqueous_complexes>
!   ion exchange ref cation sorbed conc <N_ion_exchange_sites>
!   surface complexation free site conc <N_surface_sites>
!
! **************************************************************************** !
subroutine GetAuxiliaryDataSizes(reaction, num_ints, num_doubles)

  use, intrinsic :: iso_c_binding, only : c_int

  use Reaction_aux_module, only : reaction_rt_type

  class (reaction_rt_type), intent(in) :: reaction
  integer (c_int), intent(out) :: num_ints
  integer (c_int), intent(out) :: num_doubles

  num_ints = 0

  num_doubles = &
       2 * reaction%ncomp + &
       reaction%neqcplx + &
       reaction%neqionxrxn + &
       reaction%surface_complexation%nsrfcplxrxn

end subroutine GetAuxiliaryDataSizes

! **************************************************************************** !
subroutine PackAlquimiaAuxiliaryData(reaction, rt_auxvar, aux_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData

  use Reaction_aux_module, only : reaction_rt_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type

  ! function parameters
  class (reaction_rt_type), intent(in) :: reaction
  type(reactive_transport_auxvar_type), pointer, intent(in) :: rt_auxvar
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  integer (c_int) :: num_ints, num_doubles
  real (c_double), pointer :: data(:)
  integer :: i, dindex

  call GetAuxiliaryDataSizes(reaction, num_ints, num_doubles)

  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0
  ! free ion
  do i = 1, reaction%naqcomp
     dindex = dindex + 1
     data(dindex) = rt_auxvar%pri_molal(i)
  end do
  
  ! Fake these (0) for immobile species
  ! free ion
  do i = 1, reaction%immobile%nimmobile
     dindex = dindex + 1
     data(dindex) = 0.0
  end do

  ! primary activity coeff
  do i = 1, reaction%naqcomp
     dindex = dindex + 1
     data(dindex) = rt_auxvar%pri_act_coef(i)
  end do
  
  ! Fake it for immobile species
  do i = 1, reaction%immobile%nimmobile
     dindex = dindex + 1
     data(dindex) = 0.0
  end do

  ! secondary aqueous complexe activity coeffs
  do i = 1, reaction%neqcplx
     dindex = dindex + 1
     data(dindex) = rt_auxvar%sec_act_coef(i)
  end do

  ! ion exchange
  do i = 1, reaction%neqionxrxn
     dindex = dindex + 1
     data(dindex) = rt_auxvar%eqionx_ref_cation_sorbed_conc(i)
  end do

  ! equilibrium surface complexation
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     dindex = dindex + 1
     data(dindex) = rt_auxvar%srfcplxrxn_free_site_conc(i)
  end do

end subroutine PackAlquimiaAuxiliaryData

! **************************************************************************** !
subroutine UnpackAlquimiaAuxiliaryData(aux_data, reaction, rt_auxvar)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData

  use Reaction_aux_module, only : reaction_rt_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type

  ! function parameters
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  class (reaction_rt_type), intent(inout) :: reaction
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvar

  ! local variables
  integer (c_int) :: num_ints, num_doubles
  real (c_double), pointer :: data(:)
  integer :: i, dindex

  call GetAuxiliaryDataSizes(reaction, num_ints, num_doubles)

  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0

  ! free ion
  do i = 1, reaction%naqcomp
     dindex = dindex + 1
     rt_auxvar%pri_molal(i) = data(dindex)
  end do

  ! Skip immobile species, which driver model thinks are primary species
  dindex = dindex + reaction%immobile%nimmobile

  ! primary activity coeff
  do i = 1, reaction%naqcomp
     dindex = dindex + 1
     rt_auxvar%pri_act_coef(i) = data(dindex)
  end do
  
  dindex = dindex + reaction%immobile%nimmobile

  ! aqueous complexes activity coeff
  do i = 1, reaction%neqcplx
     dindex = dindex + 1
     rt_auxvar%sec_act_coef(i) = data(dindex)
  end do


  ! ion exchange
  do i = 1, reaction%neqionxrxn
     dindex = dindex + 1
     rt_auxvar%eqionx_ref_cation_sorbed_conc(i) = data(dindex)
  end do


  ! equilibrium surface complexation
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     dindex = dindex + 1
     rt_auxvar%srfcplxrxn_free_site_conc(i) = data(dindex)
  end do

end subroutine UnpackAlquimiaAuxiliaryData

! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PrintSizes(sizes)

  use AlquimiaContainers_module, only : AlquimiaSizes

  implicit none

  ! function parameters
  type (AlquimiaSizes), intent(in) :: sizes

  write (*, '(a)') "size : "
  write (*, '(a, i4)') "  num primary : ", sizes%num_primary
  write (*, '(a, i4)') "  num kinetics minerals : ", sizes%num_minerals
  write (*, '(a, i4)') "  num aqueous complexes : ", sizes%num_aqueous_complexes
  write (*, '(a, i4)') "  num surface sites : ", sizes%num_surface_sites
  write (*, '(a, i4)') "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PrintSizes


! **************************************************************************** !
subroutine PrintState(state)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaState

  implicit none

  ! function parameters
  type (AlquimiaState), intent(in) :: state

  ! local variables
  integer :: i
  real (c_double), pointer :: conc(:)
  real (c_double), pointer :: immo(:)
  real (c_double), pointer :: vofx(:)
  real (c_double), pointer :: surf(:)
  real (c_double), pointer :: sites(:)
  real (c_double), pointer :: cec(:)

  write (*, '(a)') "state : "
  write (*, '(a, 1es13.6)') "  density water : ", state%water_density
  write (*, '(a, 1es13.6)') "  porosity : ", state%porosity
  write (*, '(a, 1es13.6)') "  temperature : ", state%temperature
  write (*, '(a, 1es13.6)') "  aqueous pressure : ", state%aqueous_pressure
  write (*, '(a, i4, a)') "  total mobile (", state%total_mobile%size, ") : "
  call c_f_pointer(state%total_mobile%data, conc, (/state%total_mobile%size/))
  do i=1, state%total_mobile%size
     write (*, '(1es13.6)') conc(i)
  end do
    write (*, '(a, i4, a)') "  total immobile (", state%total_immobile%size, ") : "
  call c_f_pointer(state%total_immobile%data, immo, (/state%total_immobile%size/))
  do i=1, state%total_immobile%size
     write (*, '(1es13.6)') immo(i)
  end do
  write (*, '(a, i4, a)') "  mineral volume fraction (", state%mineral_volume_fraction%size, ") : "
  call c_f_pointer(state%mineral_volume_fraction%data, vofx, (/state%mineral_volume_fraction%size/))
  do i=1, state%mineral_volume_fraction%size
     write (*, '(1es13.6)') vofx(i)
  end do
  write (*, '(a, i4, a)') "  mineral surface area (", state%mineral_specific_surface_area%size, ") : "
  call c_f_pointer(state%mineral_specific_surface_area%data, surf, (/state%mineral_specific_surface_area%size/))
  do i=1, state%mineral_specific_surface_area%size
     write (*, '(1es13.6)') surf(i)
  end do
  write (*, '(a, i4, a)') "  surface site density (", state%surface_site_density%size, ") : "
  call c_f_pointer(state%surface_site_density%data, sites, (/state%surface_site_density%size/))
  do i=1, state%surface_site_density%size
     write (*, '(1es13.6)') sites(i)
  end do
  write (*, '(a, i4, a)') "  cation exchange capacity (", state%cation_exchange_capacity%size, ") : "
  call c_f_pointer(state%cation_exchange_capacity%data, cec, (/state%cation_exchange_capacity%size/))
  do i=1, state%cation_exchange_capacity%size
     write (*, '(1es13.6)') cec(i)
  end do
end subroutine PrintState


! **************************************************************************** !
subroutine PrintEngineFunctionality(functionality)

  use AlquimiaContainers_module, only : AlquimiaEngineFunctionality

  implicit none

  ! function parameters
  type (AlquimiaEngineFunctionality), intent(in) :: functionality

  write (*, '(a)') "functionality : "
  write (*, '(a, L1)') "  thread safe : ", functionality%thread_safe
  write (*, '(a, L1)') "  temperature dependent : ", functionality%temperature_dependent
  write (*, '(a, L1)') "  pressure dependent : ", functionality%pressure_dependent
  write (*, '(a, L1)') "  porosity update : ", functionality%porosity_update
  write (*, '(a, i4)') "  index base : ", functionality%index_base
end subroutine PrintEngineFunctionality

! **************************************************************************** !
subroutine PrintProblemMetaData(meta_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer

  use c_f_interface_module, only : c_f_string_ptr

  use AlquimiaContainers_module, only : &
       AlquimiaProblemMetaData, kAlquimiaMaxStringLength

  implicit none

  ! function parameters
  type (AlquimiaProblemMetaData), intent(in) :: meta_data

  ! local variables
  type (c_ptr), pointer :: names(:)
  character(len=kAlquimiaMaxStringLength) :: name
  integer (c_int) :: i

  write (*, '(a)') "meta_data : "
  write (*, '(a)') "  primary names : "
  call c_f_pointer(meta_data%primary_names%data, names, (/meta_data%primary_names%size/))
  do i=1, meta_data%primary_names%size
     call c_f_string_ptr(names(i), name)
     write (*, '(a)') trim(name)
  end do
end subroutine PrintProblemMetaData


! **************************************************************************** !
subroutine PrintStatus(status)

  use AlquimiaContainers_module, only : AlquimiaEngineStatus

  implicit none

  ! function parameters
  type (AlquimiaEngineStatus), intent(in) :: status

  write (*, '(a)') "status : "
  write (*, '(a)') "  num rhs evaluation  : ", status%num_rhs_evaluations
  write (*, '(a)') "  num jacobian evaluations : ", status%num_jacobian_evaluations
  write (*, '(a)') "  num newton iterations : ", status%num_newton_iterations
  write (*, '(a)') "  converged  : ", status%converged

end subroutine PrintStatus


! **************************************************************************** !
subroutine PrintTranConstraint(tran_constraint)

  use Transport_Constraint_RT_module, only : tran_constraint_rt_type

  implicit none

  ! function parameters
  class (tran_constraint_rt_type), pointer :: tran_constraint

  write (*, '(a)') "TranConstraint :"
  write (*, '(a, i4)') "    id : ", tran_constraint%id
  write (*, '(a, a)') "    name : ", tran_constraint%name
  write (*, '(a, L1)') "    requires equilibration : ", tran_constraint%equilibrate_at_each_cell
  call PrintAqueousSpeciesConstraint(tran_constraint%aqueous_species)
  call PrintMineralConstraint(tran_constraint%minerals)
!    type(mineral_constraint_type), pointer :: minerals
!    type(srfcplx_constraint_type), pointer :: surface_complexes
!    type(colloid_constraint_type), pointer :: colloids
!    type(immobile_constraint_type), pointer :: immobile_species


end subroutine PrintTranConstraint

subroutine PrintAqueousSpeciesConstraint(aqueous_species)

  use Reaction_aux_module, only : aq_species_constraint_type

  implicit none

  type (aq_species_constraint_type), pointer :: aqueous_species 

  write (*, '(a)') "    Aqueous species :"
  write (*, '(a)') aqueous_species%names
  write (*, '(a)') "    Constraint Conc :"
  write (*, '(f18.8)') aqueous_species%constraint_conc(:)
  write (*, '(a)') "    Constraint basis molarity :"
  write (*, '(f18.8)') aqueous_species%basis_molarity(:)
  write (*, '(a)') "    Constraint type :"
  write (*, '(i4)') aqueous_species%constraint_type(:)
  write (*, '(a)') "    Constraint spec id :"
  write (*, '(i4)') aqueous_species%constraint_spec_id(:)
  write (*, '(a)') "    Constraint aux string :"
  write (*, '(a)') aqueous_species%constraint_aux_string(:)

end subroutine PrintAqueousSpeciesConstraint

subroutine PrintMineralConstraint(minerals)

  use Reaction_Mineral_Aux_module, only : mineral_constraint_type

  implicit none

  type (mineral_constraint_type), pointer :: minerals

  write (*, '(a)') "    Mineral species :"
  write (*, '(a)') minerals%names
  write (*, '(a)') "    volume fraction :"
  write (*, '(f18.8)') minerals%constraint_vol_frac(:)
  write (*, '(a)') "    Constraint area :"
  write (*, '(f18.8)') minerals%constraint_area(:)
  !This member no longer exists!
  !write (*, '(a)') "    Constraint aux string :"
  !write (*, '(a)') minerals%constraint_aux_string(:)

end subroutine PrintMineralConstraint


end module PFLOTRANAlquimiaInterface_module


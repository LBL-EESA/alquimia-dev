
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
!          Sergi Molins <smolins@lbl.gov>
!

! **************************************************************************** !
!
! CrunchFlow Alquimia Inteface Wrappers
!
! Author: Sergi Molins
!
! Different fortran compilers use different name mangling conventions
! for fortran modules:
!
!    gfortran : ___modulename_MOD_procedurename
!
!    intel : _modulename_mp_procedurename_
!
!    as a consequence we can't put the alquimia interface into a
!    module and call it directly from C/C++. Instead we use
!    some simple wrapper functions.
!
! Notes:
!
!  * Function call signatures are dictated by the alquimia API!
!
!  * alquimia data structures defined in AlquimiaContainers_module
!    (alquimia_containers.F90) are dictated by the alquimia API.
!
! **************************************************************************** !


! **************************************************************************** !
subroutine Crunch_Alquimia_Setup(input_filename, hands_off, cf_engine_state, &
     sizes, functionality, status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use CrunchAlquimiaInterface_module, only : Setup

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  logical (c_bool), value, intent(in) :: hands_off
  type (c_ptr), intent(out) :: cf_engine_state
  type (AlquimiaSizes), intent(out) :: sizes
  type (AlquimiaEngineFunctionality), intent(out) :: functionality
  type (AlquimiaEngineStatus), intent(out) :: status

  call Setup(input_filename, hands_off, cf_engine_state, sizes, functionality, status)

end subroutine Crunch_Alquimia_Setup


! **************************************************************************** !
subroutine Crunch_Alquimia_Shutdown(cf_engine_state, status) bind(c)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use CrunchAlquimiaInterface_module, only : Shutdown

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaEngineStatus), intent(out) :: status

  call Shutdown(cf_engine_state, status)

end subroutine Crunch_Alquimia_Shutdown


! **************************************************************************** !
subroutine Crunch_Alquimia_ProcessCondition( &
     cf_engine_state, &
     condition, &
     properties, &
     state, &
     aux_data, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use CrunchAlquimiaInterface_module, only : ProcessCondition

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaGeochemicalCondition), intent(in) :: condition
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent (inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  call ProcessCondition(cf_engine_state, condition, properties, &
       state, aux_data, status)

end subroutine Crunch_Alquimia_ProcessCondition


! **************************************************************************** !
subroutine Crunch_Alquimia_ReactionStepOperatorSplit( &
     cf_engine_state, &
     delta_t, &
     properties, &
     state, &
     aux_data, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use CrunchAlquimiaInterface_module, only : ReactionStepOperatorSplit

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  real (c_double), value, intent(in) :: delta_t
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  call ReactionStepOperatorSplit(cf_engine_state, delta_t, &
       properties, state, aux_data, status)

end subroutine Crunch_Alquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine Crunch_Alquimia_GetAuxiliaryOutput( &
     cf_engine_state, &
     properties, &
     state, &
     aux_data, &
     aux_output, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use CrunchAlquimiaInterface_module, only : GetAuxiliaryOutput

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaAuxiliaryOutputData), intent(inout) :: aux_output
  type (AlquimiaEngineStatus), intent(out) :: status

  call GetAuxiliaryOutput(cf_engine_state, &
       properties, state, aux_data, aux_output, status)

end subroutine Crunch_Alquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine Crunch_Alquimia_GetProblemMetaData(cf_engine_state, &
     meta_data, status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use CrunchAlquimiaInterface_module, only: GetProblemMetaData

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaProblemMetaData), intent(out) :: meta_data
  type (AlquimiaEngineStatus), intent(out) :: status

  call GetProblemMetaData(cf_engine_state, meta_data, status)

end subroutine Crunch_Alquimia_GetProblemMetaData

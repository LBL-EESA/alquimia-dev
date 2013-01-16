! **************************************************************************** !
!
! PFloTran Alquimia Inteface Wrappers
!
! Author: Benjamin Andre
!
! Notes:
!
!  * Function call signatures are dictated by the alquimia API!
!
!  * alquimia data structures defined in alquimia_containers.h90 are
!    dictated by the alquimia API.
!
!  * (bja) 2012-12 - different fortran compilers use different name
!    mangling conventions for fortran modules:
!
!    gfortran : ___modulename_MOD_procedurename
!
!    intel : _modulename_mp_procedurename_
!
!    as a consequence we can't put the alquimia interface into a
!    module and call it directly from C/C++. Instead we use
!    some simple wrapper functions.
!
! **************************************************************************** !


! **************************************************************************** !
subroutine PFloTran_Alquimia_Setup(input_filename, pft_engine_state, sizes) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (c_ptr), intent(out) :: pft_engine_state
  type (alquimia_sizes_f), intent(out) :: sizes

  call Setup(input_filename, pft_engine_state, sizes)

end subroutine PFloTran_Alquimia_Setup


! **************************************************************************** !
subroutine PFloTran_Alquimia_Shutdown(pft_engine_state) bind(c)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state

  call Shutdown(pft_engine_state)

end subroutine PFloTran_Alquimia_Shutdown


! **************************************************************************** !
subroutine PFloTran_Alquimia_ProcessCondition(pft_engine_state, condition, &
     sizes, state, status) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (alquimia_condition_f), intent(in) :: condition
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_state_f), intent(inout) :: state
  type (alquimia_engine_status_f), intent(out) :: status

  call ProcessCondition(pft_engine_state, condition, sizes, state, status)

end subroutine PFloTran_Alquimia_ProcessCondition


! **************************************************************************** !
subroutine PFloTran_Alquimia_ReactionStepOperatorSplit(pft_engine_state, &
     delta_t, material_properties, state, aux_data, status) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  real (c_double), intent(in) :: delta_t
  type (alquimia_material_properties_f), intent(in) :: material_properties
  type (alquimia_state_f), intent(inout) :: state
  type (alquimia_auxiliary_data_f), intent(inout) :: aux_data
  type (alquimia_engine_status_f), intent(out) :: status

  call ReactionStepOperatorSplit(pft_engine_state, delta_t, &
       material_properties, state, aux_data, status)

end subroutine PFloTran_Alquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetAuxiliaryOutput(pft_engine_state) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state

  call GetAuxiliaryOutput(pft_engine_state)

end subroutine PFloTran_Alquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetEngineMetaData(pft_engine_state, &
     sizes, meta_data) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_meta_data_f), intent(out) :: meta_data

  call GetEngineMetaData(pft_engine_state, sizes, meta_data)

end subroutine PFloTran_Alquimia_GetEngineMetaData

! **************************************************************************** !
subroutine PFloTran_Alquimia_GetPrimaryNameFromIndex(pft_engine_state, &
  primary_index, primary_name) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  integer (c_int), intent(in) :: primary_index
  character(kind=c_char), dimension(*), intent(out) :: primary_name

  call GetPrimaryNameFromIndex(pft_engine_state, primary_index, primary_name)

end subroutine PFloTran_Alquimia_GetPrimaryNameFromIndex


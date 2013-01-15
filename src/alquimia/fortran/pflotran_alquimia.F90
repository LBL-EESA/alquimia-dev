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
subroutine PFloTran_Alquimia_Setup(pft_internal_state, input_filename, sizes) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (alquimia_sizes_f), intent(inout) :: sizes
  type (c_ptr), intent(inout) :: pft_internal_state

  call Setup(pft_internal_state, input_filename, sizes)

end subroutine PFloTran_Alquimia_Setup


! **************************************************************************** !
subroutine PFloTran_Alquimia_Shutdown(pft_internal_state) bind(c)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_internal_state

  call Shutdown(pft_internal_state)

end subroutine PFloTran_Alquimia_Shutdown


! **************************************************************************** !
subroutine PFloTran_Alquimia_ProcessCondition(pft_internal_state, condition, &
     sizes, state) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_internal_state
  type (alquimia_condition_f), intent(in) :: condition
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_state_f), intent(inout) :: state

  call ProcessCondition(pft_internal_state, condition, sizes, state)

end subroutine PFloTran_Alquimia_ProcessCondition


! **************************************************************************** !
subroutine PFloTran_Alquimia_ReactionStepOperatorSplit(pft_internal_state) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_internal_state

  call ReactionStepOperatorSplit(pft_internal_state)

end subroutine PFloTran_Alquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetAuxiliaryOutput(pft_internal_state) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_internal_state

  call GetAuxiliaryOutput(pft_internal_state)

end subroutine PFloTran_Alquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetEngineMetaData(pft_internal_state, &
     sizes, meta_data) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_internal_state
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_meta_data_f), intent(out) :: meta_data

  call GetEngineMetaData(pft_internal_state, sizes, meta_data)

end subroutine PFloTran_Alquimia_GetEngineMetaData

! **************************************************************************** !
subroutine PFloTran_Alquimia_GetPrimaryNameFromIndex(pft_internal_state, &
  primary_index, primary_name) bind(C)

  use, intrinsic :: iso_c_binding

  use PFloTranAlquimiaInterface_module

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_internal_state
  integer (c_int), intent(in) :: primary_index
  character(kind=c_char), dimension(*), intent(out) :: primary_name

  call GetPrimaryNameFromIndex(pft_internal_state, primary_index, primary_name)

end subroutine PFloTran_Alquimia_GetPrimaryNameFromIndex


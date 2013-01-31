! **************************************************************************** !
!
! PFloTran Alquimia Inteface Wrappers
!
! Author: Benjamin Andre
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
subroutine PFloTran_Alquimia_Setup(input_filename, pft_engine_state, &
     sizes, status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use PFloTranAlquimiaInterface_module

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (c_ptr), intent(out) :: pft_engine_state
  type (AlquimiaSizes), intent(out) :: sizes
  type (AlquimiaEngineStatus), intent(out) :: status

  call Setup(input_filename, pft_engine_state, sizes, status)

end subroutine PFloTran_Alquimia_Setup


! **************************************************************************** !
subroutine PFloTran_Alquimia_Shutdown(pft_engine_state, status) bind(c)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use PFloTranAlquimiaInterface_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaEngineStatus), intent(out) :: status

  call Shutdown(pft_engine_state, status)

end subroutine PFloTran_Alquimia_Shutdown


! **************************************************************************** !
subroutine PFloTran_Alquimia_ProcessCondition( &
     pft_engine_state, &
     condition, &
     material_properties, &
     state, &
     aux_data, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use PFloTranAlquimiaInterface_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaGeochemicalCondition), intent(in) :: condition
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent (inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  call ProcessCondition(pft_engine_state, condition, material_properties, &
       state, aux_data, status)

end subroutine PFloTran_Alquimia_ProcessCondition


! **************************************************************************** !
subroutine PFloTran_Alquimia_ReactionStepOperatorSplit( &
     pft_engine_state, &
     delta_t, &
     material_properties, &
     state, &
     aux_data, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use PFloTranAlquimiaInterface_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  real (c_double), intent(in) :: delta_t
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  call ReactionStepOperatorSplit(pft_engine_state, delta_t, &
       material_properties, state, aux_data, status)

end subroutine PFloTran_Alquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetAuxiliaryOutput( &
     pft_engine_state, &
     material_properties, &
     state, &
     aux_data, &
     aux_output, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use PFloTranAlquimiaInterface_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaAuxiliaryOutputData), intent(inout) :: aux_output
  type (AlquimiaEngineStatus), intent(out) :: status

  call GetAuxiliaryOutput(pft_engine_state, &
       material_properties, state, aux_data, aux_output, status)

end subroutine PFloTran_Alquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetEngineMetaData(pft_engine_state, &
     meta_data, status) bind(C)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module
  use PFloTranAlquimiaInterface_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaMetaData), intent(out) :: meta_data
  type (AlquimiaEngineStatus), intent(out) :: status

  call GetEngineMetaData(pft_engine_state, meta_data, status)

end subroutine PFloTran_Alquimia_GetEngineMetaData



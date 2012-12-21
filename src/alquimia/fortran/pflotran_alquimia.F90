! **************************************************************************** !
!
! Alquimia inteface
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PFloTranAlquimia_Setup(input_filename, sizes) bind(C)

#include "alquimia_containers.h90"

  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (alquimia_sizes_f), intent(inout) :: sizes
  integer :: len
  integer :: i

  print *, "PFloTranAlquimia_Setup() : "
  len = 0
  do 
     if (input_filename(len+1) == C_NULL_CHAR) exit
     len = len + 1
  end do
  
  print *, "Reading : ", (input_filename(i), i=1,len)

  sizes%num_primary = 3
  sizes%num_kinetic_minerals = 1
  sizes%num_aqueous_complexes = 6
  sizes%num_surface_sites = 0
  sizes%num_ion_exchange_sites = 0
  call PFloTranAlquimia_PrintSizes(sizes)

end subroutine PFloTranAlquimia_Setup


! **************************************************************************** !
subroutine PFloTranAlquimia_ProcessCondition(condition, sizes, state) bind(C)

  use c_interface_module

#include "alquimia_containers.h90"

  type (alquimia_condition_f), intent(in) :: condition
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_state_f), intent(inout) :: state
  type (alquimia_constraint_f), pointer :: local_constraints(:)
  type (alquimia_constraint_f), pointer :: constraint
  character (ALQUIMIA_MAX_STRING_LENGTH) :: name
  character (ALQUIMIA_MAX_STRING_LENGTH) :: constraint_type
  character (ALQUIMIA_MAX_STRING_LENGTH) :: associated_species
  real (c_double) :: constraint_value
  integer :: i

  print *, "Fortran process condition : "
  call  c_f_string(condition%name, name)
  print *, "  condition name : ", trim(name)
  write (*, '(a i3)') "     num constraints : ", condition%num_constraints
  call c_f_pointer(condition%constraints, local_constraints, (/condition%num_constraints/))
  do i = 1, condition%num_constraints
     !call c_f_pointer(local_constraints(i), constraint)
     call c_f_string(local_constraints(i)%primary_species, name)
     call c_f_string(local_constraints(i)%constraint_type, constraint_type)
     call c_f_string(local_constraints(i)%associated_species, associated_species)
     constraint_value = local_constraints(i)%value
     write (*, '(a a a)', advance='no') "        ", trim(name), " : "
     write (*, '(a a a)', advance='no') trim(constraint_type), " ", trim(associated_species)
     write (*, '(a f6.2)') " ", constraint_value
  end do
  
end subroutine PFloTranAlquimia_ProcessCondition


! **************************************************************************** !
subroutine PFloTranAlquimia_ReactionStepOperatorSplit() bind(C)

#include "alquimia_containers.h90"

  print *, "Fortran reaction step operator split."
end subroutine PFloTranAlquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine PFloTranAlquimia_GetAuxiliaryOutput() bind(C)

#include "alquimia_containers.h90"

  print *, "Fortran get auxiliary output."
end subroutine PFloTranAlquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine PFloTranAlquimia_GetEngineMetaData(sizes, metadata) bind(C)

#include "alquimia_containers.h90"

  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_metadata_f), intent(out) :: metadata
  integer (c_int), pointer :: primary_indices(:)
  type (c_ptr), pointer :: names(:)
  character (kind=c_char), pointer :: name
  character (len=256) :: f_name
  integer :: i

  print *, "Fortran get engine metadata."

  metadata%thread_safe = .true.
  metadata%temperature_dependent = .false.
  metadata%pressure_dependent = .false.
  metadata%porosity_update = .true.
  metadata%operator_splitting = .true.
  metadata%global_implicit = .false.
  metadata%index_base = 1

  ! associate indices with the C memory
  if (c_associated(metadata%primary_indices)) then
     call c_f_pointer(metadata%primary_indices, primary_indices, (/sizes%num_primary/))
  else
     ! error...
     print *, "c primary indicies not associated"
  end if
  if (.not. associated(primary_indices)) then
     ! error
     print *, "f primary indices not associated"
  end if

  primary_indices(1) = 1
  primary_indices(2) = 2
  primary_indices(3) = 3

!
! NOTE(bja): I can't figure out how to get arrays of strings passed
! back and forth between C and fortran. Actually I don't understand
! arrays of strings in fortran. For now we'll just require C to loop
! through each string and call GetPrimaryNameFromIndex...
!

!!$  ! associate 'names' with the C array of pointers to pointers
!!$  if (c_associated(metadata%primary_names)) then
!!$     call c_f_pointer(metadata%primary_names, names, (/sizes%num_primary/))
!!$  else
!!$     print *, "c names not associated"
!!$  endif
!!$
!!$  if (.not. associated(names)) then
!!$     print *, "f names not associated"
!!$  end if
!!$
!!$  do i=1,sizes%num_primary
!!$     ! now associate names(i) pointer to an array of chars with a string?
!!$     if (c_associated(names(i))) then
!!$        call c_f_pointer(names(i), name, (/ALQUIMIA_MAX_STRING_LENGTH/))
!!$     else
!!$        print *, "c name not associated"        
!!$     end if
!!$     if (.not. associated(name)) then
!!$        print *, "f name not associated"
!!$     end if
!!$     write(f_name, "(a,i1,a)"),"species_", i, C_NULL_CHAR
!!$     name = trim(f_name)
!!$     print *, trim(f_name), "==", name(:)
!!$  end do

  call PFloTranAlquimia_PrintMetaData(sizes, metadata)

end subroutine PFloTranAlquimia_GetEngineMetaData

subroutine PFloTranAlquimia_GetPrimaryNameFromIndex(primary_index, primary_name) bind(C)
#include "alquimia_containers.h90"
  integer (c_int), intent(in) :: primary_index
  character(kind=c_char), dimension(*), intent(out) :: primary_name 
  character (len=5) :: names(3) = (/character(len=5) :: "H+   ","Ca++ ", "HCO3-"/)
  integer(4) :: i
  character (len=5) :: name
  print *, "Fortran GetPrimaryNameFromIndex() :"
  print *, "primary index = ", primary_index
  ! NOTE(bja): should be possible to do this w/o character by character copying?
  name = names(primary_index)
  do i = 1, 5
     primary_name(i:i) = name(i:i)
  end do
  primary_name(6:6) = C_NULL_CHAR
end subroutine PFloTranAlquimia_GetPrimaryNameFromIndex

! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PFloTranAlquimia_PrintSizes(sizes)

#include "alquimia_containers.h90"

  type (alquimia_sizes_f), intent(in) :: sizes
  print *, "size : "
  print *, "  num primary : ", sizes%num_primary
  print *, "  num kinetics minerals : ", sizes%num_kinetic_minerals
  print *, "  num aqueous complexes : ", sizes%num_aqueous_complexes
  print *, "  num surface sites : ", sizes%num_surface_sites
  print *, "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PFloTranAlquimia_PrintSizes


! **************************************************************************** !
subroutine PFloTranAlquimia_PrintMetadata(sizes, metadata)

#include "alquimia_containers.h90"

  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_metadata_f), intent(in) :: metadata
  integer (c_int), pointer :: indices(:)
  type (c_ptr), pointer :: names(:)
  character (c_char), pointer :: name
  integer (c_int) :: i
  print *, "metadata : "
  print *, "  thread safe : ", metadata%thread_safe
  print *, "  temperature dependent : ", metadata%temperature_dependent
  print *, "  pressure dependent : ", metadata%pressure_dependent
  print *, "  porosity update : ", metadata%porosity_update
  print *, "  index base : ", metadata%index_base
  print *, "  primary indices : "
  call c_f_pointer(metadata%primary_indices, indices, (/sizes%num_primary/))
  do i=1, sizes%num_primary
     print *, indices(i)
  end do
  call c_f_pointer(metadata%primary_names, names, (/sizes%num_primary/))
  do i=1, sizes%num_primary
     call c_f_pointer(names(i), name, 256)
     print *, names
  end do
end subroutine PFloTranAlquimia_PrintMetadata


! **************************************************************************** !
subroutine PFloTranAlquimia_PrintStatus(status) bind(C)

#include "alquimia_containers.h90"

  type (alquimia_engine_status_f), intent(in) :: status
  print *, "status : "
  print *, "  num rhs evaluation  : ", status%num_rhs_evaluations
  print *, "  num jacobian evaluations : ", status%num_jacobian_evaluations
  print *, "  num newton iterations : ", status%num_newton_iterations
  print *, "  converged  : ", status%converged

end subroutine PFloTranAlquimia_PrintStatus

! **************************************************************************** !
!
! Alquimia inteface
!
! Notes:
!
!  * (bja) 2012-12 - different fortran compilers use different name
!    mangling conventions for fortran modules:
!
!    gfortran : ___modulename_module_MOD_procedurename_
!
!    intel : _modulename_mp_procedurename_
!
!    as a consequence we can't put the alquimia interface into a
!    module unless we muddy the calling side using the preprocessor to
!    contol the function names. It doesn't seem worth it at this
!    time...
!
! **************************************************************************** !

module LocalState
  ! This is the data structure that stores the persistent data for
  ! pflotran, (e.g. reaction network).
  !
  ! We pass it back and forth with c as a void pointer so we don't
  ! have to a global variable.
  !
  ! It is NOT part of the alquimia interface, and the driver code
  ! should not use or depend on it!


  use Constraint_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Input_module

  type :: pflotran_internal_state
     ! NOTE(bja): these are fortran pointers, so this struct can not
     ! be unpacked on the c side!
     type (reaction_type), pointer :: reaction
     type (input_type), pointer :: input
     type (option_type), pointer :: option
     type (global_auxvar_type), pointer :: global_auxvars
     type (reactive_transport_auxvar_type), pointer :: rt_auxvars
     type(tran_constraint_list_type), pointer :: transport_constraints
     type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  end type pflotran_internal_state
end module LocalState

! **************************************************************************** !
subroutine PFloTran_Alquimia_Setup(pft_internal_state, input_filename, sizes) bind(C)

  use c_interface_module

  ! pflotran
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Database_module
  use Option_module
  use Input_module

  use LocalState
  use BatchChem

#include "alquimia_containers.h90"

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (alquimia_sizes_f), intent(inout) :: sizes
  type (c_ptr), intent(inout) :: pft_internal_state

  type(pflotran_internal_state), pointer :: internal_state

  ! local variables
  PetscErrorCode :: ierr
  PetscBool :: option_found
  character(len=ALQUIMIA_MAX_STRING_LENGTH) :: string
  character(len=ALQUIMIA_MAX_STRING_LENGTH) :: filename_out
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars

  integer :: len
  integer :: i

  print *, "PFloTran_Alquimia_Setup() : "

  ! setup pflotran's option object, including mpi
  option => OptionCreate()
  option%fid_out = OUT_UNIT

  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, option%global_commsize, ierr)
  call MPI_Comm_group(MPI_COMM_WORLD, option%global_group, ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group

  ! set the pflotran input file name
  call c_f_string(input_filename,   option%input_filename)
  print *, "  Reading : ", trim(option%input_filename)

  input => InputCreate(IN_UNIT, option%input_filename, option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  !
  ! manual initialization...
  !
  option%nphase = 1


  !
  ! initialize chemistry
  !
  call BatchChemInitializeReactions(option, input, reaction)

  !
  ! create the storage containers
  !
  ! NOTE(bja) : batch chem --> one cell

  ! global_auxvars --> cell by cell temperature, pressure, saturation, density
  allocate(global_auxvars)
  call GlobalAuxVarInit(global_auxvars, option)

  ! rt_auxvars --> cell by cell chemistry data
  allocate(rt_auxvars)
  call RTAuxVarInit(rt_auxvars, reaction, option)

  ! assign default state values, not really needed?
  global_auxvars%pres = option%reference_pressure
  global_auxvars%temp = option%reference_temperature
  ! global_auxvars%den_kg = option%reference_water_density
  ! NOTE(bja): option%ref_density = 0.0, so we set it manually.
  global_auxvars%den_kg = option%reference_water_density
  global_auxvars%sat = option%reference_saturation  

  !
  ! Grab sizes info that the driver needs to finish initializing
  !
  sizes%num_primary = reaction%ncomp
  sizes%num_kinetic_minerals = reaction%mineral%nkinmnrl
  sizes%num_aqueous_complexes = reaction%neqcplx
  sizes%num_surface_sites = -1
  sizes%num_ion_exchange_sites = reaction%neqionxrxn

  !call PFloTran_Alquimia_PrintSizes(sizes)

  !
  ! save pflotran's persistent data to a struct so the driver can
  ! store it for us
  !
  allocate(internal_state)
  internal_state%reaction => reaction
  internal_state%option => option
  internal_state%input => input
  internal_state%global_auxvars => global_auxvars
  internal_state%rt_auxvars => rt_auxvars
  if (.not. c_associated(pft_internal_state)) then
     print *, "pft internal state is null"
  endif
  pft_internal_state = c_loc(internal_state)

end subroutine PFloTran_Alquimia_Setup


! **************************************************************************** !
subroutine PFloTran_Alquimia_Shutdown(pft_internal_state) bind(c)

  use c_interface_module
  use LocalState

#include "alquimia_containers.h90"

  type (c_ptr), intent(inout) :: pft_internal_state
  type(pflotran_internal_state), pointer :: internal_state

  print *, "PFloTran_Alquimia_Shutdown() : "

  call c_f_pointer(pft_internal_state, internal_state)

  call TranConstraintCouplerDestroy(internal_state%constraint_coupler)
  call TranConstraintDestroyList(internal_state%transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(internal_state%rt_auxvars)
  !call GlobalAuxVarDestroy(internal_state%global_auxvars)
  call ReactionDestroy(internal_state%reaction)
  call InputDestroy(internal_state%input)
  call OptionDestroy(internal_state%option)

end subroutine PFloTran_Alquimia_Shutdown


! **************************************************************************** !
subroutine PFloTran_Alquimia_ProcessCondition(pft_internal_state, condition, &
     sizes, state) bind(C)

  use c_interface_module
  use LocalState

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
  PetscReal :: porosity
  integer :: i

  type (c_ptr), intent(inout) :: pft_internal_state

  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  write (*, '(a)') "PFloTran_Alquimia_ProcessCondition() : "

  ! copy the state info into pflotran's auxvars
  internal_state%global_auxvars%den_kg = state%density_water
  internal_state%global_auxvars%sat = state%saturation
  porosity = state%porosity
  internal_state%global_auxvars%temp = state%temperature
  internal_state%global_auxvars%pres = state%aqueous_pressure

  ! process the condition
  call  c_f_string(condition%name, name)
  write (*, '(a a)') "processing : ", trim(name)
  if (condition%num_constraints > 0) then
     ! the driver is supplying the constraint data
     write (*, '(a)') "NOTE: driver supplied conditions are not implemeted."
     write (*, '(a a)') "  echoing condition: ", trim(name)
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
  else
     ! the driver just supplied a name, so we check for a constraint
     ! with that name in the pflotran input file and use that.
     write (*, '(a a)') "Looking for pflotran constraint : ", trim(name)
  end if

end subroutine PFloTran_Alquimia_ProcessCondition


! **************************************************************************** !
subroutine PFloTran_Alquimia_ReactionStepOperatorSplit(pft_internal_state) bind(C)

  use LocalState

#include "alquimia_containers.h90"

  type (c_ptr), intent(inout) :: pft_internal_state

  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  print *, "PFloTran_Alquimia_ReactionStepOperatorSplit() :"
end subroutine PFloTran_Alquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetAuxiliaryOutput(pft_internal_state) bind(C)

  use LocalState

#include "alquimia_containers.h90"

  type (c_ptr), intent(inout) :: pft_internal_state

  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  print *, "PFloTran_Alquimia_GetAuxiliaryOutput() :"
end subroutine PFloTran_Alquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine PFloTran_Alquimia_GetEngineMetaData(pft_internal_state, &
     sizes, metadata) bind(C)

  use LocalState

#include "alquimia_containers.h90"

  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_metadata_f), intent(out) :: metadata
  integer (c_int), pointer :: primary_indices(:)
  integer :: i, num_primary
  type (c_ptr), intent(inout) :: pft_internal_state
  type(pflotran_internal_state), pointer :: internal_state

  !print *, "PFloTran_Alquimia_GetEngineMetadata() :"

  call c_f_pointer(pft_internal_state, internal_state)

  num_primary = internal_state%reaction%ncomp

  ! TODO(bja) : can we extract this from pflotran without hardcoding?
  metadata%thread_safe = .true.
  metadata%temperature_dependent = .true.
  metadata%pressure_dependent = .true.
  metadata%porosity_update = internal_state%reaction%update_porosity
  metadata%operator_splitting = .true.
  metadata%global_implicit = .false.
  metadata%index_base = 1

  ! associate indices with the C memory
  if (c_associated(metadata%primary_indices)) then
     call c_f_pointer(metadata%primary_indices, primary_indices, (/num_primary/))
  else
     ! error...
     print *, "c primary indicies not associated"
  end if
  if (.not. associated(primary_indices)) then
     ! error
     print *, "f primary indices not associated"
  end if

  ! NOTE(bja) : if the order in reaction%primary_species_names always
  ! is not correct we need to modifiy this...
  do i=1, num_primary
     primary_indices(i) = i
  enddo

!
! NOTE(bja): I can't figure out how to get arrays of strings passed
! back and forth between C and fortran. Actually I don't understand
! arrays of strings in fortran. For now we'll just require C to loop
! through each string and call GetPrimaryNameFromIndex...
!

  !call PFloTran_Alquimia_PrintMetaData(sizes, metadata)

end subroutine PFloTran_Alquimia_GetEngineMetaData

! **************************************************************************** !
subroutine PFloTran_Alquimia_GetPrimaryNameFromIndex(pft_internal_state, &
  primary_index, primary_name) bind(C)

  use c_interface_module
  use LocalState

#include "alquimia_containers.h90"

  integer (c_int), intent(in) :: primary_index
  character(kind=c_char), dimension(*), intent(out) :: primary_name
  character (len=ALQUIMIA_MAX_WORD_LENGTH), pointer :: primary_list(:)
  type (c_ptr), intent(inout) :: pft_internal_state
  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  primary_list => internal_state%reaction%primary_species_names

  !print *, "PFloTran_Alquimia_GetPrimaryNameFromIndex() :"
  !print *, "primary index = ", primary_index

  call f_c_string_chars(trim(primary_list(primary_index)), &
       primary_name, ALQUIMIA_MAX_STRING_LENGTH)

end subroutine PFloTran_Alquimia_GetPrimaryNameFromIndex

! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PFloTran_Alquimia_PrintSizes(sizes)

#include "alquimia_containers.h90"

  type (alquimia_sizes_f), intent(in) :: sizes
  print *, "size : "
  print *, "  num primary : ", sizes%num_primary
  print *, "  num kinetics minerals : ", sizes%num_kinetic_minerals
  print *, "  num aqueous complexes : ", sizes%num_aqueous_complexes
  print *, "  num surface sites : ", sizes%num_surface_sites
  print *, "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PFloTran_Alquimia_PrintSizes


! **************************************************************************** !
subroutine PFloTran_Alquimia_PrintMetadata(sizes, metadata)

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
end subroutine PFloTran_Alquimia_PrintMetadata


! **************************************************************************** !
subroutine PFloTran_Alquimia_PrintStatus(status) bind(C)

#include "alquimia_containers.h90"

  type (alquimia_engine_status_f), intent(in) :: status
  print *, "status : "
  print *, "  num rhs evaluation  : ", status%num_rhs_evaluations
  print *, "  num jacobian evaluations : ", status%num_jacobian_evaluations
  print *, "  num newton iterations : ", status%num_newton_iterations
  print *, "  converged  : ", status%converged

end subroutine PFloTran_Alquimia_PrintStatus

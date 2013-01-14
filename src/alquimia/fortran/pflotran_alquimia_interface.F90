! **************************************************************************** !
!
! PFloTran Alquimia Inteface module
!
! Notes:
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
!    pflotran_alquimia.F90!
!
! **************************************************************************** !

module PFloTranAlquimiaInterface_module

  use, intrinsic :: iso_c_binding

  use Constraint_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Input_module

  implicit none

  public :: Setup, &
       Shutdown, &
       ProcessCondition, &
       ReactionStepOperatorSplit, &
       GetAuxiliaryOutput, &
       GetEngineMetaData, &
       GetPrimaryNameFromIndex

#include "alquimia_containers.h90"

  private :: InitializePFloTranReactions, &
       ReadPFloTranConstraints, &
       ProcessPFloTranConstraint, &
       CopyAuxVarsToState


  type :: pflotran_internal_state
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
     type (reaction_type), pointer :: reaction
     type (option_type), pointer :: option
     type (global_auxvar_type), pointer :: global_auxvars
     type (reactive_transport_auxvar_type), pointer :: rt_auxvars
     type(tran_constraint_list_type), pointer :: transport_constraints
     type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  end type pflotran_internal_state

contains


! **************************************************************************** !
!
! Public interface routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine Setup(pft_internal_state, input_filename, sizes)

  use c_interface_module

  ! pflotran
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Database_module
  use Option_module
  use Input_module

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
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler 

  integer :: len
  integer :: i

  write (*, '(a)') "PFloTran_Alquimia_Setup() : "

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
  call InitializePFloTranReactions(option, input, reaction)

  !
  ! create the various storage containers pflotran will need later
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

  ! NOTE(bja) : constraint coupler is only needed for printing the
  ! processed constraints to pflotran.out. But destroying the coupler
  ! destroys the auxvars as well, so we need to keep it around long
  ! term.
  allocate(constraint_coupler)
  constraint_coupler => TranConstraintCouplerCreate(option)

  !
  ! Read the constraints so we can finish using the input file.
  !
  call ReadPFloTranConstraints(option, input, reaction, &
     global_auxvars, rt_auxvars, transport_constraints)

  ! close the input file because we don't need it any more
  call InputDestroy(input)

  !
  ! save pflotran's persistent data to a struct so the driver can
  ! store it for us
  !
  allocate(internal_state)
  internal_state%reaction => reaction
  internal_state%option => option
  internal_state%global_auxvars => global_auxvars
  internal_state%rt_auxvars => rt_auxvars
  internal_state%constraint_coupler => constraint_coupler
  internal_state%transport_constraints => transport_constraints

  if (.not. c_associated(pft_internal_state)) then
     print *, "pft internal state is null"
  endif
  pft_internal_state = c_loc(internal_state)

  !
  ! Grab sizes info that the driver needs to finish initializing
  !
  sizes%num_primary = reaction%ncomp
  sizes%num_kinetic_minerals = reaction%mineral%nkinmnrl
  sizes%num_aqueous_complexes = reaction%neqcplx
  sizes%num_surface_sites = -1
  sizes%num_ion_exchange_sites = reaction%neqionxrxn

  call PrintSizes(sizes)

  write (*, '(a)') "PFloTran_Alquimia_Setup() : successful."

end subroutine Setup


! **************************************************************************** !
subroutine Shutdown(pft_internal_state)

  use c_interface_module

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
  call OptionDestroy(internal_state%option)

end subroutine Shutdown


! **************************************************************************** !
subroutine ProcessCondition(pft_internal_state, condition, &
     sizes, state)

  use c_interface_module

  use String_module
  use Constraint_module

#include "alquimia_containers.h90"

#include "definitions.h"
#include "finclude/petsclog.h"

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
  type(tran_constraint_type), pointer :: tran_constraint

  type (c_ptr), intent(inout) :: pft_internal_state
  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  write (*, '(a)') "PFloTran_Alquimia_ProcessCondition() : "

  ! copy the driver's state info into pflotran's auxvars
  internal_state%global_auxvars%den_kg = state%density_water
  internal_state%global_auxvars%sat = state%saturation
  porosity = state%porosity
  internal_state%global_auxvars%temp = state%temperature
  internal_state%global_auxvars%pres = state%aqueous_pressure

  !
  ! process the condition
  !
  call  c_f_string(condition%name, name)
  write (*, '(a a)') "processing : ", trim(name)

  if (condition%num_constraints > 0) then
     ! the driver is supplying the constraint data, so we need to
     ! construct an object in pflotran's internal format.
     tran_constraint => TranConstraintCreate(internal_state%option)
     ! TODO(bja) : manual setup goes here...
     tran_constraint%name = trim(name)
     tran_constraint%requires_equilibration = PETSC_TRUE

     call TranConstraintAddToList(tran_constraint, internal_state%transport_constraints)
     ! TODO(bja) : the rest of this if block eventually goes away....
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
     tran_constraint => internal_state%transport_constraints%first
     do
        if (associated(tran_constraint)) then
           ! check the name of this constraint
           if (StringCompareIgnoreCase(tran_constraint%name, name)) then
              ! found the constraint we are looking for, bail from the
              ! loop with the current pointer
              exit
           else
              ! check the next constraint
              tran_constraint => tran_constraint%next
           end if
        else
           ! end of the list (or empty list) without out finding a match.
           write (*, '(a a)') "Could not find pflotran constraint : ", trim(name)
           ! TODO(bja) : else report an error to the driver
           exit
        end if
     end do
  end if

  if (associated(tran_constraint)) then
     ! tran_constraint should be valid. Now we can ask pflotran to
     ! process it...
     call ProcessPFloTranConstraint( &
          internal_state%option, &
          internal_state%reaction, &
          internal_state%global_auxvars, &
          internal_state%rt_auxvars, &
          tran_constraint, &
          internal_state%constraint_coupler)
     ! now repack the processed constraint data into the alquimia
     ! struct for the driver.
     call CopyAuxVarsToState( &
          internal_state%reaction, &
          internal_state%global_auxvars, &
          internal_state%rt_auxvars, &
          state)
  end if

end subroutine ProcessCondition


! **************************************************************************** !
subroutine ReactionStepOperatorSplit(pft_internal_state)

#include "alquimia_containers.h90"

  type (c_ptr), intent(inout) :: pft_internal_state

  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  print *, "PFloTran_Alquimia_ReactionStepOperatorSplit() :"
end subroutine ReactionStepOperatorSplit


! **************************************************************************** !
subroutine GetAuxiliaryOutput(pft_internal_state)

#include "alquimia_containers.h90"

  type (c_ptr), intent(inout) :: pft_internal_state

  type(pflotran_internal_state), pointer :: internal_state

  call c_f_pointer(pft_internal_state, internal_state)

  print *, "PFloTran_Alquimia_GetAuxiliaryOutput() :"
end subroutine GetAuxiliaryOutput


! **************************************************************************** !
subroutine GetEngineMetaData(pft_internal_state, &
     sizes, metadata)

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

end subroutine GetEngineMetaData

! **************************************************************************** !
subroutine GetPrimaryNameFromIndex(pft_internal_state, &
  primary_index, primary_name)

  use c_interface_module

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

end subroutine GetPrimaryNameFromIndex


! **************************************************************************** !
!
! Private work routines
!
! **************************************************************************** !

subroutine InitializePFloTranReactions(option, input, reaction)

  use Reaction_module
  use Reaction_Aux_module
  use Database_module
  use Option_module
  use Input_module
  use String_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string

  ! check for a chemistry block in the  input file
  string = "CHEMISTRY"
  call InputFindStringInFile(input, option, string)
  if (.not.InputError(input)) then
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

end subroutine InitializePFloTranReactions

! **************************************************************************** !
!
! ReadPFloTranConstraints
!
! We are just reading the data from the input file. No processing is
! done here because we don't know yet if we are using these.
!
! **************************************************************************** !

subroutine ReadPFloTranConstraints(option, input, reaction, &
     global_auxvars, rt_auxvars, transport_constraints)

  use Reaction_module
  use Reaction_Aux_module
  use Database_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Input_module
  use Constraint_module
  use String_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  type(option_type), pointer, intent(in) :: option
  type(input_type), pointer, intent(inout) :: input
  type(reaction_type), pointer, intent(inout) :: reaction
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvars
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvars
  type(tran_constraint_list_type), pointer, intent(inout) :: transport_constraints

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_type), pointer :: tran_constraint

  allocate(transport_constraints)
  call TranConstraintInitList(transport_constraints)

  ! look through the input file
  rewind(input%fid)        
  do
    call InputReadFlotranString(input, option)
    if (InputError(input)) exit

    call InputReadWord(input, option, word, PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call printMsg(option)

    select case(trim(card))
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input, option, tran_constraint%name, PETSC_TRUE)
        call InputErrorMsg(input, option, 'constraint', 'name') 
        call printMsg(option, tran_constraint%name)
        call TranConstraintRead(tran_constraint, reaction, input, option)
        call TranConstraintAddToList(tran_constraint, transport_constraints)
        nullify(tran_constraint)

      case default
         ! do nothing
    end select
  enddo

end subroutine ReadPFloTranConstraints

subroutine ProcessPFloTranConstraint(option, reaction, &
     global_auxvars, rt_auxvars, tran_constraint, constraint_coupler)
  use Reaction_module
  use Reaction_Aux_module
  use Database_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Constraint_module
  use Option_module
  use String_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  type(option_type), pointer, intent(in) :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  PetscBool :: use_prev_soln_as_guess
  PetscInt :: num_iterations

  !
  ! process constraints
  !
  num_iterations = 0
  use_prev_soln_as_guess = PETSC_FALSE

  if (.not. associated(tran_constraint)) then
     ! TODO(bja) : report error
  end if
  ! initialize constraints
  option%io_buffer = "initializing constraint : " // tran_constraint%name
  call printMsg(option)
  call ReactionProcessConstraint(reaction, &
       tran_constraint%name, &
       tran_constraint%aqueous_species, &
       tran_constraint%minerals, &
       tran_constraint%surface_complexes, &
       tran_constraint%colloids, &
       tran_constraint%biomass, &
       option)

  ! equilibrate
  option%io_buffer = "equilibrate constraint : " // tran_constraint%name
  call printMsg(option)
  call ReactionEquilibrateConstraint(rt_auxvars, global_auxvars, reaction, &
       tran_constraint%name, &
       tran_constraint%aqueous_species, &
       tran_constraint%minerals, &
       tran_constraint%surface_complexes, &
       tran_constraint%colloids, &
       tran_constraint%biomass, &
       option%reference_porosity, &
       num_iterations, &
       use_prev_soln_as_guess, &
       option)

  ! link the constraint to the constraint coupler so we can print it
  constraint_coupler%constraint_name = tran_constraint%name
  constraint_coupler%aqueous_species => tran_constraint%aqueous_species
  constraint_coupler%minerals => tran_constraint%minerals
  constraint_coupler%surface_complexes => tran_constraint%surface_complexes
  constraint_coupler%colloids => tran_constraint%colloids
  constraint_coupler%global_auxvar => global_auxvars
  constraint_coupler%rt_auxvar => rt_auxvars     
  call ReactionPrintConstraint(constraint_coupler, reaction, option)


end subroutine ProcessPFloTranConstraint


! **************************************************************************** !
subroutine CopyAuxVarsToState(reaction, &
     global_auxvars, rt_auxvars, state)

  use Reaction_aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

#include "alquimia_containers.h90"

  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  type(reaction_type), pointer :: reaction
  type (alquimia_state_f), intent(inout) :: state
  real (c_double), pointer :: larray(:)
  integer :: i, phase_index

  phase_index = 1

  call c_f_pointer(state%total_primary, larray, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
        larray(i) = rt_auxvars%total(i, phase_index)
  end do

  call c_f_pointer(state%free_ion, larray, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
        larray(i) = rt_auxvars%pri_molal(i)
  end do

  call c_f_pointer(state%total_sorbed, larray, (/reaction%neqsorb/))
  do i = 1, reaction%neqsorb
        larray(i) = rt_auxvars%total_sorb_eq(i)
  end do

end subroutine CopyAuxVarsToState

! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PrintSizes(sizes)

#include "alquimia_containers.h90"

  type (alquimia_sizes_f), intent(in) :: sizes
  print *, "size : "
  print *, "  num primary : ", sizes%num_primary
  print *, "  num kinetics minerals : ", sizes%num_kinetic_minerals
  print *, "  num aqueous complexes : ", sizes%num_aqueous_complexes
  print *, "  num surface sites : ", sizes%num_surface_sites
  print *, "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PrintSizes


! **************************************************************************** !
subroutine PrintMetadata(sizes, metadata)

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
end subroutine PrintMetadata


! **************************************************************************** !
subroutine PrintStatus(status)

#include "alquimia_containers.h90"

  type (alquimia_engine_status_f), intent(in) :: status
  print *, "status : "
  print *, "  num rhs evaluation  : ", status%num_rhs_evaluations
  print *, "  num jacobian evaluations : ", status%num_jacobian_evaluations
  print *, "  num newton iterations : ", status%num_newton_iterations
  print *, "  converged  : ", status%converged

end subroutine PrintStatus

end module PFloTranAlquimiaInterface_module

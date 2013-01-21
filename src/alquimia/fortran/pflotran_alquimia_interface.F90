! **************************************************************************** !
!
! PFloTran Alquimia Inteface module
!
! Author: Benjamin Andre
!
! Notes:
!
!  * Public function call signatures, including intent, are dictated
!    by the alquimia API.
!
!  * alquimia data structures defined in alquimia_containers.h90 are
!    dictated by the alquimia API.
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
!    pflotran_alquimia.F90!
!
! **************************************************************************** !

module PFloTranAlquimiaInterface_module

  ! pflotran modules
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Constraint_module

  implicit none

  public :: Setup, &
       Shutdown, &
       ProcessCondition, &
       ReactionStepOperatorSplit, &
       GetAuxiliaryOutput, &
       GetEngineMetaData, &
       GetPrimaryNameFromIndex

  private :: InitializePFloTranReactions, &
       ReadPFloTranConstraints, &
       ProcessPFloTranConstraint, &
       CopyAlquimiaToAuxVars, &
       CopyAuxVarsToAlquimia

  integer(kind=8), private, parameter :: integrity_check_value = &
       b"0100001001100101011011100100000101101110011001000111001001100101"

  type, private :: pflotran_engine_state
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
     type (option_type), pointer :: option
     type (reaction_type), pointer :: reaction
     type (reactive_transport_auxvar_type), pointer :: rt_auxvars
     type (global_auxvar_type), pointer :: global_auxvars
     type(tran_constraint_list_type), pointer :: transport_constraints
     type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  end type pflotran_engine_state

contains


! **************************************************************************** !
!
! Public interface routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine Setup(input_filename, pft_engine_state, sizes, status)
!  NOTE: Function signature is dictated by the alquimia API.
!
!  NOTE: Assumes that MPI_Init() and / or PetscInitialize() have already
!    been called by the driver

  use, intrinsic :: iso_c_binding

  use c_interface_module

  ! pflotran
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Input_module

  implicit none

#include "alquimia_containers.h90"

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (c_ptr), intent(out) :: pft_engine_state
  type (alquimia_sizes_f), intent(out) :: sizes
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  type(pflotran_engine_state), pointer :: engine_state
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
  write (*, '(a, a)') "  Reading : ", trim(option%input_filename)

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
  option%liquid_phase = 1
  option%reference_water_density = 998.2


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
  global_auxvars%den_kg = option%reference_water_density
  global_auxvars%sat = option%reference_saturation  

  ! NOTE(bja) : constraint coupler is only needed for printing the
  ! processed constraints to pflotran.out. But destroying the coupler
  ! destroys the auxvars as well, so we need to keep it around long
  ! term.
  constraint_coupler => TranConstraintCouplerCreate(option)

  !
  ! Read the constraints so we can finish using the input file.
  !
  call ReadPFloTranConstraints(option, input, reaction, transport_constraints)

  ! close the input file because we don't need it any more
  call InputDestroy(input)

  !
  ! save pflotran's persistent data to a struct so the driver can
  ! store it for us
  !
  allocate(engine_state)
  engine_state%integrity_check = integrity_check_value
  engine_state%option => option
  engine_state%reaction => reaction
  engine_state%rt_auxvars => rt_auxvars
  engine_state%global_auxvars => global_auxvars
  engine_state%constraint_coupler => constraint_coupler
  engine_state%transport_constraints => transport_constraints

  pft_engine_state = c_loc(engine_state)

  !
  ! Grab sizes info that the driver needs to finish initializing
  !
  sizes%num_primary = reaction%ncomp
  sizes%num_kinetic_minerals = reaction%mineral%nkinmnrl
  sizes%num_aqueous_complexes = reaction%neqcplx
  sizes%num_surface_sites = reaction%surface_complexation%nsrfcplxrxn
  sizes%num_ion_exchange_sites = reaction%neqionxrxn

  !call PrintSizes(sizes)

  status%error = ALQUIMIA_NO_ERROR
  call f_c_string_ptr("Alquimia::PFloTran::Setup() : successful.", &
       status%message, ALQUIMIA_MAX_STRING_LENGTH)
  
end subroutine Setup


! **************************************************************************** !
subroutine Shutdown(pft_engine_state, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  type(pflotran_engine_state), pointer :: engine_state

  !write (*, '(a)') "PFloTranAlquimiaInterface::Shutdown() : "

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = ALQUIMIA_ERROR_ENGINE_INTEGRITY
     call f_c_string_ptr("ERROR: pointer to engine state is not valid!", &
          status%message, ALQUIMIA_MAX_STRING_LENGTH)
     return
  end if

  call TranConstraintCouplerDestroy(engine_state%constraint_coupler)
  call TranConstraintDestroyList(engine_state%transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(engine_state%rt_auxvars)
  !call GlobalAuxVarDestroy(engine_state%global_auxvars)
  call ReactionDestroy(engine_state%reaction)
  call OptionDestroy(engine_state%option)

  deallocate(engine_state)
  nullify(engine_state)

  status%error = ALQUIMIA_NO_ERROR

end subroutine Shutdown


! **************************************************************************** !
subroutine ProcessCondition(pft_engine_state, condition, material_properties, &
     state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use String_module
  use Constraint_module

  implicit none

#include "alquimia_containers.h90"

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (alquimia_condition_f), intent(in) :: condition
  type (alquimia_material_properties_f), intent(in) :: material_properties
  type (alquimia_state_f), intent(inout) :: state
  type (alquimia_auxiliary_data_f), intent (inout) :: aux_data
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  type (alquimia_constraint_f), pointer :: local_constraints(:)
  type (alquimia_constraint_f), pointer :: constraint
  character (ALQUIMIA_MAX_STRING_LENGTH) :: name
  character (ALQUIMIA_MAX_STRING_LENGTH) :: constraint_type
  character (ALQUIMIA_MAX_STRING_LENGTH) :: associated_species
  real (c_double) :: constraint_value
  PetscReal :: porosity, volume
  integer :: i
  type(tran_constraint_type), pointer :: tran_constraint
  type(pflotran_engine_state), pointer :: engine_state

  !write (*, '(a)') "PFloTranAlquimiaInterface::ProcessCondition() : "

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = ALQUIMIA_ERROR_ENGINE_INTEGRITY
     call f_c_string_ptr("ERROR: pointer to engine state is not valid!", &
          status%message, ALQUIMIA_MAX_STRING_LENGTH)
     return
  end if

  ! NOTE(bja): do NOT call CopyAlquimiaToAuxVars() here because most
  ! of that data is uninitialized in alquimia!
  engine_state%global_auxvars%den_kg = state%density_water
  engine_state%global_auxvars%sat = state%saturation
  engine_state%global_auxvars%temp = state%temperature
  engine_state%global_auxvars%pres = state%aqueous_pressure
  porosity = state%porosity

  !
  ! process the condition
  !
  call  c_f_string(condition%name, name)
  write (*, '(a, a)') "processing : ", trim(name)

  if (condition%num_constraints > 0) then
     ! the driver is supplying the constraint data, so we need to
     ! construct an object in pflotran's internal format.
     tran_constraint => TranConstraintCreate(engine_state%option)
     tran_constraint%name = trim(name)
     tran_constraint%requires_equilibration = PETSC_TRUE
     ! TODO(bja) : remaining manual setup goes here...

     call TranConstraintAddToList(tran_constraint, engine_state%transport_constraints)
     ! TODO(bja) : the rest of this if block eventually goes away,
     ! just keeping it around for now as a reminder about how to
     ! access the constraint info....
     write (*, '(a)') "NOTE: driver supplied conditions are not implemeted."
     write (*, '(a, a)') "  echoing condition: ", trim(name)
     write (*, '(a, i3)') "     num constraints : ", condition%num_constraints
     call c_f_pointer(condition%constraints, local_constraints, (/condition%num_constraints/))
     do i = 1, condition%num_constraints
        !call c_f_pointer(local_constraints(i), constraint)
        call c_f_string(local_constraints(i)%primary_species, name)
        call c_f_string(local_constraints(i)%constraint_type, constraint_type)
        call c_f_string(local_constraints(i)%associated_species, associated_species)
        constraint_value = local_constraints(i)%value
        write (*, '(a, a, a)', advance='no') "        ", trim(name), " : "
        write (*, '(a, a, a)', advance='no') trim(constraint_type), " ", trim(associated_species)
        write (*, '(a, f6.2)') " ", constraint_value
     end do
  else
     ! the driver just supplied a name, so we check for a constraint
     ! with that name in the pflotran input file and use that.
     write (*, '(a, a)') "Looking for pflotran constraint : ", trim(name)
     tran_constraint => engine_state%transport_constraints%first
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
           write (*, '(a, a)') "Could not find pflotran constraint : ", trim(name)
           ! TODO(bja) : else report an error to the driver
           exit
        end if
     end do
  end if

  if (associated(tran_constraint)) then
     ! tran_constraint should be valid. Now we can ask pflotran to
     ! process it...
     call ProcessPFloTranConstraint( &
          engine_state%option, &
          engine_state%reaction, &
          engine_state%global_auxvars, &
          engine_state%rt_auxvars, &
          tran_constraint, &
          engine_state%constraint_coupler)
     ! now repack the processed constraint data into the alquimia
     ! struct for the driver.
     call CopyAuxVarsToAlquimia( &
          engine_state%reaction, &
          engine_state%global_auxvars, &
          engine_state%rt_auxvars, &
          state, aux_data)
  end if
  status%error = ALQUIMIA_NO_ERROR
end subroutine ProcessCondition


! **************************************************************************** !
subroutine ReactionStepOperatorSplit(pft_engine_state, &
     delta_t, material_properties, state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  real (c_double), intent(in) :: delta_t
  type (alquimia_material_properties_f), intent(in) :: material_properties
  type (alquimia_state_f), intent(inout) :: state
  type (alquimia_auxiliary_data_f), intent(inout) :: aux_data
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  type(pflotran_engine_state), pointer :: engine_state
  PetscReal :: porosity, volume

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = ALQUIMIA_ERROR_ENGINE_INTEGRITY
     call f_c_string_ptr("ERROR: pointer to engine state is not valid!", &
          status%message, ALQUIMIA_MAX_STRING_LENGTH)
     return
  end if

  !write (*, '(a)') "F_PFloTranAlquimiaInterface::ReactionStepOperatorSplit() :"

  !call PrintState(engine_state%reaction, state)

  call CopyAlquimiaToAuxVars(state, aux_data, material_properties, &
       engine_state%reaction, engine_state%global_auxvars, engine_state%rt_auxvars, &
       porosity, volume)

  call CopyAuxVarsToAlquimia( &
       engine_state%reaction, &
       engine_state%global_auxvars, &
       engine_state%rt_auxvars, &
       state, aux_data)

  status%error = ALQUIMIA_NO_ERROR

end subroutine ReactionStepOperatorSplit



! **************************************************************************** !
subroutine GetAuxiliaryOutput(pft_engine_state, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  type(pflotran_engine_state), pointer :: engine_state

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = ALQUIMIA_ERROR_ENGINE_INTEGRITY
     call f_c_string_ptr("ERROR: pointer to engine state is not valid!", &
          status%message, ALQUIMIA_MAX_STRING_LENGTH)
     return
  end if

  write (*, '(a)') "PFloTran_Alquimia_GetAuxiliaryOutput() :"
  status%error = ALQUIMIA_NO_ERROR
end subroutine GetAuxiliaryOutput


! **************************************************************************** !
subroutine GetEngineMetaData(pft_engine_state, sizes, meta_data, status)
  !  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_meta_data_f), intent(out) :: meta_data
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  integer (c_int), pointer :: primary_indices(:)
  integer :: i, num_primary
  type(pflotran_engine_state), pointer :: engine_state

  !write (*, '(a)') "PFloTran_Alquimia_GetEngineMetaData() :"

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = ALQUIMIA_ERROR_ENGINE_INTEGRITY
     call f_c_string_ptr("ERROR: pointer to engine state is not valid!", &
          status%message, ALQUIMIA_MAX_STRING_LENGTH)
     return
  end if

  num_primary = engine_state%reaction%ncomp

  ! TODO(bja) : can we extract this from pflotran without hardcoding?
  meta_data%thread_safe = .true.
  meta_data%temperature_dependent = .true.
  meta_data%pressure_dependent = .true.
  meta_data%porosity_update = engine_state%reaction%update_porosity
  meta_data%operator_splitting = .true.
  meta_data%global_implicit = .false.
  meta_data%index_base = 1

  ! associate indices with the C memory
  if (c_associated(meta_data%primary_indices)) then
     call c_f_pointer(meta_data%primary_indices, primary_indices, (/num_primary/))
  else
     ! error...
     write (*, '(a)') "c primary indicies not associated"
  end if
  if (.not. associated(primary_indices)) then
     ! error
     write (*, '(a)') "f primary indices not associated"
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

  ! NOTE(bja) : meta_data%primary_names() is empty for this call!
  !call PrintMetaData(sizes, meta_data)
  status%error = 0
end subroutine GetEngineMetaData


! **************************************************************************** !
subroutine GetPrimaryNameFromIndex(pft_engine_state, &
     primary_index, primary_name, status)
!  NOTE: not officially part of the alquimia API. Eventually this should
!    go away once we have passing arrays of strings from C to fortran
!    and back

  use, intrinsic :: iso_c_binding

  use c_interface_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  integer (c_int), intent(in) :: primary_index
  character(kind=c_char), dimension(*), intent(out) :: primary_name
  type (alquimia_engine_status_f), intent(out) :: status

  ! local variables
  character (len=ALQUIMIA_MAX_WORD_LENGTH), pointer :: primary_list(:)
  type(pflotran_engine_state), pointer :: engine_state

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = ALQUIMIA_ERROR_ENGINE_INTEGRITY
     call f_c_string_ptr("ERROR: pointer to engine state is not valid!", &
          status%message, ALQUIMIA_MAX_STRING_LENGTH)
     return
  end if

  primary_list => engine_state%reaction%primary_species_names

  !write (*, '(a)') "PFloTran_Alquimia_GetPrimaryNameFromIndex() :"
  !write (*, '(a)') "primary index = ", primary_index

  call f_c_string_chars(trim(primary_list(primary_index)), &
       primary_name, ALQUIMIA_MAX_STRING_LENGTH)

  status%error = ALQUIMIA_NO_ERROR

end subroutine GetPrimaryNameFromIndex


! **************************************************************************** !
!
! Private work routines
!
! **************************************************************************** !

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

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  type(input_type), pointer, intent(in) :: input
  type(reaction_type), pointer, intent(out) :: reaction

  ! local variables
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
subroutine ReadPFloTranConstraints(option, input, reaction, transport_constraints)
!  NOTE: We are just reading the data from the input file. No processing
!    is done here because we don't know yet if we are using these.

  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Input_module
  use Constraint_module
  use String_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  type(input_type), pointer, intent(inout) :: input
  type(reaction_type), pointer, intent(inout) :: reaction
  type(tran_constraint_list_type), pointer, intent(inout) :: transport_constraints

  ! local variables
  PetscBool :: debug = PETSC_FALSE
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

    if (debug) then
       option%io_buffer = 'pflotran card:: ' // trim(card)
       call printMsg(option)
    end if

    select case(trim(card))
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input, option, tran_constraint%name, PETSC_TRUE)
        call InputErrorMsg(input, option, 'constraint', 'name') 
        if (debug) then
           call printMsg(option, tran_constraint%name)
        end if
        call TranConstraintRead(tran_constraint, reaction, input, option)
        call TranConstraintAddToList(tran_constraint, transport_constraints)
        nullify(tran_constraint)

      case default
         ! do nothing
    end select
  enddo

end subroutine ReadPFloTranConstraints

! **************************************************************************** !
subroutine ProcessPFloTranConstraint(option, reaction, &
     global_auxvars, rt_auxvars, tran_constraint, constraint_coupler)
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Constraint_module
  use Option_module
  use String_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  type(option_type), pointer, intent(in) :: option
  type(reaction_type), pointer, intent(inout) :: reaction
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvars
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvars
  type(tran_constraint_type), pointer, intent(inout) :: tran_constraint
  type(tran_constraint_coupler_type), pointer, intent(inout) :: constraint_coupler 

  ! local variables
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
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
subroutine CopyAlquimiaToAuxVars(state, aux_data, material_prop, &
  reaction, global_auxvars, rt_auxvars, porosity, volume)

  use, intrinsic :: iso_c_binding

  use Reaction_aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (alquimia_state_f), intent(in) :: state
  type (alquimia_auxiliary_data_f), intent(in) :: aux_data
  type (alquimia_material_properties_f), intent(in) :: material_prop
  type(reaction_type), pointer, intent(inout) :: reaction
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvars
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvars
  PetscReal, intent(out) :: porosity
  PetscReal, intent(out) :: volume

  ! local variables
  real (c_double), pointer :: local_array(:)
  integer :: i, phase_index

  !write (*, '(a)') "PFloTran_Alquimia_CopyAlquimiaToAuxVars() :"

  phase_index = 1

  !
  ! state
  !
  global_auxvars%den_kg = state%density_water
  global_auxvars%sat = state%saturation
  global_auxvars%temp = state%temperature
  global_auxvars%pres = state%aqueous_pressure

  porosity = state%porosity
  volume = material_prop%volume

  !
  ! primary aqueous
  !
  call c_f_pointer(state%total_primary, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvars%total(i, phase_index) = local_array(i)
  end do

  call c_f_pointer(state%free_ion, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvars%pri_molal(i) = local_array(i)
  end do

  call c_f_pointer(aux_data%primary_activity_coeff, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvars%pri_act_coef(i) = local_array(i)
  end do

  !
  ! aqueous complexes
  !
  call c_f_pointer(aux_data%secondary_activity_coeff, local_array, &
       (/reaction%neqcplx/))
  do i = 1, reaction%neqcplx
     rt_auxvars%sec_act_coef(i) = local_array(i)
  end do

  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     rt_auxvars%mnrl_volfrac(i) = local_array(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     rt_auxvars%mnrl_area(i) = local_array(i)
  end do

  ! isotherms

  ! ion exchange

  ! surface complexation


end subroutine CopyAlquimiaToAuxVars

! **************************************************************************** !
subroutine CopyAuxVarsToAlquimia(reaction, global_auxvars, rt_auxvars, &
     state, aux_data)

  use, intrinsic :: iso_c_binding

  use Reaction_aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type(reaction_type), pointer, intent(in) :: reaction
  type(global_auxvar_type), pointer, intent(in) :: global_auxvars
  type(reactive_transport_auxvar_type), pointer, intent(in) :: rt_auxvars
  type (alquimia_state_f), intent(inout) :: state
  type (alquimia_auxiliary_data_f), intent(inout) :: aux_data

  ! local variables
  real (c_double), pointer :: local_array(:)
  integer :: i, phase_index

  phase_index = 1

  !write (*, '(a)') "PFloTran_Alquimia_CopyAuxVarsToAlquimia() :"

  ! TODO(bja) : state... if pressure/temp/porosity updates are allowed...?

  !
  ! primary aqueous species
  !
  call c_f_pointer(state%total_primary, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     local_array(i) = rt_auxvars%total(i, phase_index)
  end do

  call c_f_pointer(state%free_ion, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     local_array(i) = rt_auxvars%pri_molal(i)
  end do

  call c_f_pointer(aux_data%primary_activity_coeff, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     local_array(i) = rt_auxvars%pri_act_coef(i)
  end do

  !
  ! secondary aqueous complexes
  !
  call c_f_pointer(aux_data%secondary_activity_coeff, local_array, &
       (/reaction%neqcplx/))
  do i = 1, reaction%neqcplx
     local_array(i) = rt_auxvars%sec_act_coef(i)
  end do

  !
  ! sorbed
  !
  call c_f_pointer(state%total_sorbed, local_array, (/reaction%neqsorb/))
  do i = 1, reaction%neqsorb
     local_array(i) = rt_auxvars%total_sorb_eq(i)
  end do

  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     local_array(i) = rt_auxvars%mnrl_volfrac(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     local_array(i) = rt_auxvars%mnrl_area(i)
  end do

  !
  ! ion exchange
  !
  call c_f_pointer(aux_data%ion_exchange_ref_cation_conc, local_array, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     local_array(i) = rt_auxvars%eqionx_ref_cation_sorbed_conc(i)
  end do


  !
  ! equilibrium surface complexation
  !
  call c_f_pointer(aux_data%surface_complex_free_site_conc, local_array, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     local_array(i) = rt_auxvars%srfcplxrxn_free_site_conc(i)
  end do


end subroutine CopyAuxVarsToAlquimia

! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PrintSizes(sizes)

  use, intrinsic :: iso_c_binding

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (alquimia_sizes_f), intent(in) :: sizes

  write (*, '(a)') "size : "
  write (*, '(a, i4)') "  num primary : ", sizes%num_primary
  write (*, '(a, i4)') "  num kinetics minerals : ", sizes%num_kinetic_minerals
  write (*, '(a, i4)') "  num aqueous complexes : ", sizes%num_aqueous_complexes
  write (*, '(a, i4)') "  num surface sites : ", sizes%num_surface_sites
  write (*, '(a, i4)') "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PrintSizes


! **************************************************************************** !
subroutine PrintState(reaction, state)

  use, intrinsic :: iso_c_binding

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (reaction_type), intent(in) :: reaction
  type (alquimia_state_f), intent(in) :: state

  ! local variables
  integer :: i
  real (c_double), pointer :: conc(:)

  write (*, '(a)') "state : "
  write (*, '(a, 1es13.6)') "  density water : ", state%density_water
  write (*, '(a, 1es13.6)') "  saturation : ", state%saturation
  write (*, '(a, 1es13.6)') "  porosity : ", state%porosity
  write (*, '(a, 1es13.6)') "  temperature : ", state%temperature
  write (*, '(a, 1es13.6)') "  aqueous pressure : ", state%aqueous_pressure
  write (*, '(a, i4, a)') "  total primary (", reaction%naqcomp, ") : "
  call c_f_pointer(state%total_primary, conc, (/reaction%naqcomp/))
  do i=1, reaction%naqcomp
     write (*, '(1es13.6)') conc(i)
  end do
end subroutine PrintState


! **************************************************************************** !
subroutine PrintMetaData(sizes, meta_data)

  use, intrinsic :: iso_c_binding

  use c_interface_module

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (alquimia_sizes_f), intent(in) :: sizes
  type (alquimia_meta_data_f), intent(in) :: meta_data

  ! local variables
  integer (c_int), pointer :: indices(:)
  type (c_ptr), pointer :: names(:)
  character(len=ALQUIMIA_MAX_STRING_LENGTH) :: name
  integer (c_int) :: i

  write (*, '(a)') "meta_data : "
  write (*, '(a, L1)') "  thread safe : ", meta_data%thread_safe
  write (*, '(a, L1)') "  temperature dependent : ", meta_data%temperature_dependent
  write (*, '(a, L1)') "  pressure dependent : ", meta_data%pressure_dependent
  write (*, '(a, L1)') "  porosity update : ", meta_data%porosity_update
  write (*, '(a, i4)') "  index base : ", meta_data%index_base
  write (*, '(a)') "  primary indices : "
  call c_f_pointer(meta_data%primary_indices, indices, (/sizes%num_primary/))
  do i=1, sizes%num_primary
     write (*, '(i4)') indices(i)
  end do
  call c_f_pointer(meta_data%primary_names, names, (/sizes%num_primary/))
  do i=1, sizes%num_primary
     call c_f_string_ptr(names(i), name)
     write (*, '(a)') trim(name)
  end do
end subroutine PrintMetaData


! **************************************************************************** !
subroutine PrintStatus(status)

  use, intrinsic :: iso_c_binding

  implicit none

#include "alquimia_containers.h90"

  ! function parameters
  type (alquimia_engine_status_f), intent(in) :: status

  write (*, '(a)') "status : "
  write (*, '(a)') "  num rhs evaluation  : ", status%num_rhs_evaluations
  write (*, '(a)') "  num jacobian evaluations : ", status%num_jacobian_evaluations
  write (*, '(a)') "  num newton iterations : ", status%num_newton_iterations
  write (*, '(a)') "  converged  : ", status%converged

end subroutine PrintStatus

end module PFloTranAlquimiaInterface_module

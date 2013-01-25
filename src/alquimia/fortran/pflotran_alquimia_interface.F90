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
       GetEngineMetaData

  private :: InitializePFloTranReactions, &
       ReadPFloTranConstraints, &
       ProcessPFloTranConstraint, &
       CopyAlquimiaToAuxVars, &
       CopyAuxVarsToAlquimia

  integer(kind=8), private, parameter :: integrity_check_value = 5784429996817932654_8
  !integer(kind=int64), private, parameter :: integrity_check_value = &
  !     b"0101000001000110011011000110111101010100011100100110000101101110"

  type, private :: PFloTranEngineState
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
     type (reactive_transport_auxvar_type), pointer :: rt_auxvar
     type (global_auxvar_type), pointer :: global_auxvar
     type(tran_constraint_list_type), pointer :: transport_constraints
     type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  end type PFloTranEngineState

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

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Input_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (c_ptr), intent(out) :: pft_engine_state
  type (AlquimiaSizes), intent(out) :: sizes
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFloTranEngineState), pointer :: engine_state
  PetscErrorCode :: ierr
  PetscBool :: option_found
  character(len=kAlquimiaMaxStringLength) :: string
  character(len=kAlquimiaMaxStringLength) :: filename_out
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(global_auxvar_type), pointer :: global_auxvar
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
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

  ! global_auxvar --> cell by cell temperature, pressure, saturation, density
  allocate(global_auxvar)
  call GlobalAuxVarInit(global_auxvar, option)

  ! rt_auxvar --> cell by cell chemistry data
  allocate(rt_auxvar)
  call RTAuxVarInit(rt_auxvar, reaction, option)

  ! assign default state values, not really needed?
  global_auxvar%pres = option%reference_pressure
  global_auxvar%temp = option%reference_temperature
  global_auxvar%den_kg = option%reference_water_density
  global_auxvar%sat = option%reference_saturation  

  ! NOTE(bja) : constraint coupler is only needed for printing the
  ! processed constraints to pflotran.out. But destroying the coupler
  ! destroys the auxvar as well, so we need to keep it around long
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
  engine_state%rt_auxvar => rt_auxvar
  engine_state%global_auxvar => global_auxvar
  engine_state%constraint_coupler => constraint_coupler
  engine_state%transport_constraints => transport_constraints

  pft_engine_state = c_loc(engine_state)

  !
  ! Grab sizes info that the driver needs to finish initializing
  !
  sizes%num_primary = reaction%ncomp
  sizes%num_sorbed = 0 ! FIXME(bja): need some sort of if (using_sorption)...
  sizes%num_kinetic_minerals = reaction%mineral%nkinmnrl
  sizes%num_aqueous_complexes = reaction%neqcplx
  sizes%num_surface_sites = reaction%surface_complexation%nsrfcplxrxn
  sizes%num_ion_exchange_sites = reaction%neqionxrxn

  !call PrintSizes(sizes)

  status%error = kAlquimiaNoError
  call f_c_string_ptr("Alquimia::PFloTran::Setup() : successful.", &
       status%message, kAlquimiaMaxStringLength)
  
end subroutine Setup


! **************************************************************************** !
subroutine Shutdown(pft_engine_state, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use AlquimiaContainers_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFloTranEngineState), pointer :: engine_state

  !write (*, '(a)') "PFloTranAlquimiaInterface::Shutdown() : "

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  call TranConstraintCouplerDestroy(engine_state%constraint_coupler)
  call TranConstraintDestroyList(engine_state%transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(engine_state%rt_auxvar)
  !call GlobalAuxVarDestroy(engine_state%global_auxvar)
  call ReactionDestroy(engine_state%reaction)
  call OptionDestroy(engine_state%option)

  deallocate(engine_state)
  nullify(engine_state)

  status%error = kAlquimiaNoError

end subroutine Shutdown


! **************************************************************************** !
subroutine ProcessCondition(pft_engine_state, condition, material_properties, &
     state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use AlquimiaContainers_module

  ! pflotran
  use String_module
  use Constraint_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaCondition), intent(in) :: condition
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent (inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type (AlquimiaConstraint), pointer :: local_constraints(:)
  type (AlquimiaConstraint), pointer :: constraint
  character (kAlquimiaMaxStringLength) :: name
  character (kAlquimiaMaxStringLength) :: constraint_type
  character (kAlquimiaMaxStringLength) :: associated_species
  real (c_double) :: constraint_value
  PetscReal :: porosity, volume
  integer :: i
  type(tran_constraint_type), pointer :: tran_constraint
  type(PFloTranEngineState), pointer :: engine_state

  !write (*, '(a)') "PFloTranAlquimiaInterface::ProcessCondition() : "

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  ! NOTE(bja): do NOT call CopyAlquimiaToAuxVars() here because most
  ! of that data is uninitialized in alquimia!
  engine_state%global_auxvar%den_kg = state%water_density
  engine_state%global_auxvar%sat = state%saturation
  engine_state%global_auxvar%temp = state%temperature
  engine_state%global_auxvar%pres = state%aqueous_pressure
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
          engine_state%global_auxvar, &
          engine_state%rt_auxvar, &
          tran_constraint, &
          engine_state%constraint_coupler)
     ! now repack the processed constraint data into the alquimia
     ! struct for the driver.
     call CopyAuxVarsToAlquimia( &
          engine_state%reaction, &
          engine_state%global_auxvar, &
          engine_state%rt_auxvar, &
          porosity, &
          state, aux_data)
  end if
  status%error = kAlquimiaNoError
end subroutine ProcessCondition


! **************************************************************************** !
subroutine ReactionStepOperatorSplit(pft_engine_state, &
     delta_t, material_properties, state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  real (c_double), intent(in) :: delta_t
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFloTranEngineState), pointer :: engine_state
  PetscReal :: porosity, volume, vol_frac_prim
  PetscReal :: tran_xx(state%total_primary%size)
  PetscInt :: i, phase_index

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !write (*, '(a)') "F_PFloTranAlquimiaInterface::ReactionStepOperatorSplit() :"

  !call PrintState(engine_state%reaction, state)

  call CopyAlquimiaToAuxVars(state, aux_data, material_properties, &
       engine_state%reaction, engine_state%global_auxvar, engine_state%rt_auxvar, &
       porosity, volume)

  ! copy total primaries into dummy transport variable
  phase_index = 1
  
  do i = 1, state%total_primary%size
     tran_xx(i) = engine_state%rt_auxvar%total(i, phase_index)
  enddo

  vol_frac_prim = 1.0

  engine_state%option%tran_dt = delta_t

  call RReact(engine_state%rt_auxvar, engine_state%global_auxvar, &
       tran_xx, volume, porosity, &
       status%num_newton_iterations, &
       engine_state%reaction, engine_state%option, vol_frac_prim)

  call RUpdateSolution(engine_state%rt_auxvar, engine_state%global_auxvar, &
       engine_state%reaction, engine_state%option)

  call CopyAuxVarsToAlquimia( &
       engine_state%reaction, &
       engine_state%global_auxvar, &
       engine_state%rt_auxvar, &
       porosity, &
       state, aux_data)

  status%error = kAlquimiaNoError

end subroutine ReactionStepOperatorSplit



! **************************************************************************** !
subroutine GetAuxiliaryOutput( &
     pft_engine_state, &
     material_properties, &
     state, &
     aux_data, &
     aux_output, &
     status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use AlquimiaContainers_module

  ! pflotran
  use Mineral_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaAuxiliaryOutputData), intent(inout) :: aux_output
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(PFloTranEngineState), pointer :: engine_state
  integer :: i, ph_index
  real (c_double), pointer :: local_array(:)
  PetscReal :: porosity, volume

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !write (*, '(a)') "PFloTranAlquimiaInterface::GetAuxiliaryOutput() :"

  ! NOTE(bja): right now, all info from the previous reaction step is
  ! still in the auxvars. We are assuming that the driver has called
  ! GetAuxiliaryOutput immediately after reaction step...!

  !call CopyAlquimiaToAuxVars(state, aux_data, material_properties, &
  !     engine_state%reaction, engine_state%global_auxvar, engine_state%rt_auxvar, &
  !     porosity, volume)

  ph_index = engine_state%reaction%species_idx%h_ion_id
  ! FIXME(bja): this violates the no geochemistry calculations in alquimia rule.
  aux_output%pH = -log10(engine_state%rt_auxvar%pri_act_coef(ph_index) * &
       engine_state%rt_auxvar%pri_molal(ph_index))

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

  status%error = kAlquimiaNoError
end subroutine GetAuxiliaryOutput


! **************************************************************************** !
subroutine GetEngineMetaData(pft_engine_state, meta_data, status)
  !  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use AlquimiaContainers_module

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  type (AlquimiaMetaData), intent(out) :: meta_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  integer (c_int), pointer :: local_indices(:)
  type (c_ptr), pointer :: name_list(:)
  character (len=kAlquimiaMaxWordLength), pointer :: pflotran_names(:)
  character (c_char), pointer :: name
  integer :: i, list_size
  type(PFloTranEngineState), pointer :: engine_state

  !write (*, '(a)') "PFloTran_Alquimia_GetEngineMetaData() :"

  call c_f_pointer(pft_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  ! TODO(bja) : can we extract this from pflotran without hardcoding?
  meta_data%thread_safe = .true.
  meta_data%temperature_dependent = .true.
  meta_data%pressure_dependent = .true.
  meta_data%porosity_update = engine_state%reaction%update_porosity
  meta_data%operator_splitting = .true.
  meta_data%global_implicit = .false.
  meta_data%index_base = 1

  !
  ! copy primary indices and names
  !

  if (meta_data%primary_names%size /= engine_state%reaction%ncomp) then
     write (*, '(a, i3, a, i3, a)') "meta_data%primary_names%size (", &
          meta_data%primary_names%size, ") != pflotran%reaction%ncomp(", &
          engine_state%reaction%ncomp, ")"
  end if
  list_size = meta_data%primary_names%size

  ! associate indices with the C memory
  call c_f_pointer(meta_data%primary_indices%data, local_indices, (/list_size/))
  ! NOTE(bja) : if the order in reaction%primary_species_names always
  ! is not correct we need to modifiy this...
  do i = 1, list_size
     local_indices(i) = i
  enddo

  pflotran_names => engine_state%reaction%primary_species_names

  call c_f_pointer(meta_data%primary_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name, kAlquimiaMaxStringLength)
     call f_c_string_chars(trim(pflotran_names(local_indices(i))), &
          name, kAlquimiaMaxStringLength)     
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

  ! associate indices with the C memory
  call c_f_pointer(meta_data%mineral_indices%data, local_indices, (/list_size/))
  ! NOTE(bja) : if the order in reaction%mineral_names always
  ! is not correct we need to modifiy this...
  do i = 1, list_size
     local_indices(i) = i
  enddo

  pflotran_names => engine_state%reaction%mineral%mineral_names

  call c_f_pointer(meta_data%mineral_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name, kAlquimiaMaxStringLength)
     call f_c_string_chars(trim(pflotran_names(local_indices(i))), &
          name, kAlquimiaMaxStringLength)     
  end do

  status%error = 0
end subroutine GetEngineMetaData


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
     global_auxvar, rt_auxvar, tran_constraint, constraint_coupler)
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
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvar
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvar
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
  call ReactionEquilibrateConstraint(rt_auxvar, global_auxvar, reaction, &
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
  constraint_coupler%global_auxvar => global_auxvar
  constraint_coupler%rt_auxvar => rt_auxvar     
  call ReactionPrintConstraint(constraint_coupler, reaction, option)


end subroutine ProcessPFloTranConstraint


! **************************************************************************** !
subroutine CopyAlquimiaToAuxVars(state, aux_data, material_prop, &
  reaction, global_auxvar, rt_auxvar, porosity, volume)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none

  ! function parameters
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaMaterialProperties), intent(in) :: material_prop
  type(reaction_type), pointer, intent(inout) :: reaction
  type(global_auxvar_type), pointer, intent(inout) :: global_auxvar
  type(reactive_transport_auxvar_type), pointer, intent(inout) :: rt_auxvar
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
  global_auxvar%den_kg(1) = state%water_density
  global_auxvar%sat(1) = state%saturation
  global_auxvar%temp(1) = state%temperature
  global_auxvar%pres(1) = state%aqueous_pressure

  porosity = state%porosity
  volume = material_prop%volume

  !
  ! primary aqueous
  !
  call c_f_pointer(state%total_primary%data, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvar%total(i, phase_index) = local_array(i)
  end do

  call c_f_pointer(state%free_ion%data, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvar%pri_molal(i) = local_array(i)
  end do

  call c_f_pointer(aux_data%primary_activity_coeff%data, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     rt_auxvar%pri_act_coef(i) = local_array(i)
  end do

  !
  ! aqueous complexes
  !
  call c_f_pointer(aux_data%secondary_activity_coeff%data, local_array, &
       (/reaction%neqcplx/))
  do i = 1, reaction%neqcplx
     rt_auxvar%sec_act_coef(i) = local_array(i)
  end do

  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction%data, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     rt_auxvar%mnrl_volfrac(i) = local_array(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area%data, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     rt_auxvar%mnrl_area(i) = local_array(i)
  end do

  !
  ! ion exchange, CEC only present in reaction, not aux_vars?
  !
  call c_f_pointer(state%cation_exchange_capacity%data, local_array, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     reaction%eqionx_rxn_CEC(i) = local_array(i)
  end do

  call c_f_pointer(aux_data%ion_exchange_ref_cation_conc%data, local_array, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     rt_auxvar%eqionx_ref_cation_sorbed_conc(i) = local_array(i)
  end do


  !
  ! equilibrium surface complexation
  !
  call c_f_pointer(state%surface_site_density%data, local_array, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     reaction%surface_complexation%srfcplxrxn_site_density(i) = local_array(i)
  end do

  call c_f_pointer(aux_data%surface_complex_free_site_conc%data, local_array, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     rt_auxvar%srfcplxrxn_free_site_conc(i) = local_array(i)
  end do

  !
  ! isotherms, TODO(bja): copy out of material_props into reaction (not auxvar)
  !


end subroutine CopyAlquimiaToAuxVars

! **************************************************************************** !
subroutine CopyAuxVarsToAlquimia(reaction, global_auxvar, rt_auxvar, &
     porosity, state, aux_data)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module

  ! pflotran
  use Reaction_aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none

  ! function parameters
  type(reaction_type), pointer, intent(in) :: reaction
  type(global_auxvar_type), pointer, intent(in) :: global_auxvar
  type(reactive_transport_auxvar_type), pointer, intent(in) :: rt_auxvar
  PetscReal, intent(in) :: porosity
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  real (c_double), pointer :: local_array(:)
  integer :: i, phase_index

  phase_index = 1 ! TODO(bja): grab from pflotran?

  !write (*, '(a)') "PFloTran_Alquimia_CopyAuxVarsToAlquimia() :"

  !
  ! state
  !
  state%water_density = global_auxvar%den_kg(1)
  state%saturation = global_auxvar%sat(1)
  state%temperature = global_auxvar%temp(1)
  state%aqueous_pressure = global_auxvar%pres(1)

  state%porosity = porosity

  !
  ! primary aqueous species
  !
  call c_f_pointer(state%total_primary%data, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     local_array(i) = rt_auxvar%total(i, phase_index)
  end do

  call c_f_pointer(state%free_ion%data, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     local_array(i) = rt_auxvar%pri_molal(i)
  end do

  call c_f_pointer(aux_data%primary_activity_coeff%data, local_array, (/reaction%naqcomp/))
  do i = 1, reaction%naqcomp
     local_array(i) = rt_auxvar%pri_act_coef(i)
  end do

  !
  ! secondary aqueous complexes
  !
  call c_f_pointer(aux_data%secondary_activity_coeff%data, local_array, &
       (/reaction%neqcplx/))
  do i = 1, reaction%neqcplx
     local_array(i) = rt_auxvar%sec_act_coef(i)
  end do

  !
  ! sorbed
  !
  call c_f_pointer(state%total_sorbed%data, local_array, (/reaction%neqsorb/))
  do i = 1, reaction%neqsorb
     local_array(i) = rt_auxvar%total_sorb_eq(i)
  end do

  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction%data, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     local_array(i) = rt_auxvar%mnrl_volfrac(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area%data, local_array, &
       (/reaction%mineral%nkinmnrl/))
  do i = 1, reaction%mineral%nkinmnrl
     local_array(i) = rt_auxvar%mnrl_area(i)
  end do

  !
  ! ion exchange, CEC only present in reaction, not aux_vars?
  !
  call c_f_pointer(state%cation_exchange_capacity%data, local_array, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     local_array(i) = reaction%eqionx_rxn_CEC(i)
  end do

  call c_f_pointer(aux_data%ion_exchange_ref_cation_conc%data, local_array, &
       (/reaction%neqionxrxn/))
  do i = 1, reaction%neqionxrxn
     local_array(i) = rt_auxvar%eqionx_ref_cation_sorbed_conc(i)
  end do


  !
  ! equilibrium surface complexation, site density in reaction, not aux_vars
  !
  call c_f_pointer(state%surface_site_density%data, local_array, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     local_array(i) = reaction%surface_complexation%srfcplxrxn_site_density(i)
  end do

  call c_f_pointer(aux_data%surface_complex_free_site_conc%data, local_array, &
       (/reaction%surface_complexation%nsrfcplxrxn/))
  do i = 1, reaction%surface_complexation%nsrfcplxrxn
     local_array(i) = rt_auxvar%srfcplxrxn_free_site_conc(i)
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

  use AlquimiaContainers_module

  implicit none

  ! function parameters
  type (AlquimiaSizes), intent(in) :: sizes

  write (*, '(a)') "size : "
  write (*, '(a, i4)') "  num primary : ", sizes%num_primary
  write (*, '(a, i4)') "  num kinetics minerals : ", sizes%num_kinetic_minerals
  write (*, '(a, i4)') "  num aqueous complexes : ", sizes%num_aqueous_complexes
  write (*, '(a, i4)') "  num surface sites : ", sizes%num_surface_sites
  write (*, '(a, i4)') "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PrintSizes


! **************************************************************************** !
subroutine PrintState(state)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module

  implicit none

  ! function parameters
  type (AlquimiaState), intent(in) :: state

  ! local variables
  integer :: i
  real (c_double), pointer :: conc(:)

  write (*, '(a)') "state : "
  write (*, '(a, 1es13.6)') "  density water : ", state%water_density
  write (*, '(a, 1es13.6)') "  saturation : ", state%saturation
  write (*, '(a, 1es13.6)') "  porosity : ", state%porosity
  write (*, '(a, 1es13.6)') "  temperature : ", state%temperature
  write (*, '(a, 1es13.6)') "  aqueous pressure : ", state%aqueous_pressure
  write (*, '(a, i4, a)') "  total primary (", state%total_primary%size, ") : "
  call c_f_pointer(state%total_primary%data, conc, (/state%total_primary%size/))
  do i=1, state%total_primary%size
     write (*, '(1es13.6)') conc(i)
  end do
end subroutine PrintState


! **************************************************************************** !
subroutine PrintMetaData(meta_data)

  use, intrinsic :: iso_c_binding

  use c_interface_module

  use AlquimiaContainers_module

  implicit none

  ! function parameters
  type (AlquimiaMetaData), intent(in) :: meta_data

  ! local variables
  integer (c_int), pointer :: indices(:)
  type (c_ptr), pointer :: names(:)
  character(len=kAlquimiaMaxStringLength) :: name
  integer (c_int) :: i

  write (*, '(a)') "meta_data : "
  write (*, '(a, L1)') "  thread safe : ", meta_data%thread_safe
  write (*, '(a, L1)') "  temperature dependent : ", meta_data%temperature_dependent
  write (*, '(a, L1)') "  pressure dependent : ", meta_data%pressure_dependent
  write (*, '(a, L1)') "  porosity update : ", meta_data%porosity_update
  write (*, '(a, i4)') "  index base : ", meta_data%index_base
  write (*, '(a)') "  primary indices : "
  call c_f_pointer(meta_data%primary_indices%data, indices, (/meta_data%primary_indices%size/))
  do i=1, meta_data%primary_indices%size
     write (*, '(i4)') indices(i)
  end do
  call c_f_pointer(meta_data%primary_names%data, names, (/meta_data%primary_names%size/))
  do i=1, meta_data%primary_names%size
     call c_f_string_ptr(names(i), name)
     write (*, '(a)') trim(name)
  end do
end subroutine PrintMetaData


! **************************************************************************** !
subroutine PrintStatus(status)

  use, intrinsic :: iso_c_binding

  use AlquimiaContainers_module

  implicit none

  ! function parameters
  type (AlquimiaEngineStatus), intent(in) :: status

  write (*, '(a)') "status : "
  write (*, '(a)') "  num rhs evaluation  : ", status%num_rhs_evaluations
  write (*, '(a)') "  num jacobian evaluations : ", status%num_jacobian_evaluations
  write (*, '(a)') "  num newton iterations : ", status%num_newton_iterations
  write (*, '(a)') "  converged  : ", status%converged

end subroutine PrintStatus

end module PFloTranAlquimiaInterface_module

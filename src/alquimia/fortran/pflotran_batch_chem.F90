module BatchChem

  implicit none

  private

#include "definitions.h"

  public :: BatchChemInitializeReactions, &
            BatchChemReadConstraints, &
            BatchChemProcessConstraint, &
            BatchChemCopyAuxVarsToState

contains

subroutine BatchChemInitializeReactions(option, input, reaction)

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

end subroutine BatchChemInitializeReactions

! **************************************************************************** !
!
! BatchChemReadConstraints
!
! We are just reading the data from the input file. No processing is
! done here because we don't know yet if we are using these.
!
! **************************************************************************** !

subroutine BatchChemReadConstraints(option, input, reaction, &
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

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_list_type), pointer :: transport_constraints

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

end subroutine BatchChemReadConstraints


subroutine BatchChemProcessConstraint(option, reaction, &
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

  type(option_type), pointer :: option
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

end subroutine BatchChemProcessConstraint

! **************************************************************************** !
subroutine BatchChemCopyAuxVarsToState(reaction, &
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

end subroutine BatchChemCopyAuxVarsToState

end module BatchChem


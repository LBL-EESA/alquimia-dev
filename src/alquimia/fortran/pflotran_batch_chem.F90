module BatchChem

  implicit none

  private

#include "definitions.h"

  public :: BatchChemInitializeReactions, &
            BatchChemProcessConstraints

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


subroutine BatchChemProcessConstraints(option, input, reaction, &
     global_auxvars, rt_auxvars, transport_constraints, constraint_coupler)

  use Reaction_module
  use Reaction_Aux_module
  use Database_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Constraint_module
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
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  PetscBool :: use_prev_soln_as_guess
  PetscInt :: num_iterations

  ! initialize memory
  allocate(transport_constraints)
  call TranConstraintInitList(transport_constraints)
  allocate(constraint_coupler)
  constraint_coupler => TranConstraintCouplerCreate(option)


  !
  ! read the constraints...
  !

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

  !
  ! process constraints
  !
  num_iterations = 0
  use_prev_soln_as_guess = PETSC_FALSE
  tran_constraint => transport_constraints%first
  ! NOTE(bja): we only created one set of global and rt auxvars, so if
  ! there is more than one constratint in the input file, they will be
  ! over written.
  do 
     if (.not. associated(tran_constraint)) exit
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

     ! link the constraint to the constraint coupler
     constraint_coupler%constraint_name = tran_constraint%name
     constraint_coupler%aqueous_species => tran_constraint%aqueous_species
     constraint_coupler%minerals => tran_constraint%minerals
     constraint_coupler%surface_complexes => tran_constraint%surface_complexes
     constraint_coupler%colloids => tran_constraint%colloids
     constraint_coupler%global_auxvar => global_auxvars
     constraint_coupler%rt_auxvar => rt_auxvars
     
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
     call ReactionPrintConstraint(constraint_coupler, reaction, option)
     tran_constraint => tran_constraint%next
  enddo

end subroutine BatchChemProcessConstraints


end module BatchChem


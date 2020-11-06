
!
! Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
! through Lawrence Berkeley National Laboratory (subject to receipt of any 
! required approvals from the U.S. Dept. of Energy).  All rights reserved.
! 
! Alquimia is available under a BSD license. See LICENSE.txt for more
! information.
!
! If you have questions about your rights to use or distribute this software, 
! please contact Berkeley Lab's Technology Transfer and Intellectual Property 
! Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
! 
! NOTICE.  This software was developed under funding from the U.S. Department 
! of Energy.  As such, the U.S. Government has been granted for itself and 
! others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
! license in the Software to reproduce, prepare derivative works, and perform 
! publicly and display publicly.  Beginning five (5) years after the date 
! permission to assert copyright is obtained from the U.S. Department of Energy, 
! and subject to any subsequent five (5) year renewals, the U.S. Government is 
! granted for itself and others acting on its behalf a paid-up, nonexclusive, 
! irrevocable, worldwide license in the Software to reproduce, prepare derivative
! works, distribute copies to the public, perform publicly and display publicly, 
! and to permit others to do so.
! 
! Authors: Benjamin Andre <bandre@lbl.gov>
!

! **************************************************************************** !
!
! CruncFlow Alquimia Inteface module
!
! Authors: 
!        interface: Sergi Molins <smolins@lbl.gov>
!        CrunchFlow: Carl Steefel
!        template (plfotran_alquimia_interface): Benjamin Andre
!
! Notes:
!
!  * Public function call signatures, including intent, are dictated
!    by the alquimia API.
!
!  * alquimia data structures defined in the AlquimiaContainers_module
!    (alquimia_containers.F90) are dictated by the alquimia API.
!
!  * All function calls involving crunchflow native data structures must
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
!    crunch_alquimia_wrappers.F90!
!
!  * (smr) 2014-08 - crunchflow variables are global and accessed 
!    through modules - not thread safe!
!
! **************************************************************************** !

module CrunchAlquimiaInterface_module

  ! crunchflow modules if needed
  use crunchtype

  implicit none

  public :: &
       Setup, &
       Shutdown, &
       ProcessCondition, &
       ReactionStepOperatorSplit, &
       GetAuxiliaryOutput, &
       GetProblemMetaData

  private :: &
       SetupCrunchOptions, &
       SetEngineFunctionality, &
       SetAlquimiaSizes, &
       ProcessCrunchConstraint, &
       ConvertAlquimiaConditionToCrunch, &
       CopyAlquimiaToAuxVars, &
       CopyAuxVarsToAlquimia, &
       GetAuxiliaryDataSizes, &
       PackAlquimiaAuxiliaryData, &
       UnpackAlquimiaAuxiliaryData, &
!       PrintTranConstraint, &
!       PrintAqueousSpeciesConstraint, &
       PrintSizes, &
       PrintState, &
       PrintEngineFunctionality, &
       PrintProblemMetaData, &
       PrintStatus, &
       PrintAuxiliaryData

  integer(kind=8), private, parameter :: integrity_check_value = 1793265457844299968_8 ! _8 is an 8-byte integer constant, make different from other engines

  ! TO DO
  ! add this so that CF can keep track of MPI parameters and print only if rank == 0
  !type, private :: option_type
  !
  !  PetscMPIInt :: global_comm             ! MPI_COMM_WORLD
  !  PetscMPIInt :: global_rank             ! rank in MPI_COMM_WORLD
  !  PetscMPIInt :: global_commsize         ! size of MPI_COMM_WORLD
  !  PetscMPIInt :: global_group            ! id of group for MPI_COMM_WORLD

  !  PetscMPIInt :: mycomm                  ! PETSC_COMM_WORLD
  !  PetscMPIInt :: myrank                  ! rank in PETSC_COMM_WORLD
  !  PetscMPIInt :: mycommsize              ! size of PETSC_COMM_WORLD
  !  PetscMPIInt :: mygroup                 ! id of group for PETSC_COMM_WORLD
  !  PetscMPIInt :: io_rank                 ! I/O rank (=0)

  !end type option_type

  type, private :: CrunchEngineState
     ! This is the data structure that stores the persistent data for
     ! crunchflow, (e.g. reaction network).
     !
     ! We pass it back and forth with c as a void pointer so we don't
     ! have to a global variable.
     !
     ! It is NOT part of the alquimia interface, and the driver code
     ! should not use or depend on it!
     !
     ! NOTE(bja): these are fortran pointers, so this struct can not
     ! be unpacked on the c side!
     !
     ! NOTE(smr): crunchflow uses global variable defined in modules
     ! most variables are accessed that way 
     ! mostly only sizes, flags and pointer indexes are part of CF's engine_state
     !
     integer(kind=8)          :: integrity_check
     integer(i4b)             :: ncomp
     integer(i4b)             :: nspec
     integer(i4b)             :: nkin
     integer(i4b)             :: nrct
     integer(i4b)             :: ngas
     integer(i4b)             :: npot
     integer(i4b)             :: nexchange
     integer(i4b)             :: nexch_sec
     integer(i4b)             :: nsurf
     integer(i4b)             :: nsurf_sec
     integer(i4b)             :: ndecay
     integer(i4b)             :: ikin
     integer(i4b)             :: neqn
     integer(i4b)             :: nretard
     integer(i4b)             :: nx
     integer(i4b)             :: ny
     integer(i4b)             :: nz
     integer(i4b)             :: igamma
     integer(i4b)             :: ikph
     integer(i4b)             :: jpor
     real(dp)                 :: corrmax
     real(dp)                 :: deltmin
     real(dp)                 :: time
     logical                   :: hands_off
  end type CrunchEngineState
  
contains


! **************************************************************************** !
!
! Public interface routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine Setup(input_filename, hands_off, cf_engine_state, sizes, functionality, status)
!  NOTE: Function signature is dictated by the alquimia API.
!
!  NOTE: Assumes that MPI_Init() and / or PetscInitialize() have already
!    been called by the driver

  use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_bool

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow global variables
  use crunchtype
  use params, only: mls
  use medium, only: isaturate

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  logical (c_bool), value, intent(in) :: hands_off
  type (c_ptr), intent(out) :: cf_engine_state
  type (AlquimiaSizes), intent(out) :: sizes
  type (AlquimiaEngineFunctionality), intent(out) :: functionality
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(CrunchEngineState), pointer :: engine_state  
  character(len=kAlquimiaMaxStringLength) :: string
  !
  ! these are pointers to crunchflow's variables that go to CF's state
  !
  ! geochemical sytem sizes, domain sizes
  !
  integer(i4b)          :: ncomp
  integer(i4b)          :: nspec
  integer(i4b)          :: nkin
  integer(i4b)          :: nrct
  integer(i4b)          :: ngas
  integer(i4b)          :: npot
  integer(i4b)          :: nexchange
  integer(i4b)          :: nexch_sec
  integer(i4b)          :: nsurf
  integer(i4b)          :: nsurf_sec
  integer(i4b)          :: ndecay
  integer(i4b)          :: ikin
  integer(i4b)          :: neqn
  integer(i4b)          :: nretard
  integer(i4b)          :: nx
  integer(i4b)          :: ny
  integer(i4b)          :: nz
  !
  ! flags, pointer indexes, parameters
  !
  integer(i4b)          :: igamma ! gamma coeff update type
  integer(i4b)          :: ikph   ! points to pH
  integer(i4b)          :: jpor   ! porosity update type (0 is constant)
  real(dp)              :: corrmax ! max change in unknown in log units
  real(dp)              :: deltmin ! minimum delta time
  !
  ! this are dummy arguments required by start98's signature
  ! but not needed for now as part of cf_engine_state
  !
  integer(i4b) :: ipath  ! reaction_path calculations
  integer(i4b) :: ikmast ! points to master species
  integer(i4b) :: ikO2   ! points to O2
  integer(i4b) :: nstop  ! # output times for CF
  integer(i4b) :: nseries ! time series print for CF   
  integer(i4b) :: str_mon ! time keeping for CF
  integer(i4b) :: str_day
  integer(i4b) :: str_hr
  integer(i4b) :: str_min
  integer(i4b) :: str_sec
  integer(i4b) :: NumInputFiles = 1    ! CF restart control param
  integer(i4b) :: InputFileCounter = 1 ! CF restart control param
  character(len=mls) :: data1 ! filename database file  
  character(len=2*mls) :: ltitle ! in out for start98
  real(dp) :: tstep    ! time step control
  real(dp) :: delt     ! delta time
  real(dp) :: ttol     ! time tolerance
  logical, parameter :: alquimia=.true. ! alquimia flag
! end of required variables for Crunchflow's start98 signature

! not crunch-related  
  integer :: len
  integer :: i
 
  integer(i4b) :: jx
  integer(i4b) :: jy
  integer(i4b) :: jz
  
  ! write (*, '(a)') "CrunchFlow_Alquimia_Setup() : "

  ! setup CrunchFlow options
  call SetupCrunchOptions(input_filename) !, inputfilename) !, option)

  ! initialize parameters for CrunchFlow (comment out some things for now, likely not needed)
  !!FirstCall = .TRUE.

  !!trans = 'N'
  !!solve_flag = .false.
  !!steady = .false.
  !!os3dpetsc = .FALSE.
  !!pi = DACOS(-1.0D0)
  ncomp = 0
  ngas = 0
  nspec = 0
  ndecay = 0
  nsurf = 0
  nexchange = 0
  nexch_sec = 0
  nsurf_sec = 0
  nrct = 0
  nkin = 0
  ikin = 0
  !!NumSourceTerms = 0
  isaturate = 0  
  !!dtmaxcour = 0.0   
  !!iprnt = 0     

  ! initialize CrunchFlow
  CALL StartTope(ncomp,nspec,nkin,nrct,ngas,npot,nx,ny,nz,data1,ipath,igamma,  &
               ikmast,ikph,iko2,ltitle,tstep,delt,deltmin,ttol,jpor,ikin,nstop,       &
               corrmax,nseries,nexchange,nexch_sec,nsurf,nsurf_sec,ndecay,str_mon,    &
               str_day,str_hr,str_min,str_sec,NumInputFiles,InputFileCounter)

  ! number of unknowns
  neqn = ncomp + nexchange + nsurf + npot
    
  ! allocate vars for OS3D, gases too if needed
  call AllocateOS3D(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
  if (isaturate == 1) call AllocateGasesOS3D(nx,ny,nz,ncomp)

  ! allocate all other needed variables
  call AllocateALL(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)

  ! Alquimia sizes
  call SetAlquimiaSizes(ncomp, nspec, nkin, nrct, ngas, &
                        nexchange, nsurf, ndecay, npot, &
                        nretard, ikin, sizes)

  ! Engine functionality
  ! note: no need to pass anything for now
  call SetEngineFunctionality(jpor,functionality)
  !
  ! save crunch's persistent data to a struct so the driver can store it for us
  ! crunch uses global variables so there is no need for anything to be stored for now
  !
  allocate(engine_state)
  ! integrity check value
  engine_state%integrity_check = integrity_check_value
  ! geochemical system sizes
  engine_state%ncomp = ncomp
  engine_state%nspec = nspec
  engine_state%nkin = nkin
  engine_state%nrct = nrct
  engine_state%ngas = ngas
  engine_state%npot = npot
  engine_state%nexchange = nexchange
  engine_state%nexch_sec = nexch_sec
  engine_state%nsurf = nsurf
  engine_state%nsurf_sec = nsurf_sec
  engine_state%ndecay = ndecay
  engine_state%ikin = ikin
  engine_state%neqn = neqn
  engine_state%nretard = nretard
  ! domain size (start98-alquimia sets it to nx=ny=nz=1)
  engine_state%nx = nx
  engine_state%ny = ny
  engine_state%nz = nz
  ! flags, pointer indexes, parameters
  engine_state%igamma = igamma
  engine_state%ikph = ikph
  engine_state%jpor = jpor
  engine_state%corrmax = corrmax
  engine_state%deltmin = deltmin
  engine_state%time = 0.0d0 ! hardwired -- this disables lag time for microbial growth

  engine_state%hands_off = hands_off
  
  cf_engine_state = c_loc(engine_state)

  status%error = kAlquimiaNoError
  call f_c_string_ptr("Alquimia::Crunch::Setup() : successful.", &
       status%message, kAlquimiaMaxStringLength)
  
end subroutine Setup


! **************************************************************************** !
subroutine Shutdown(cf_engine_state, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow
  use crunchtype

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(CrunchEngineState), pointer :: engine_state

  !write (*, '(a)') "CrunchAlquimiaInterface::Shutdown() : "

  call c_f_pointer(cf_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

! nothing to destroy for crunchflow for now

  deallocate(engine_state)
  nullify(engine_state)

  status%error = kAlquimiaNoError

end subroutine Shutdown


! **************************************************************************** !
subroutine ProcessCondition(cf_engine_state, condition, properties, &
     state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  use c_f_interface_module, only : c_f_string_ptr, f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow
  use crunchtype
  
  use runtime, only: nchem
  use concentration, only: condlabel, stmp, jinit

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaGeochemicalCondition), intent(in) :: condition
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent (inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  real (c_double), pointer :: data(:)
  character (kAlquimiaMaxStringLength) :: name
  real (c_double) :: constraint_value
  !!!real(dp) :: porosity, volume
  integer :: i
  type(CrunchEngineState), pointer :: engine_state
  logical, parameter :: copy_auxdata = .false.

  integer(i4b) :: found
  integer(i4b) :: nco
  integer(i4b) :: nlen
  !
  ! geochemical sytem sizes, domain sizes
  !
  integer(i4b)          :: ncomp
  integer(i4b)          :: nspec
  integer(i4b)          :: nkin
  integer(i4b)          :: nrct
  integer(i4b)          :: ngas
  integer(i4b)          :: npot
  integer(i4b)          :: nexchange
  integer(i4b)          :: nexch_sec
  integer(i4b)          :: nsurf
  integer(i4b)          :: nsurf_sec
  integer(i4b)          :: ndecay
  integer(i4b)          :: ikin
  integer(i4b)          :: neqn
  integer(i4b)          :: nretard

  !write (*, '(a)') "CrunchAlquimiaInterface::ProcessCondition() : "

  call c_f_pointer(cf_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  ! geochemical system sizes
  ncomp = engine_state%ncomp
  nspec = engine_state%nspec  
  nkin = engine_state%nkin 
  nrct = engine_state%nrct 
  ngas = engine_state%ngas 
  npot = engine_state%npot 
  nexchange = engine_state%nexchange
  nexch_sec = engine_state%nexch_sec
  nsurf = engine_state%nsurf
  nsurf_sec = engine_state%nsurf_sec
  ndecay = engine_state%ndecay
  ikin = engine_state%ikin
  neqn = engine_state%neqn
  nretard = engine_state%nretard

  ! NOTE(bja): the data stored in alquimia's aux_data is uninitialized
  ! at this point, so don't want to copy it! (copy_auxdata = false)
  !!call CopyAlquimiaToAuxVars(copy_auxdata, state, aux_data, material_properties, &
  !!                           ncomp, nspec, nkin, nrct, ngas, nexchange, nsurf, ndecay, npot, nretard)
  !
  ! process the condition
  !
  call  c_f_string_ptr(condition%name, name)
  !! note: need to replace with crunchflow print message call for rank 0!!
  !!engine_state%option%io_buffer = "processing : " // trim(name)
  !!call printMsg(engine_state%option)


  nlen = LEN(name)
  call majuscules(name,nlen)

  if (condition%aqueous_constraints%size > 0) then
     ! the driver is supplying the constraint data, so we need to
     ! construct a geochemical condition in crunchflow's internal format.
 
     ! check whether this condition has already been added to list
     ! prepend 'alq_' to name coming from driver (this is how they are
     ! stored to avoid duplication of names -- easier than checking)
     found = 0
     do_crunchConditions: do nco = 1,nchem

       if (trim(condlabel(nco)) == 'alq_'//trim(adjustl(name))) then         
         
         found = nco
         exit do_crunchConditions
       end if

     end do do_crunchConditions

!    the condition has not be added yet, go ahead and add it now
     if (found == 0) then
       call ConvertAlquimiaConditionToCrunch(state,properties, &
                                             ncomp,nspec,nrct,nkin, &
                                             ngas,nsurf,nexchange, &
                                             ikin,nexch_sec,nsurf_sec, &
                                             condition)

        nco = nchem
        found = nco

     end if

     !! (smr) fixme: not supported yet
!!     status%error = kAlquimiaErrorUnsupportedFunctionality
!!     call f_c_string_ptr("ERROR: crunchflow interface does not support externally supplied geochemical conditions!", &
!!          status%message, kAlquimiaMaxStringLength)
!!     return
  
  else
     ! the driver just supplied a name, so we check for a constraint
     ! with that name in the crunchflow input file and use that.
     
     !! note: need to replace with crunchflow print message call for rank 0!!
     !!engine_state%option%io_buffer = "Looking for crunchflow constraint : " // trim(name)
     !!call printMsg(engine_state%option)
     
     found = 0
     do_crunchCond: do nco = 1,nchem     

       if (trim(condlabel(nco)) == trim(name)) then         
         found = nco
         exit do_crunchCond
       end if

     end do do_crunchCond

     if (found == 0) then
        
        ! end of the list without out finding a match.
         call f_c_string_ptr("INPUT ERROR: Could not find crunchflow constraint : "//trim(name), &
              status%message, kAlquimiaMaxStringLength)
              status%error = kAlquimiaErrorUnknownConstraintName
         return

     end if 

  end if

  ! move here for now, we need to know now what condition we are dealing with
  jinit(1,1,1) = nco
  call CopyAlquimiaToAuxVars(copy_auxdata, engine_state%hands_off, &
                             state, aux_data, properties, &
                             ncomp, nspec, nkin, nrct, ngas, nexchange, nsurf, ndecay, npot, nretard)

  if (allocated(stmp)) DEALLOCATE(stmp) !! to revisit ?
  allocate(stmp(ncomp))
  stmp = 0.0d0

! call to crunchflow subroutines are packaged in ProcessCrunchConstraint for clarity
!
  call ProcessCrunchConstraint(engine_state, nco)

  if (allocated(stmp)) DEALLOCATE(stmp) !! to revisit ?
!
! repack the processed constraint data into the alquimia struct for the driver
!
  call CopyAuxVarsToAlquimia(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             nretard, state, aux_data )
  status%error = kAlquimiaNoError

end subroutine ProcessCondition


! **************************************************************************** !
subroutine ReactionStepOperatorSplit(cf_engine_state, &
     delta_t, properties, state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_double, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow
  use crunchtype
  use params, only: secyr
  use runtime, only: iterat, Duan
  use medium, only: isaturate
  use concentration, only: ulab, &
                           sp10, sp, spold, &
                           s, sn, &
                           spex, spexold, &
                           spsurf10, spsurf, spsurfold, &
                           xgram, xgramOld, stmp
  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  real (c_double), value, intent(in) :: delta_t
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(CrunchEngineState), pointer :: engine_state
  
  ! local pointer to crunchflow's state
  !
  ! geochemical sytem sizes, domain sizes
  !
  integer(i4b)          :: ncomp
  integer(i4b)          :: nspec
  integer(i4b)          :: nkin
  integer(i4b)          :: nrct
  integer(i4b)          :: ngas
  integer(i4b)          :: npot
  integer(i4b)          :: nexchange
  integer(i4b)          :: nexch_sec
  integer(i4b)          :: nsurf
  integer(i4b)          :: nsurf_sec
  integer(i4b)          :: ndecay
  integer(i4b)          :: ikin
  integer(i4b)          :: neqn
  integer(i4b)          :: nretard
  integer(i4b)          :: nx
  integer(i4b)          :: ny
  integer(i4b)          :: nz
  !
  ! flags, pointer indexes, parameters
  !
  integer(i4b)          :: igamma ! gamma coeff update type
  integer(i4b)          :: ikph   ! points to pH
  integer(i4b)          :: jpor   ! porosity update type (0 is constant)
  real(dp)              :: corrmax ! max unknown update log scale
  real(dp)              :: deltmin ! minimum delta time
  real(dp)              :: time    ! time
  
  ! other local
  real(dp) :: porosity, volume, vol_frac_prim    
  integer(i4b) :: jx, jy, jz
  integer(i4b) :: i, is, ix, ik
  logical, parameter :: copy_auxdata = .true.

  ! crunchflow local
  integer(i4b)                :: newtmax
  integer(i4b)                :: icvg
  integer(i4b)                :: ineg
  real(dp)                    :: AqueousToBulk
  real(dp)                    :: delt
  real(dp)                    :: dtnewest
  
  call c_f_pointer(cf_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !write (*, '(a)') "F_CrunchAlquimiaInterface::ReactionStepOperatorSplit() :"
  !write(*,*)'------- GOING IN -----------------------------------------'
  !call PrintState(state)

  ! 1st: GET ENGINE STATE FROM ALQUIMIA AND UPDATE CRUNCHFLOW VARIABLES WITH ALQUIMIA STATE
  
  ! geochemical system sizes
  ncomp = engine_state%ncomp
  nspec = engine_state%nspec  
  nkin = engine_state%nkin 
  nrct = engine_state%nrct 
  ngas = engine_state%ngas 
  npot = engine_state%npot 
  nexchange = engine_state%nexchange
  nexch_sec = engine_state%nexch_sec
  nsurf = engine_state%nsurf
  nsurf_sec = engine_state%nsurf_sec
  ndecay = engine_state%ndecay
  ikin = engine_state%ikin
  neqn = engine_state%neqn
  nretard = engine_state%nretard
  ! domain size (start98-alquimia sets it to nx=ny=nz=1)
  nx = engine_state%nx
  ny = engine_state%ny
  nz = engine_state%nz
  ! flags, pointer indexes, parameters
  igamma = engine_state%igamma
  jpor = engine_state%jpor
  corrmax = engine_state%corrmax
  deltmin = engine_state%deltmin
  time = engine_state%time

  !call PrintAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
  !                        nexchange, nsurf, ndecay, npot, &
  !                        aux_data)


  call CopyAlquimiaToAuxVars(copy_auxdata,   engine_state%hands_off, &
                             state, aux_data, properties, &
                             ncomp, nspec, nkin, nrct, ngas, nexchange, nsurf, ndecay, npot, nretard)

  IF (nexchange > 0) THEN
    CALL UpdateExchanger(nx,ny,nz,nexchange)
  END IF

! set delt variable in CrunchFlow units
! IN= seconds, CF=years
  delt = delta_t / secyr

! this is allocated inside the calculations, not as part of setup
  if (allocated(stmp)) DEALLOCATE(stmp)

! 2ND: CALL OS3D_NEWTON

! straight from CrunchFlow's CrunchFlow.f90 subroutine, except where noted ------------------------
    newtmax = 0
    loop4001: DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx          
          CALL keqcalc2(ncomp,nrct,nspec,ngas,nsurf_sec,jx,jy,jz)
          IF (igamma == 3) THEN

            IF (Duan) THEN
              CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
            ELSE
              CALL gamma(ncomp,nspec,jx,jy,jz)
            END IF

          END IF
          
          CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)

          IF (ncomp == 1 .AND. ulab(1) == 'Tracer') THEN
             sp(1,jx,jy,jz) = DLOG(sn(1,jx,jy,jz))
             sp10(1,jx,jy,jz) = sn(1,jx,jy,jz)
          END IF

          CALL os3d_newton(ncomp,nspec,nkin,nrct,ngas,ikin,           &
            nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,neqn,igamma,   & 
            delt,corrmax,jx,jy,jz,iterat,icvg,nx,ny,nz,time,AqueousToBulk)
          
          CALL SpeciesLocal(ncomp,nspec,jx,jy,jz)
          CALL totconc(ncomp,nspec,jx,jy,jz)
    
          IF (iterat > newtmax) THEN
            newtmax = iterat
          END IF

!!        (smr) seems useless
!!          if (iterat > 10) then
!!            continue
!!          end if

          IF (icvg == 1) THEN
            EXIT loop4001
          END IF

        END DO
      END DO
    END DO loop4001

!  3RD: CHECK FOR SUCCESS
    
    IF (icvg == 1) THEN

!! (smr) time step reduction should be done by the driver 
!! (smr) we only need to tell alquimia that this solve did not converge
!! (smr) and recover appropriate initial guesses

!!      (smr) leave out
!!      ddtold = delt
!!      delt = delt/10.0
!!      itskip = 1

      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            DO i = 1,ncomp+nspec
              sp(i,jx,jy,jz) = spold(i,jx,jy,jz)
              sp10(i,jx,jy,jz) = EXP(sp(i,jx,jy,jz))
            END DO
            DO ix = 1,nexchange
              spex(ix,jx,jy,jz) = spexold(ix,jx,jy,jz)
            END DO
            DO is = 1,nsurf+nsurf_sec
              spsurf(is,jx,jy,jz) = spsurfold(is,jx,jy,jz)
              spsurf10(is,jx,jy,jz) = EXP(spsurfold(is,jx,jy,jz))         
            END DO
            CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
          END DO
        END DO
      END DO

!!      CALL species(ncomp,nspec,nx,ny,nz)
!!      CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)

      IF (igamma == 2) THEN
        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
            IF (Duan) THEN
              CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
            ELSE
              CALL gamma(ncomp,nspec,jx,jy,jz)
            END IF
            END DO
          END DO
        END DO
      END IF

! (smr) this section is not needed
!!      WRITE(*,*)
!!      WRITE(*,*) '***** NO CONVERGENCE OF NEWTON ITERATIONS IN OS3D *****'
!!      WRITE(*,*) '                REDUCING TIME STEP'
!!      WRITE(*,5086) ddtold*OutputTimeScale
!!      WRITE(*,5085) delt*OutputTimeScale
!!      WRITE(*,*)
!!      GO TO 4000   ! Loop back to start the time step over
      
    ELSE
      iterat = newtmax
    END IF   !  End of IF block for case where no Newton convergence

! end -- straight from CrunchFlow's CrunchFlow.f90 subroutine, except where noted ----------------- 

! 4TH: UPDATE MINERAL VOLUME FRACTIONS, REACTIVE SURFACE AREA, AND POROSITY

  if (icvg == 1) then
  ! not converged, do nothing

    ! report alquimia status 
    status%error = kAlquimiaNoError
    status%converged = .false.
    status%num_rhs_evaluations = iterat
    status%num_jacobian_evaluations = iterat
    status%num_newton_iterations = iterat

  else 

    CALL mineral_update(nx,ny,nz,nrct,delt,dtnewest,ineg,jpor,deltmin)

    if (ineg == 1) then
    ! time step cut to adjust to complete mineral depletion
    ! delt = dtnewest

      write(*,*)'time step cut required'

      ! recover initial guesses
      DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          DO ik = 1,ncomp
            sp(ik,jx,jy,jz) = spold(ik,jx,jy,jz)
            sp10(ik,jx,jy,jz) = EXP(sp(ik,jx,jy,jz))
          END DO
          DO ix = 1,nexchange
            spex(ix,jx,jy,jz) = spexold(ix,jx,jy,jz)
          END DO
          DO is = 1,nsurf
            spsurf(is,jx,jy,jz) = spsurfold(is,jx,jy,jz)
          END DO
          CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
          IF (isaturate == 1) THEN
            CALL gases(ncomp,ngas,jx,jy,jz)
          END IF
        END DO
      END DO
      END DO

      ! (smr) leave out for now
      ! CALL species(ncomp,nspec,nx,ny,nz) 
      ! CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)

    ! --> (smr) to do
    ! --> give error msg: reduce time step to adjust to mineral dissolution
    ! --> suggest dtnewest
    
    else
    
      ! actual success
    
      ! report alquimia status
      status%error = kAlquimiaNoError
      status%converged = .true.
      status%num_rhs_evaluations = iterat
      status%num_jacobian_evaluations = iterat
      status%num_newton_iterations = iterat

    end if
  
  end if

  ! update old concentrations to be used in next time step as initial guesses
  ! straight from CrunchFlow.f90 (done there before transport step)
  DO jz = 1,nz
     DO jy = 1,ny
        DO jx = 1,nx
           CALL oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)  
           CALL oldkd(ncomp,jx,jy,jz)  ! smr
           IF (isaturate == 1) THEN
              CALL oldcongas(ncomp,ngas,jx,jy,jz)
           END IF
           CALL oldsurf(ncomp,nsurf,nsurf_sec,jx,jy,jz)
        END DO
     END DO
  END DO

  CALL xmass(nx,ny,nz,ncomp,nspec)
  xgramOld = xgram
  ! end straight from CrunchFlow.f90

  ! update density
  CALL density(jx,jy,jz)
  !!roOld = ro

  ! send state variables back to Alquimia state with solution of this solve
  call CopyAuxVarsToAlquimia(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             nretard, state, aux_data )

  !write(*,*)'------- GOING OUT ----------------------------------------'
  !call PrintState(state)
  !call PrintAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
  !                        nexchange, nsurf, ndecay, npot, &
  !                        aux_data)

  return

end subroutine ReactionStepOperatorSplit


! **************************************************************************** !
subroutine GetAuxiliaryOutput( &
     cf_engine_state, &
     properties, &
     state, &
     aux_data, &
     aux_output, &
     status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_double, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow
  use params, only: secyr, clg
  use concentration, only: sp10, sp, &
                           gam
  use mineral, only: dppt, &
                     silog

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaProperties), intent(in) :: properties
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaAuxiliaryOutputData), intent(inout) :: aux_output
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type(CrunchEngineState), pointer :: engine_state
  integer :: i, ikph
  real (c_double), pointer :: local_array(:)
  !!PetscReal :: porosity, volume
  logical, parameter :: copy_auxdata = .true.

  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1

  call c_f_pointer(cf_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !write (*, '(a)') "CrunchAlquimiaInterface::GetAuxiliaryOutput() :"

  ikph = engine_state%ikph
  ! FIXME(bja): this violates the no geochemistry calculations in alquimia rule.
  if (ikph /= 0) then
     aux_output%pH = -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz) )/clg
  else
     ! FIXME(bja, 2013-07) need a meaningful N/A value
     aux_output%pH = -100.d0
  end if

  !
  ! mineral data
  !
  call c_f_pointer(aux_output%mineral_reaction_rate%data, local_array, &
       (/aux_output%mineral_reaction_rate%size/))
    do i = 1, aux_output%mineral_reaction_rate%size
     local_array(i) = dppt(i,jx,jy,jz) / secyr
  end do

  call c_f_pointer(aux_output%mineral_saturation_index%data, local_array, &
       (/aux_output%mineral_saturation_index%size/))
  do i = 1, aux_output%mineral_saturation_index%size
  ! this assumes single cell problem, thus silog contains result from cell 1,1,1
  ! '1' hardwires this to first pathway (np) for i-th mineral, which is ok.
     local_array(i) = silog(1,i) 
  end do

  !
  ! primary data
  !
  call c_f_pointer(aux_output%primary_free_ion_concentration%data, local_array, &
       (/aux_output%primary_free_ion_concentration%size/))
  do i = 1, aux_output%primary_free_ion_concentration%size
     local_array(i) = sp10(i,jx,jy,jz)
  end do

  call c_f_pointer(aux_output%primary_activity_coeff%data, local_array, &
       (/aux_output%primary_activity_coeff%size/))
  do i = 1, aux_output%primary_activity_coeff%size
     local_array(i) = gam(i,jx,jy,jz)
  end do

  !
  ! secondary aqueous complex data
  !
  call c_f_pointer(aux_output%secondary_free_ion_concentration%data, local_array, &
       (/aux_output%secondary_free_ion_concentration%size/))
  do i = 1, aux_output%secondary_free_ion_concentration%size
     local_array(i) = sp10(engine_state%ncomp+i,jx,jy,jz)
  end do


  call c_f_pointer(aux_output%secondary_activity_coeff%data, local_array, &
       (/aux_output%secondary_activity_coeff%size/))
  do i = 1, aux_output%secondary_activity_coeff%size
     local_array(i) = gam(engine_state%ncomp+i,jx,jy,jz)
  end do

  status%error = kAlquimiaNoError
end subroutine GetAuxiliaryOutput


! **************************************************************************** !
subroutine GetProblemMetaData(cf_engine_state, meta_data, status)
  !  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_int, c_char, c_f_pointer

  use c_f_interface_module, only : f_c_string_ptr, f_c_string_chars

  use AlquimiaContainers_module

  use concentration, only: ulab, &
                           namsurf, &
                           namexc, &
                           distrib
  use mineral, only: umin

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaProblemMetaData), intent(out) :: meta_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type (c_ptr), pointer :: name_list(:)
  character (c_char), pointer :: name
  integer :: i, list_size, j
  integer(c_int), pointer :: idata(:)
  type(CrunchEngineState), pointer :: engine_state

  !write (*, '(a)') "Crunch_Alquimia_GetEngineMetaData() :"

  call c_f_pointer(cf_engine_state, engine_state)
  if (engine_state%integrity_check /= integrity_check_value) then
     status%error = kAlquimiaErrorEngineIntegrity
     call f_c_string_ptr("DEV_ERROR: pointer to engine state is not valid!", &
          status%message, kAlquimiaMaxStringLength)
     return
  end if

  !
  ! copy primary indices and names
  !

  write(*,*)'meta sizes: ',meta_data%primary_names%size 

  if (meta_data%primary_names%size /= engine_state%ncomp) then
     write (*, '(a, i3, a, i3, a)') "meta_data%primary_names%size (", &
          meta_data%primary_names%size, ") != crunchflow%ncomp(", &
          engine_state%ncomp, ")"
      stop
  end if
  list_size = meta_data%primary_names%size

  ! ulab : name of primary species

  call c_f_pointer(meta_data%primary_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(ulab(i)), &
          name, kAlquimiaMaxStringLength)     
  end do

!
! positivity constraints
!

  if (meta_data%positivity%size /= engine_state%ncomp) then
     write (*, '(a, i3, a, i3, a)') "meta_data%positivity%size (", &
          meta_data%positivity%size, ") != crunchflow%ncomp(", &
          engine_state%ncomp, ")"
      stop
  end if
  list_size = meta_data%positivity%size

  call c_f_pointer(meta_data%positivity%data, idata, (/list_size/))
  do i = 1, list_size
      if (i == engine_state%ikph) then
!       H+ component can be negative
        idata(i) = 0
      else
        idata(i) = 1 
      end if
  end do

  !
  ! copy mineral indices and names
  !

  if (meta_data%mineral_names%size /= engine_state%nrct) then
     write (*, '(a, i3, a, i3, a)') "meta_data%mineral_names%size (", &
          meta_data%mineral_names%size, ") != crunchflow%nrct(", &
          engine_state%nrct, ")"
  end if
  list_size = meta_data%mineral_names%size

  ! umin : name of minerals

  call c_f_pointer(meta_data%mineral_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(umin(i)), &
          name, kAlquimiaMaxStringLength)     
  end do

  !
  ! surface sites
  !
  if (meta_data%surface_site_names%size /= engine_state%nsurf) then
     write (*, '(a, i3, a, i3, a)') "meta_data%surface_site_names%size (", &
          meta_data%surface_site_names%size, ") != crunchflow%nsurf(", &
          engine_state%nsurf, ")"
  end if

  list_size = meta_data%surface_site_names%size

  ! namsurf : name of surface sites

  call c_f_pointer(meta_data%surface_site_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(namsurf(i)), &
          name, kAlquimiaMaxStringLength)     
  end do


  !
  ! ion exchange
  !
  if (meta_data%ion_exchange_names%size /= engine_state%nexchange) then
     write (*, '(a, i3, a, i3, a)') "meta_data%ion_exchange_names%size (", &
          meta_data%ion_exchange_names%size, ") != crunchflow%nexchange(", &
          engine_state%nexchange, ")"
  end if

  list_size = meta_data%ion_exchange_names%size

  ! namexc : name of exchange sites

  call c_f_pointer(meta_data%ion_exchange_names%data, name_list, (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars(trim(namexc(i)), &
          name, kAlquimiaMaxStringLength)     
  end do


  !
  ! isotherm indices
  !
  ! list_size = 0 for now
  ! CF does not track individual Kd species, can only tell by distrib
  !
  list_size = meta_data%isotherm_species_names%size
  call c_f_pointer(meta_data%isotherm_species_names%data, name_list, &
       (/list_size/))
  i = 0
  do j = 1, engine_state%ncomp
     
     if(distrib(j) > 0.0d0) then       
       i = i + 1
       if (i > engine_state%nretard) then
         write (*, '(a, i3, a, i3, a)') "meta_data%isotherm_species_names%size (", &
         meta_data%isotherm_species_names%size, ") != crunchflow%nretard(", &
         engine_state%nretard, ")"
       end if
       call c_f_pointer(name_list(i), name)
       call f_c_string_chars( &
            trim(ulab(j)), &  ! name of component, note index j
            name, kAlquimiaMaxStringLength)
     end if

  end do

  status%error = 0
end subroutine GetProblemMetaData


! **************************************************************************** !
!
! Private work routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine SetupCrunchOptions(input_filename) !, inputfilename)

  use, intrinsic :: iso_c_binding, only : c_char

  use c_f_interface_module, only : c_f_string_chars

  use crunchtype
  use params, only : mls

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  !type (option_type), intent(inout) :: option

  ! local variables
  CHARACTER (LEN=mls) :: inputfilename
  !PetscErrorCode :: ierr
  integer(i4b) :: ls
  
  ! set the crunchflow input file name
  call c_f_string_chars(input_filename, inputfilename)
  
  CALL stringlen(inputfilename,ls)
   
  !! needs to be rewritten for multiple processes
  open(666,file='PestControl.ant',status='unknown')
  write(666,'(a)')inputfilename(1:ls)
  close(666)

  !
  ! mpi
  !
  ! note: will need to let CF know that only rank==0 can print stuff
  !
  !option%global_comm = MPI_COMM_WORLD
  !call MPI_Comm_rank(MPI_COMM_WORLD, option%global_rank, ierr)
  !call MPI_Comm_size(MPI_COMM_WORLD, option%global_commsize, ierr)
  !call MPI_Comm_group(MPI_COMM_WORLD, option%global_group, ierr)
  !option%mycomm = option%global_comm
  !option%myrank = option%global_rank
  !option%mycommsize = option%global_commsize
  !option%mygroup = option%global_group

  !
  ! no access to output file name
  !

end subroutine SetupCrunchOptions

! **************************************************************************** !
subroutine SetEngineFunctionality(jpor,functionality)

  use AlquimiaContainers_module, only : AlquimiaEngineFunctionality

  use temperature, only: RunIsothermal
    
  implicit none

  ! function parameters
  integer(i4b), intent(in) :: jpor
  type (AlquimiaEngineFunctionality), intent(out) :: functionality

  functionality%thread_safe = .false. ! (smr) crunch global variables do not allow for thread safe engine
  if (RunIsothermal) then
    functionality%temperature_dependent = .false.
  else   
    functionality%temperature_dependent = .true.
  end if
  functionality%pressure_dependent = .false.
  if (jpor == -1) then
     functionality%porosity_update = .false.
  else 
     functionality%porosity_update = .true.
  end if
  functionality%operator_splitting = .true.
  functionality%global_implicit = .false.
  functionality%index_base = 1

end subroutine SetEngineFunctionality

! **************************************************************************** !
subroutine SetAlquimiaSizes(ncomp, nspec, nkin, nrct, ngas, &
                            nexchange, nsurf, ndecay, npot, &
                            nretard, ikin, sizes)

! note(smr): zeros denote unavailable processes for now, aug 8 2014

  use AlquimiaContainers_module, only : AlquimiaSizes
  
  use crunchtype
  use concentration, only: distrib

  implicit none

  ! function parameters
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot
  integer(i4b), intent(in) :: ikin
  integer(i4b), intent(inout) :: nretard
  type (AlquimiaSizes), intent(out) :: sizes

  integer(i4b)             :: i

  sizes%num_primary = ncomp
!!  if (nsurf > 0 .or. nexchange > 0 .or. nretard > 0) then
!!     sizes%num_sorbed = ncomp ! as per alquimia convention
!!  else
!!     sizes%num_sorbed = 0
!!  end if
  sizes%num_minerals = nrct ! nkin
  sizes%num_aqueous_complexes = nspec
  sizes%num_surface_sites = nsurf
  sizes%num_ion_exchange_sites = nexchange
  ! number of retardation species are not tracked explicitly in CF outside start98
  ! the array distrib (declared in concentration.f90) is allocated for ncomp
  ! with distrib(i) = kd (L/Kg) if species i is retarded, and distrib(i) = 0 if not 
  nretard = 0
  do i=1,ncomp
   if (distrib(i) > 0.0d0)  then
     nretard = nretard + 1     
   end if    
  end do
  sizes%num_isotherm_species = nretard
  if (nsurf > 0 .or. nexchange > 0 .or. nretard > 0) then
     sizes%num_sorbed = ncomp ! as per alquimia convention
  else
     sizes%num_sorbed = 0
  end if
  sizes%num_aqueous_kinetics = ikin
  call GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             sizes%num_aux_integers, sizes%num_aux_doubles)

  !call PrintSizes(sizes)

end subroutine SetAlquimiaSizes

! **************************************************************************** !
subroutine ProcessCrunchConstraint(engine_state, nco)

  use crunchtype
  use runtime, only: nchem, DensityModule, Duan
  use params, only : mls
  use concentration, only: ulab, sptmp10, sptmp, guess, &                            
                           spextmp10, spextmp, &
                           guess_surf, spsurftmp, spsurftmp10, &
                           gamtmp, MeanSalt, spgastmp10, spgastmp, &
                           spcond, spcond10, spcondgas, spcondgas10, &
                           scond, stmp, wtaq, &
                           spcondex, spcondex10, &
                           spcondsurf, spcondsurf10, &
                           LogPotentialInit, jinit, &
                           sp10, sp, s, sn, &
                           spgas10, spgas, &
                           vrSave, vrInitial, &
                           exchangesites, spex, spex10, totexch, &
                           spsurf10, spsurf, &
                           GasPressureTotal, GasPressureTotalInit, &
                           xgram, xgramOld, &
                           ctot, c_surf, itype, equilibrate
  use mineral, only: LogPotential_tmp, iedl, &
                     volfx, volin, &
                     VolumeLastTimeStep, &
                     area, areain, &
                     LogPotential !mintype
  use medium, only: AqueousToBulkCond, SaturationCond, constantpor, &
                    porin, por, porOld, porcond, isaturate                   
  use temperature, only: t, tempcond, rocond, ro
  use transport, only: activecell

  implicit none

  ! function parameters
  type(CrunchEngineState), intent(inout) :: engine_state  
  integer(i4b), intent(in)    :: nco
  
  ! local variables
  
  ! local pointer to crunchflow's state
  !
  ! geochemical sytem sizes, domain sizes
  !
  integer(i4b)           :: ncomp
  integer(i4b)           :: nspec
  integer(i4b)           :: nkin
  integer(i4b)           :: nrct
  integer(i4b)           :: ngas
  integer(i4b)           :: npot
  integer(i4b)           :: nexchange
  integer(i4b)           :: nexch_sec
  integer(i4b)           :: nsurf
  integer(i4b)           :: nsurf_sec
  integer(i4b)           :: ndecay
  integer(i4b)           :: ikin
  integer(i4b)           :: neqn
  integer(i4b)           :: nx
  integer(i4b)           :: ny
  integer(i4b)           :: nz
  !
  ! flags, pointer indexes, parameters
  !
  integer(i4b)           :: igamma ! gamma coeff update type
  integer(i4b)           :: ikph   ! points to pH
  integer(i4b)           :: jpor   ! porosity update type (0 is constant)
    
  ! actual local variables
  ! 
  integer(i4b)           :: i, ix, is, ik, nex  
  integer(i4b)           :: jx, jy, jz
  integer(i4b)           :: k, kk, npt, ns
  integer(i4b)           :: ipest    ! totally useless
  integer(i4b)           :: PestUnit ! totally useless

  real(dp)               :: tempc
  real(dp)               :: portemp
  real(dp)               :: denmol
  real(dp)               :: tk
  real(dp)               :: MeanSaltConcentration
  real(dp)               :: MassFraction
  real(dp)               :: sum
  real(dp)               :: convert

  character (len=mls)    :: namtemp
  
  logical                :: pest = .false.


! use pointers to get CF state, and be able to recycle code from CF as is

! geochemical system sizes
  ncomp = engine_state%ncomp
  nspec = engine_state%nspec  
  nkin = engine_state%nkin 
  nrct = engine_state%nrct 
  ngas = engine_state%ngas 
  npot = engine_state%npot 
  nexchange = engine_state%nexchange
  nexch_sec = engine_state%nexch_sec
  nsurf = engine_state%nsurf
  nsurf_sec = engine_state%nsurf_sec
  ndecay = engine_state%ndecay
  ikin = engine_state%ikin
  neqn = engine_state%neqn
  ! domain size (start98-alquimia sets it to nx=ny=nz=1)
  nx = engine_state%nx
  ny = engine_state%ny
  nz = engine_state%nz
  ! flags, pointer indexes, parameters
  igamma = engine_state%igamma
  ikph = engine_state%ikph
  jpor = engine_state%jpor

! allocate work variable
  ALLOCATE(AqueousToBulkCond(nchem))

! 1st: SPECIATE AND EQUILIBRATE CONDITION 'nco'

! this is straight from CrunchFlow start98, except where noted with (smr) -------------------------
  tempc = tempcond(nco)
  
  portemp = porcond(nco)

  CALL keqcalc2_init(ncomp,nrct,nspec,ngas,nsurf_sec,tempc)
  
  DO i = 1,ncomp
    namtemp = ulab(i)
    sptmp10(i) = guess(i,nco)
    sptmp(i) = DLOG(sptmp10(i))
  END DO
  DO ix = 1,nexchange
    spextmp10(ix) = 1.0
    spextmp(ix) = DLOG(spextmp10(ix))
  END DO
  DO is = 1,nsurf
    spsurftmp10(is) = guess_surf(is,nco)
    spsurftmp(is) = DLOG(spsurftmp10(is))
  END DO
  DO ik = 1,ncomp+nspec
    IF (ulab(ik) == 'H2O') THEN
      gamtmp(ik) = 0.00
    ELSE
      gamtmp(ik) = 0.0
    END IF
  END DO

  LogPotential_tmp = 0.0
  if (ngas > 0) then ! (smr) if
  spgastmp = -100.0d0
  spgastmp10 = 1.0D-35
  end if ! (smr) if

  CALL species_init(ncomp,nspec)
  CALL gases_init(ncomp,ngas,tempc)
  CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nchem)
  CALL exchange_init(ncomp,nspec,nexchange,nexch_sec,nchem)
  CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nco)
  CALL totgas_init(ncomp,nspec,ngas)
  CALL totsurf_init(ncomp,nsurf,nsurf_sec)
  CALL totexchange_init(ncomp,nexchange,nexch_sec,nco)
  
  ! leave this out (smr)
  !IF(.NOT. RunningPest) THEN
  !  WRITE(*,*)
  !  WRITE(*,*) ' --> Initializing Condition:  ',dumstring(1:ls)
  !  WRITE(*,*)
  !END IF
	
  CALL equilib_co2(ncomp,nspec,nrct,ngas,nsurf,igamma,ikph,  &
       nco,nexchange,nexch_sec,nsurf_sec,npot,neqn,tempc,portemp,    &
       DensityModule,ipest,PestUnit,pest)
 
!! Convert from bars pressure to n/V (mol/m^3) by converting to Pascals, then dividing by RT
  tk = tempc + 273.15d0
  denmol = 1.e05/(8.314*tk)   ! P/RT = n/V, with pressure converted from bars to Pascals
 
  if (ngas > 0) then !smr (if)
  spgastmp10 = spgastmp10*denmol
  spgastmp = DLOG(spgastmp10)
  end if !smr

  DO ik = 1,ncomp+nspec
    spcond(ik,nco) = sptmp(ik)
    spcond10(ik,nco) = sptmp10(ik)
  END DO
  DO kk = 1,ngas
    spcondgas(kk,nco) = spgastmp(kk)
    spcondgas10(kk,nco) = spgastmp10(kk)
  END DO
  
  DO i = 1,ncomp
    scond(i,nco) = stmp(i)
  END DO

  IF (DensityModule /= 'temperature') THEN
!   Calculate the correction for the mass fraction of water:  kg_solution/kg_water
    MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*scond(MeanSalt(1),nco) +   &
            wtaq(MeanSalt(2))*scond(MeanSalt(2),nco)) 
    MassFraction = 1.0/(1.0 + MeanSaltConcentration)
  ELSE
    MassFraction = 1.0
  END IF

  AqueousToBulkCond(nco) = rocond(nco)*SaturationCond(nco)*porcond(nco)*MassFraction

  DO nex = 1,nexchange+nexch_sec
    spcondex(nex,nco) = spextmp(nex) 
    spcondex10(nex,nco) = spextmp10(nex)
  END DO
  DO is = 1,nsurf
    spcondsurf(is,nco) = spsurftmp(is)
    spcondsurf10(is,nco) = spsurftmp10(is)
    IF (iedl(is) == 0) THEN     !!  Electrostatic option
       LogPotentialInit(is,nco) = LogPotential_tmp(is)
    END IF
  END DO
  DO ns = 1,nsurf_sec
    is = ns + nsurf
    spcondsurf(is,nco) = spsurftmp(is)
    spcondsurf10(is,nco) = spsurftmp10(is)
  END DO
! end: this is straight from CrunchFlow start98 ---------------------------------------------------


! 2nd: ALQUIMIA IS A SINGLE CELL MODEL - LINK CONDITION to cell 1 1 1 

! 'nco' points to the condition that applies to the alquimia cell
  jx=1
  jy=1
  jz=1
  activecell(jx,jy,jz) = 1 ! set to '1' to indicate to crunchflow it needs to do calculations with this cell
  jinit(jx,jy,jz) = nco ! jinit points to this condition, needed for variables that refer back to the condition
                        ! e.g. in the calculations during ReactionStepOperatorSplit, e.g. areain, volin
  
! 3rd: TRANSFER CONDITION SPECIATION RESULTS TO PERMANENT VARIABLES
! note: TO DO -- might have to overwrite some variables with alquimia state, e.g. volfx, area, 

! this is straight from CrunchFlow start98, except where noted with (smr) -------------------------
      DO ik = 1,ncomp+nspec
        sp10(ik,jx,jy,jz) = spcond10(ik,jinit(jx,jy,jz))
        sp(ik,jx,jy,jz) = spcond(ik,jinit(jx,jy,jz))
      END DO
      DO i = 1,ncomp
        s(i,jx,jy,jz) = scond(i,jinit(jx,jy,jz))
        sn(i,jx,jy,jz) = scond(i,jinit(jx,jy,jz))
      END DO

!!      DO ix = 1,nexchange+nexch_sec
!!        spex10(ix,jx,jy,jz) = spcondex10(ix,jinit(jx,jy,jz))
!!        spex(ix,jx,jy,jz) = spcondex(ix,jinit(jx,jy,jz))
!!      END DO
!!      DO is = 1,nsurf
!!        spsurf(is,jx,jy,jz) = spcondsurf(is,jinit(jx,jy,jz))
!!      END DO
!!      DO is = 1,nsurf+nsurf_sec
!!        spsurf10(is,jx,jy,jz) = spcondsurf10(is,jinit(jx,jy,jz))
!!      END DO

      DO kk = 1,ngas
        spgas10(kk,jx,jy,jz) = spcondgas10(kk,jinit(jx,jy,jz))
        spgas(kk,jx,jy,jz) = spcondgas(kk,jinit(jx,jy,jz))
      END DO
      sum = 0.0
      DO k = 1,nrct
        volfx(k,jx,jy,jz) = volin(k,jinit(jx,jy,jz))
        VolumeLastTimeStep(k,jx,jy,jz) = volfx(k,jx,jy,jz)
        area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))
  !      if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz)
        sum = sum + volfx(k,jx,jy,jz)
      END DO

      IF (Duan) THEN
!!      Save Residual Volume for CO2 fugacity calculation
        vrSave(jx,jy,jz) = vrInitial(jinit(jx,jy,jz))
      END IF

      IF (constantpor /= 0.0) THEN
         porin(jx,jy,jz) = constantpor
         por(jx,jy,jz) = constantpor
         porOld(jx,jy,jz) = por(jx,jy,jz)
      END IF

!(smr)      IF (jpor == 2) THEN                  !! Read porosity from file
!(smr)        porin(jx,jy,jz) = work3(jx,jy,jz)
!(smr)      ELSE
        porin(jx,jy,jz) = porcond(jinit(jx,jy,jz))
!(smr)      END IF
!!      satliq(jx,jy,jz) = SaturationCond(jinit(jx,jy,jz))

      por(jx,jy,jz) = porin(jx,jy,jz)
      porOld(jx,jy,jz) = por(jx,jy,jz)
      IF (por(jx,jy,jz) <= 0.0) THEN
        WRITE(*,*)
        WRITE(*,*) '  You have specified a porosity < 0'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

!!CIS      IF (por(jx,jy,jz) > 1.0) THEN
!!CIS        WRITE(*,*)
!!CIS        WRITE(*,*) '  You have specified a porosity > 1'
!!CIS        WRITE(*,*)
!!CIS        READ(*,*)
!!CIS        STOP
!!CIS      END IF

      t(jx,jy,jz) = tempc
      CALL density(jx,jy,jz)

!  Change the site concentrations to sites per bulk volume porous medium

      IF (DensityModule /= 'temperature') THEN
!       Calculate the correction for the mass fraction of water:  kg_solution/kg_water
        MeanSaltConcentration = 0.001d0*(wtaq(MeanSalt(1))*s(MeanSalt(1),jx,jy,jz) +   &
            wtaq(MeanSalt(2))*s(MeanSalt(2),jx,jy,jz)) 
        MassFraction = 1.0d0/(1.0d0 + MeanSaltConcentration)
      ELSE
        MassFraction = 1.0d0
      END IF

!!      convert = ro(jx,jy,jz)*porin(jx,jy,jz)*MassFraction
      convert = ro(jx,jy,jz)*porcond(jinit(jx,jy,jz))*SaturationCond(jinit(jx,jy,jz))*MassFraction

      DO ix = 1,nexchange
        exchangesites(ix,jx,jy,jz) = convert*totexch(ix,jinit(jx,jy,jz)) ! Now in equivalents/m3 por. med.
!!        exchangesites(ix,jx,jy,jz) = totexch(ix,jinit(jx,jy,jz)) ! Already in equivalents/m3 por. med.
      END DO
      
      do ix = 1,nexchange
        spex(ix,jx,jy,jz) = spcondex(ix,jinit(jx,jy,jz))
        spex10(ix,jx,jy,jz) = convert*spcondex10(ix,jinit(jx,jy,jz))  ! Now in eq/m3 por. med.
      end do
      DO ix = 1,nexch_sec
        spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,jinit(jx,jy,jz))  ! Now in eq/m3 por. med.
      END DO

      DO is = 1,nsurf+nsurf_sec
        spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,jinit(jx,jy,jz))
      END DO
      DO is = 1,nsurf
        spsurf(is,jx,jy,jz) = LOG(convert*spcondsurf10(is,jinit(jx,jy,jz)))
      END DO

!     also temperature
      t(jx,jy,jz) = tempcond(jinit(jx,jy,jz))

! end: this is straight from CrunchFlow start98 ---------------------------------------------------

! some more initializations, from CrunchFlow.f90 this time
  !DO jz = 1,nz; DO jy = 1,ny; DO jx = 1,nx
        DO npt = 1,npot
          LogPotential(npt,jx,jy,jz) = LogPotentialInit(npt,jinit(jx,jy,jz))
        END DO
  !END DO; END DO; END DO

  IF (Duan) THEN
    !DO jz = 1,nz; DO jy = 1,ny; DO jx = 1,nx
          GasPressureTotal(jx,jy,jz) = GasPressureTotalInit(jinit(jx,jy,jz))
    !END DO; END DO; END DO
  END IF
! end: from CrunchFlow.f90

! initialize old concentrations
  call exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
  CALL oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)  
  CALL oldkd(ncomp,jx,jy,jz)  
  IF (isaturate == 1) THEN
    CALL oldcongas(ncomp,ngas,jx,jy,jz)
  END IF
  CALL oldsurf(ncomp,nsurf,nsurf_sec,jx,jy,jz)

! intialize xmass old for this condition
  CALL xmassNodeByNode(jx,jy,jz,ncomp,nspec)
  xgramOld = xgram

! deallocate
  IF (ALLOCATED(AqueousToBulkCond)) THEN
    DEALLOCATE(AqueousToBulkCond)
  END IF

  return

end subroutine ProcessCrunchConstraint


! **************************************************************************** !
subroutine ConvertAlquimiaConditionToCrunch(state,properties, &
                                            ncomp,nspec,nrct,nkin, &
                                            ngas,nsurf,nexchange, &
                                            ikin,nexch_sec,nsurf_sec, &
                                            alquimia_condition)
  use, intrinsic :: iso_c_binding

  use c_f_interface_module, only : c_f_string_ptr, CaseInsensitiveStrcmp
  use AlquimiaContainers_module

  use runtime, only: nchem, DensityModule
  use concentration, only: ncon, condlabel, ctot, itype, guess, equilibrate, &
                           gaspp, OneOverMassFraction, conversion, icec, cec, &
                           totexch, c_surf, guess_surf, ulab, MeanSalt, iexchange, &
                           kexch, ksurf, wtaq, SolidDensity
  use temperature, only: tempcond, rocond
  use medium, only: porcond, SaturationCond
  use mineral, only: volin, areain, iarea, specific, volmol, wtmin, &
                     site_density,umin
  use params, only: mls, mchem

  implicit none
 
! function parameters
  type (AlquimiaState), intent(in)                             :: state
  type (AlquimiaProperties), intent(in)                     ::properties
  type (AlquimiaGeochemicalCondition), intent(in) :: alquimia_condition

  integer(i4b), intent(in)          :: ncomp
  integer(i4b), intent(in)          :: nspec
  integer(i4b), intent(in)          :: nrct
  integer(i4b), intent(in)          :: nkin
  integer(i4b), intent(in)          :: ngas
  integer(i4b), intent(in)          :: nsurf
  integer(i4b), intent(in)          :: nexchange
  integer(i4b), intent(in)          :: ikin
  integer(i4b), intent(in)          :: nexch_sec
  integer(i4b), intent(in)          :: nsurf_sec

! local variables
  integer :: i, ix, j, k, ks
  real (c_double), pointer :: data(:)
  real (dp) :: RoSolution
  real(dp) :: permole, MeanSaltConcentration, MineralMolality
  real(dp) :: tempc, rotemp
  character (kAlquimiaMaxStringLength) :: name, constraint_type
  character (kAlquimiaMaxStringLength) :: associated_species
  character (len=mls) :: dumstring
  type (AlquimiaAqueousConstraint), pointer :: alq_aqueous_constraints(:)
  type (AlquimiaMineralConstraint), pointer :: alq_mineral_constraints(:)

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: unitsflag
  REAL(DP), DIMENSION(:), ALLOCATABLE                           :: pH
  REAL(DP), DIMENSION(:), ALLOCATABLE                           :: guessph
!
! there are internal 
!
  IF (ALLOCATED(pH)) THEN
    DEALLOCATE(pH)
    ALLOCATE(pH(mchem))
  ELSE
    ALLOCATE(pH(mchem))
  ENDIF
  IF (ALLOCATED(guesspH)) THEN
    DEALLOCATE(guesspH)
    ALLOCATE(guesspH(mchem))
  ELSE
    ALLOCATE(guesspH(mchem))
  ENDIF
  IF (ALLOCATED(unitsflag)) THEN
    DEALLOCATE(unitsflag)
    ALLOCATE(unitsflag(mchem))
  ELSE
    ALLOCATE(unitsflag(mchem))
  END IF
  unitsflag = 1
!
! update number of solutions and point to it
!
  nchem = nchem +1 
!
! check whether there is enough space; there should be, but ...
!
  if (nchem > mchem) then
    write(*,*) 'Number of solutions exceeds maximum number of allocated spaces: ',mchem
    stop
  end if
! 
! reallocate to make room for new nchem
!
call reallocate(ncomp,nspec,nrct,nkin,ngas,nsurf,nexchange,ikin,nexch_sec,nsurf_sec)
!
! assign name of condition, prepend 'alq_' to indicate this comes from driver
!
  call c_f_string_ptr(alquimia_condition%name, name)
  condlabel(nchem) = 'alq_'//trim(adjustl(name))
!
! aqueous species constraints
!
  call c_f_pointer(alquimia_condition%aqueous_constraints%data, &
       alq_aqueous_constraints, (/alquimia_condition%aqueous_constraints%size/))
!
! check size
!
  if (alquimia_condition%aqueous_constraints%size /= ncomp) then
     write(*,*) 'Number of aqueous constraints ' // &
                     'does not equal the number of primary chemical ' // &
                     'components in constraint: ' // &
                      trim(condlabel(nchem))
     stop
  end if
!
! get temperature, density, porosity , saturation from AlquimiaState & Properties
!
  tempcond(nchem) = state%temperature
  rocond(nchem)      = state%water_density
  porcond(nchem)    = state%porosity  
  SaturationCond(nchem) = properties%saturation
!!  SolidDensity(nchem) = properties%solid_density
  SolidDensity(nchem) = 2650.0d0 ! hard wired until I figure the problem with memory bug
                                                         ! created by adding solid_density --> does it have to do with batch_chem????
!
! transfer constraints to crunchflow format
!
  do_aqueous_components: do i = 1, ncomp

     do_alquimia_constraints: do j = 1, alquimia_condition%aqueous_constraints%size

        call c_f_string_ptr(alq_aqueous_constraints(j)%primary_species_name, name)
        if (name == ulab(i)) then
           ! j = alquimia order, i = crunchflow order
           exit do_alquimia_constraints
        end if

     end do do_alquimia_constraints
!
!    get value
!
     ctot(i,nchem) = alq_aqueous_constraints(j)%value
!
!    defaults
!
     if (ctot(i,nchem) == 0.0) then
       ctot(i,nchem) = 1.e-30
     end if
     itype(i,nchem) = 1
     guess(i,nchem) = 1.E-06
     equilibrate(i,nchem) = .false.
     ph(nchem) = 7.0
     unitsflag(nchem) = 5   ! molarity
!
!    get constraint type and associated name if necessary
!
     call c_f_string_ptr(alq_aqueous_constraints(i)%constraint_type, constraint_type)
     call c_f_string_ptr(alq_aqueous_constraints(i)%associated_species, &
          associated_species)
!
!    assign type of constraint, and name associated with it if mineral or gas constraint
!
     if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringFree)) then
        itype(i,nchem) = 8
        ctot(i,nchem)    = alq_aqueous_constraints(j)%value
        guess(i,nchem) = alq_aqueous_constraints(j)%value
     else if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringTotalAqueous)) then
        itype(i,nchem) = 1
        guess(i,nchem) = alq_aqueous_constraints(j)%value
        ctot(i,nchem)    = alq_aqueous_constraints(j)%value

        !  as per CrunchFlow
        if_Hplus: if (ulab(i) == 'h+' .or.  ulab(i) == 'H+') then
          if (guess(i,nchem) < 1.e-15) then
            if (ctot(i,nchem) > 0.0) then
              guess(i,nchem) = 1.e-05
            else
              guess(i,nchem) = 1.e-08
            end if
          end if
        end if if_Hplus

     else if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringTotalSorbed)) then
        itype(i,nchem) = 1
!        equilibrate(i,nchem) = .true.
     else if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringPH)) then
        itype(i,nchem) = 7
        ph(nchem) = alq_aqueous_constraints(j)%value
        guess(i,nchem) = 10**(-ph(nchem))
     else if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringMineral)) then
        itype(i,nchem) = 3
        ncon(i,nchem) = trim(associated_species)
     else if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringGas)) then
        itype(i,nchem) = 4
        ncon(i,nchem) = trim(associated_species)
        gaspp(i,nchem) = alq_aqueous_constraints(j)%value
        if (gaspp(i,nchem) == 0.0) then
          gaspp(i,nchem) = 1.e-30
         end if
         ctot(i,nchem) = gaspp(i,nchem)
     else if (CaseInsensitiveStrcmp(constraint_type, kAlquimiaStringCharge)) then
        itype(i,nchem) = 2
     else
        write(*,*)'Constraint type: ' // trim(constraint_type) // &
                  ' not recognized in constraint,concentration'
     end if
  end do  do_aqueous_components
! 
! recompute units -- alquimia gives molarity 
! this comes straight from CrunchFlow find_condition
!
  tempc = tempcond(nchem)
  rotemp = rocond(nchem) 
  CALL rocalc(tempc,rotemp,nchem,DensityModule,unitsflag)
  rocond(nchem) = rotemp

  IF (unitsflag(nchem) == 5) THEN              !  Units in molarity 

    IF (DensityModule /= 'temperature') THEN
!     Calculate the correction for the mass fraction of water:  kg_solution/kg_water
      MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*ctot(MeanSalt(1),nchem) +     &
      wtaq(MeanSalt(2))*ctot(MeanSalt(2),nchem))                                       !! Mean salt in kg salt/L solution (g/cm**3)
      RoSolution = rocond(nchem)/1000.0d0                                              !! Convert fluid density from kg/m**3 to kg/L (or g/cm**3), same as salt
      OneOverMassFraction(nchem) = RoSolution/(RoSolution - MeanSaltConcentration)
      conversion(nchem) = OneOverMassFraction(nchem)/RoSolution
    ELSE
      OneOverMassFraction(nchem) = 1.0d0
      conversion(nchem) = 1000.0d0/rocond(nchem)
    END IF

    DO i = 1,ncomp
      if (itype(i,nchem) /= 7 .and. &
           itype(i,nchem) /= 4 )  then
        ctot(i,nchem) = ctot(i,nchem)*conversion(nchem)               
        guess(i,nchem) = guess(i,nchem)*conversion(nchem)
      else 
         continue
      end if
    END DO

  ELSE 
    conversion(nchem) = 1.0d0
    OneOverMassFraction(nchem) = 1.0d0
  END IF

!  Change "unitsflag" back to molality now and recompute the density
  IF (unitsflag(nchem) == 5) THEN
    unitsflag(nchem) = 1
  END IF

  CALL rocalc(tempc,rotemp,nchem,DensityModule,unitsflag)
  rocond(nchem) = rotemp
! 
! end of copied section, unit conversion
!

! 
! minerals
!
!
!  we should do it like this but too many changes for now
!  call c_f_pointer(state%mineral_volume_fraction%data, data, (/nrct/))
!  do i = 1, nrct
!     volin(i,nchem) = data(i)
!  end do

!  call c_f_pointer(state%mineral_specific_surface_area%data, data, (/nrct/))
!  do i = 1, nrct
!     areain(i,nchem) = data(i) ! m2/m3 bulk
!     iarea(i,nchem)   = 0  ! surface area provided as bulk

!     if (volin(i,nchem) /= 0.0) then
!        specific(i,nchem) = areain(i,nchem)*volmol(i)/(volin(i,nchem)*wtmin(i))
!     else
!        specific(i,nchem) = 0.0
!     end if 
!  end do

!
! minerals
!
! this will have to go if we change the original approach - but for now we'll keep it
  if (alquimia_condition%mineral_constraints%size > 0 .and. &
       alquimia_condition%mineral_constraints%size /= nrct) then
    write(*,*)'number of minerals in constraint does not match number of minerals'
    stop
  end if

  call c_f_pointer(alquimia_condition%mineral_constraints%data, &
       alq_mineral_constraints, (/alquimia_condition%mineral_constraints%size/))
  do i = 1, nrct
     volin(i,nchem) = alq_mineral_constraints(i)%volume_fraction
  end do

  do i = 1, nrct
     areain(i,nchem) = alq_mineral_constraints(i)%specific_surface_area  ! m2/m3 bulk
     iarea(i,nchem)   = 0  ! surface area provided as bulk

     if (volin(i,nchem) /= 0.0) then
        specific(i,nchem) = areain(i,nchem)*volmol(i)/(volin(i,nchem)*wtmin(i))
     else
        specific(i,nchem) = 0.0
     end if 
  end do
 !
 ! ion exchange
 !
  call c_f_pointer(state%cation_exchange_capacity%data, data, (/nexchange/))
  do ix = 1, nexchange
!!     exchangesites(ix,jx,jy,jz) = data(ix)

     IF (iexchange(ix) == 0) THEN                                          !!  Bulk exchange on the sediment

!       not allowed for now, there is no SolidSolutionRatio concept in Alquimia
        write(*,*) 'The CrunchFlow Alquimia interface does not  ' // &
                        'support bulk exchange on sediment for now. ' // &
                        'Please specify a mineral to exchange '
        stop

!!      IF (icec(ix) == 1) THEN
!!        totexch(ix,nchem) = cec(ix,nchem)*SolidSolutionRatio(nchem)
!!      ELSE
!!        IF (SolidSolutionRatio(nchem) == 0.0d0) THEN
!!          cec(ix,nchem) = 0.0d0
!!        ELSE
!!          cec(ix,nchem) = totexch(ix,nchem)/SolidSolutionRatio(nchem)
!!        END IF
!!      END IF

    ELSE IF (iexchange(ix) == 1) THEN   !  Exchange on a specific mineral

        k = kexch(ix)
!!      IF (icec(ix) == 1) THEN                                             !!  Calculate from equivalents/g mineral and the mineral volume fraction
        icec(ix) = 1
        cec(ix,nchem) = data(ix) * volmol(k) / wtmin(k) /  OneOverMassFraction(nchem)
        totexch(ix,nchem) = data(ix)*volin(k,nchem)/(SaturationCond(nchem)*porcond(nchem)*rocond(nchem))  
!!        totexch(ix,nchem) = OneOverMassFraction(nchem)*cec(ix,n<chem)*wtmin(k)*volin(k,nchem)/(volmol(k)*SaturationCond(nchem)*porcond(nchem)*rocond(nchem))  
!!      ELSE IF (icec(ix) == 0) THEN                                         !!  Direct specification of totexch (equivalents/kgw)--need to calculate CEC for later use if mineral fraction changes
!!        cec(ix,nchem) = totexch(ix,nchem)*volmol(k)*SaturationCond(nchem)*porcond(nchem)*rocond(nchem)/(wtmin(k)*volin(k,nchem)*OneOverMassFraction(nchem))
!!      END IF 

    END IF 

  end do
  !
  ! equilibrium surface complexation
  !
  call c_f_pointer(state%surface_site_density%data, data, (/nsurf/))

  DO ks = 1,nsurf
    k = ksurf(ks)
    IF (specific(k,nchem) == 0.0) THEN
      WRITE(*,*) 
      WRITE(*,*) ' Specific surface area for mineral = 0 serving as sorbate for surface complex'
      dumstring = umin(k)
      WRITE(*,*) ' --->Mineral: ', dumstring
      WRITE(*,*) ' --->In geochemical condition: ', condlabel
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    IF (volin(k,nchem) == 0.0) THEN
      WRITE(*,*) 
      WRITE(*,*) ' Specific surface areaVolume fraction for mineral = 0 serving as sorbate for surface complex'
      dumstring = umin(k)
      WRITE(*,*) ' --->Mineral: ', dumstring
      WRITE(*,*) ' --->In geochemical condition: ', condlabel
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

!   Site_density Crunch: moles site/m^2 mineral site_density 
!   Site density Alquimia: moles/m^3 bulk data
!   Specific:     m^2/g mineral
!   Wtmin:        g/mole mineral

    site_density(ks,nchem) = data(ks) / areain(k,nchem)
    permole = site_density(ks,nchem)*specific(k,nchem)*wtmin(k)    !  Mole sites/Mole mineral

!   Now convert to moles sites per kg solution
!!  volin(m^3 mineral/m^3 porous medium) /( volmol[m^3/mol] * rocond[kg/m^3 fluid] * porcond[m^3 pore/m^3 PM] * SaturationCond[m^3 fluid/m^3 pore)

!!  Moles mineral/kgw
!!    IF (MineralMoles(k,nchem) > 0.0d0 .AND. volin(k,nchem) == 0.0d0) THEN
!!      c_surf(ks,nchem) = permole*MineralMoles(k,nchem)
!!    ELSE IF (volin(k,nchem) == 0.0d0) THEN
!!      MineralMolality = voltemp(k,nchem)/( volmol(k)*rocond(nchem)*porcond(nchem)*SaturationCond(nchem) )
!!      c_surf(ks,nchem) = permole*MineralMolality
!!    ELSE
      MineralMolality = volin(k,nchem)/( volmol(k)*rocond(nchem)*porcond(nchem)*SaturationCond(nchem) )
      c_surf(ks,nchem) = permole*MineralMolality
!!    END IF

    IF (c_surf(ks,nchem) < 1.D-30) THEN
      c_surf(ks,nchem) = 1.D-30
    END IF

!    WRITE(*,*) ' k = ',k
!    WRITE(*,*) ' Site density ',site_density(ks,nchem)
!    WRITE(*,*) ' Specific area ',specific(ks,nchem)
!    WRITE(*,*) ' Wtmin = ',wtmin(k)
!    WRITE(*,*) ' Permole = ',permole
!    WRITE(*,*) ' Volin = ',volin(k,nchem)
!    WRITE(*,*) ' volmol = ',volmol(k)
  END DO
  
  DO ks = 1,nsurf
    guess_surf(ks,nchem) = c_surf(ks,nchem)
  END DO
!
! deallocate internally allocated vars
! 
  DEALLOCATE(pH)
  DEALLOCATE(guesspH)
  DEALLOCATE(unitsflag)
  return

end subroutine ConvertAlquimiaConditionToCrunch


! **************************************************************************** !
subroutine CopyAlquimiaToAuxVars(copy_auxdata, hands_off, &
                                 state, aux_data, prop, &
                                 ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot, &
                                 nretard )

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module

  ! crunchflow
  use crunchtype
  use temperature, only : t, ro
  use transport, only : satliq
  use medium, only : por
  use concentration, only : sn, &
!                            ssurfold, &
!                            sexold, &
                            exchangesites, &
                            ssurfn, &
                            distrib, &
                            SolidDensity, &
                            jinit, &
                            xgram, &
                            ratek

  use mineral, only : volfx, area, &
                                  rate0
  use params, only: secyr

  implicit none

  ! function parameters
  logical, intent(in) :: copy_auxdata
  logical, intent(in) :: hands_off
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaProperties), intent(in) :: prop
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot
  integer(i4b), intent(in) :: nretard

  ! local variables
  real (c_double), pointer :: data(:)
  integer :: i, iret
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1
  
  !write (*, '(a)') "Crunch_Alquimia_CopyAlquimiaToAuxVars() :"

  !
  ! copy state to crunchflow
  !
  ro(jx,jy,jz)            = state%water_density
  satliq(jx,jy,jz)        = prop%saturation
  t(jx,jy,jz)             = state%temperature
  !! not needed by crunch = state%aqueous_pressure
  por(jx,jy,jz)           = state%porosity  
  !! not needed by crunch = prop%volume
  !!write(*,*)'HERE dens, temp, por: ', state%water_density, state%temperature, state%porosity
  !
  ! primary aqueous
  !
  call c_f_pointer(state%total_mobile%data, data, (/ncomp/))
  do i = 1, ncomp
     sn(i,jx,jy,jz) = data(i) / ro(jx,jy,jz) / xgram(jx,jy,jz) * 1.d3  ! fixed point for nonlinear problem is in 'sn'
     !!write(*,*)'receive sn(i,jx,jy,jz): ',sn(i,jx,jy,jz),i
  end do

  ! sorbed primary
  if (nsurf > 0 .or. nexchange > 0 .or. nretard > 0) then
     call c_f_pointer(state%total_immobile%data, data, (/ncomp/))
     do i = 1, ncomp
        ! don't send anything back to crunch, rely on data stored on global vars for now
        ! crunchflow keeps these data separate so there is not way to deal with
        ! both exchange and surface simultaneously through total_immobile
        !!ssurfold(i,jx,jy,jz) = data(i)                                          
        !!sexold(i,jx,jy,jz)   = data(i)
        !!write(*,*)'receive from alquimia i, sexold(i,jx,jy,jz) ', i, sexold(i,jx,jy,jz)  
        !!write(*,*)'receive from alquimia i, ssurfold(i,jx,jy,jz) ', i, ssurfold(i,jx,jy,jz)  
     end do
  end if

  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction%data, data, (/nrct/))
  do i = 1, nkin
     volfx(i,jx,jy,jz) = data(i)
  end do

  call c_f_pointer(state%mineral_specific_surface_area%data, data, (/nrct/))
  do i = 1, nkin
     area(i,jx,jy,jz) = data(i)
  end do

  !
  ! ion exchange
  !
  call c_f_pointer(state%cation_exchange_capacity%data, data, (/nexchange/))
  do i = 1, nexchange
     exchangesites(i,jx,jy,jz) = data(i)
  end do

  !
  ! equilibrium surface complexation
  !
  call c_f_pointer(state%surface_site_density%data, data, (/nsurf/))
  do i = 1, nsurf
     ssurfn(i,jx,jy,jz) = data(i) ! CF will not use this, will use LogTotalSurface (in auxiliary for consistency)
  end do

  !
  ! in hands-off mode geochemical properties
  ! are not provided by the driver so copying them over would
  ! lose crunchflow's input file values 
  !
  if_hands_off: if (.not. hands_off) then
  
  !
  ! isotherms (smr) only linear kd model - need to convert units to L/Kg solid
  !
  call c_f_pointer(prop%isotherm_kd%data, data, (/nretard/))
  iret = 0
  do i = 1, ncomp
              
     if (distrib(i) > 0.0d0) then
       iret = iret + 1
!       distrib(i) = data(iret) * por(jx,jy,jz) / (1.0d0 - por(jx,jy,jz)) &
!                               / SolidDensity(jinit(jx,jy,jz)) / ro(jx,jy,jz)  * 1.d3       
       if (por(jx,jy,jz) == 1.0d0) then
          continue
       else
          distrib(i) = data(iret) / (1.0d0 - por(jx,jy,jz)) &
                          / SolidDensity(jinit(jx,jy,jz)) / ro(jx,jy,jz)  * 1.d3       
       end if  
     end if

  end do
  if (iret /= nretard) then
    write(*,*)'number of isotherm species is not that in crunchflows input file'
    write(*,*)'current number: ',iret,' input file number: ',nretard
    stop
  end if

  !
  ! mineral reaction rate constant
  !
  call c_f_pointer(prop%mineral_rate_cnst%data, data, &
       (/prop%mineral_rate_cnst%size/))
  do i = 1, prop%mineral_rate_cnst%size
     ! rate0 is in units of mol/m2/yr
       rate0(1,i) = data(i) *secyr         ! 1 --> hardwired for now to 1 pathway
  end do

  !
  ! aqueous kinetic reaction rate constant
  !
  call c_f_pointer(prop%aqueous_kinetic_rate_cnst%data, data, &
       (/prop%aqueous_kinetic_rate_cnst%size/))
  do i = 1, prop%aqueous_kinetic_rate_cnst%size
     ! ratek is in units of 1/yr 
       ratek(1,i) = data(i) * secyr      ! 1 -->hardwired for now to 1 pathway
  end do

  end if if_hands_off
  
  if (copy_auxdata) then
     call UnpackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                      nexchange, nsurf, ndecay, npot, &
                                      aux_data)
  end if

end subroutine CopyAlquimiaToAuxVars

! **************************************************************************** !
subroutine CopyAuxVarsToAlquimia(ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot, &
                                 nretard, state, aux_data)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module

  ! crunch
  use crunchtype
  use temperature, only : t, ro
  use transport, only : satliq
  use medium, only : por !!, porin
  use concentration, only : sn, &
                             ssurfold, &
                             sexold, &
                             skdold, &
                             exchangesites, &
                             ssurfn, &
                             jinit, &
                             SolidDensity, &
                             distrib, &
                             xgram

  use mineral, only : volfx, area

  implicit none

  ! function parameters
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot
  integer(i4b), intent(in) :: nretard
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  real (c_double), pointer :: data(:)
  integer :: i
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1
  real(dp)                :: retardation

  !write (*, '(a)') "Crunch_Alquimia_CopyAuxVarsToAlquimia() :"

  !
  ! state
  !
  state%water_density = ro(jx,jy,jz)
  state%temperature = t(jx,jy,jz)
  !!state%aqueous_pressure = crunch has nothing for you
  state%porosity = por(jx,jy,jz)
  !!write(*,*)'dens, temp, por: ', state%water_density, state%temperature, state%porosity
  !
  ! primary aqueous species
  !
  call c_f_pointer(state%total_mobile%data, data, (/ncomp/))
  do i = 1, ncomp
     data(i) = sn(i,jx,jy,jz) * ro(jx,jy,jz) * xgram(jx,jy,jz) / 1.d3
     !!write(*,*)'sn(i,jx,jy,jz): ',sn(i,jx,jy,jz)
  end do  
  !
  ! sorbed
  !
  if (nsurf > 0 .or. nexchange > 0 .or. nretard > 0) then
     call c_f_pointer(state%total_immobile%data, data, (/ncomp/))
     do i = 1, ncomp
        ! send it to alquimia for plotting etc but cannot make it back through total_immobile
        data(i) = 0.0d0
        if (nexchange > 0) then
          data(i) = data(i) + sexold(i,jx,jy,jz)
        end if
        if (nsurf > 0) then
          data(i) = data(i) + ssurfold(i,jx,jy,jz)
        end if
        if (nretard > 0) then
          data(i) = data(i) + skdold(i,jx,jy,jz)
!          Retardation = 0.001d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))/por(jx,jy,jz)
!          data(i) = data(i) + xgram(jx,jy,jz) * satliq(jx,jy,jz) * ro(jx,jy,jz) * sn(i,jx,jy,jz) * Retardation * distrib(i)           
        end if
        !! write(*,*)'i, ssurfold(i,jx,jy,jz) ', i, sexold(i,jx,jy,jz)    
     end do
  end if
  
  !
  ! minerals
  !
  call c_f_pointer(state%mineral_volume_fraction%data, data, (/nrct/))
  do i = 1, nrct
     data(i) = volfx(i,jx,jy,jz)
  end do

  call c_f_pointer(state%mineral_specific_surface_area%data, data, (/nrct/))
  do i = 1, nkin
     data(i) = area(i,jx,jy,jz)
  end do

  !
  ! ion exchange
  !
  call c_f_pointer(state%cation_exchange_capacity%data, data, (/nexchange/))
  do i = 1, nexchange
     data(i) = exchangesites(i,jx,jy,jz)
  end do

  !
  ! equilibrium surface complexation
  !
  call c_f_pointer(state%surface_site_density%data, data, (/nsurf/))
  do i = 1, nsurf
     data(i) = ssurfn(i,jx,jy,jz) 
     !! write(*,*)'ssurfn(i,jx,jy,jz): ',ssurfn(i,jx,jy,jz)
  end do

  ! NOTE(bja): isotherms are material properties, and can't be changed
  ! by chemistry. We don't need to copy theme here!

  ! NOTE(smr): reaction rates constants are properties, and can't be changed
  ! by chemistry. We don't need to copy theme here!

  call PackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot, &
                                 aux_data)

end subroutine CopyAuxVarsToAlquimia

! **************************************************************************** !
!
! Auxiliary Data helper routines
!
! We are packaging data into the auxiliary data vectors in the following order:
!
! ints: none
!
! doubles:
!   free ion conc <N_primary>
!   primary_activity_coeff <N_primary>
!   secondary_activity_coeff <N_aqueous_complexes>
!   ion exchange ref cation sorbed conc <N_ion_exchange_sites>
!   surface complexation free site conc <N_surface_sites>
!
! **************************************************************************** !
subroutine GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot, &
                                 num_ints, num_doubles)

  use, intrinsic :: iso_c_binding, only : c_int

  ! function parameters  
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot

  integer (c_int), intent(out) :: num_ints
  integer (c_int), intent(out) :: num_doubles

  num_ints = 1

  num_doubles = &
       3 * ncomp + &
       nspec + &
       2 * nexchange + &
       ncomp + &
       2 * nsurf + &
       npot + &
       2 * ncomp + &
       1

end subroutine GetAuxiliaryDataSizes

! **************************************************************************** !
subroutine PackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                     nexchange, nsurf, ndecay, npot, &
                                     aux_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData

  use crunchtype
  use concentration, only : sp10, sp, gam, &
                                           spex, spex10, &
                                           spsurf, LogTotalSurface, &
                                           sexold, ssurfold, &
                                           skdold, &
                                           jinit
  use mineral, only : LogPotential
  use medium, only : porin

  ! function parameters  
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  integer (c_int) :: num_ints, num_doubles
  integer(c_int), pointer :: idata(:)
  real (c_double), pointer :: data(:)
  integer :: i, dindex, idindex
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1

  call GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             num_ints, num_doubles)

! integers
  call c_f_pointer(aux_data%aux_ints%data, idata, (/num_ints/))
  idindex = 0

 ! jinit
  idindex = idindex + 1
  idata(idindex) = jinit(jx,jy,jz)

! doubles
  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0

  ! sp10 (primary species concs)
  do i = 1, ncomp
     dindex = dindex + 1
     data(dindex) = sp10(i,jx,jy,jz)
  end do

  ! sp (primary species log concs)
  do i = 1, ncomp
     dindex = dindex + 1
     data(dindex) = sp(i,jx,jy,jz)
  end do

  ! primary activity coeff
  do i = 1, ncomp
     dindex = dindex + 1
     data(dindex) = gam(i,jx,jy,jz)
  end do

  ! secondary aqueous complexe activity coeffs
  do i = 1, nspec
     dindex = dindex + 1
     data(dindex) = gam(i+ncomp,jx,jy,jz)
  end do

  ! ion exchange
  do i = 1, nexchange
     dindex = dindex + 1
     data(dindex) = spex(i,jx,jy,jz)
  end do  

  ! ion exchange
  do i = 1, nexchange
     dindex = dindex + 1
     data(dindex) = spex10(i,jx,jy,jz)
  end do  

  ! ion exchange
  do i = 1, ncomp
     dindex = dindex + 1
     data(dindex) = sexold(i,jx,jy,jz)
!     write(*,*)'pack sexold: ',sexold(i,jx,jy,jz)
  end do  

  ! surface primary
  do i = 1, nsurf
     dindex = dindex + 1
     data(dindex) = spsurf(i,jx,jy,jz) 
  end do

  ! total surface
  do i = 1, nsurf
     dindex = dindex + 1
     data(dindex) = LogTotalSurface(i,jx,jy,jz) 
     !! write(*,*)'pack i, LogTotalSurface(i,jx,jy,jz) ',i, LogTotalSurface(i,jx,jy,jz) 
  end do

  ! surface potential
  do i = 1, npot
     dindex = dindex + 1
     data(dindex) = LogPotential(i,jx,jy,jz) 
  end do

! surface 
  do i = 1, ncomp
     dindex = dindex + 1
     data(dindex) = ssurfold(i,jx,jy,jz) 
     !! write(*,*)'ssurfold(i,jx,jy,jz): ',ssurfold(i,jx,jy,jz)
  end do

! surface 
  do i = 1, ncomp
     dindex = dindex + 1
     data(dindex) = skdold(i,jx,jy,jz) 
  end do

! initial porosity 
  dindex = dindex + 1
  data(dindex) = porin(jx,jy,jz) 

  return

end subroutine PackAlquimiaAuxiliaryData

! **************************************************************************** !
subroutine UnpackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                       nexchange, nsurf, ndecay, npot, &
                                       aux_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData
  
  use crunchtype
  use concentration, only : sp10, sp, spold, gam, &
                            spex, spex10, spexold, &
                            spsurf, LogTotalSurface, &
                            sexold, ssurfold, & 
                            skdold, &
                            jinit
  use mineral, only : LogPotential
  use medium, only: porin

  ! function parameters  
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data

  ! local variables
  integer (c_int) :: num_ints, num_doubles
  integer (c_int), pointer :: idata(:)
  real (c_double), pointer :: data(:)
  integer :: i, dindex, idindex
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1

  call GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             num_ints, num_doubles)

! integers
  call c_f_pointer(aux_data%aux_ints%data, idata, (/num_ints/))
  idindex = 0

 ! jinit
  idindex = idindex + 1
  jinit(jx,jy,jz) = idata(idindex)

! doubles
  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0

  ! sp10 (primary species concs)
  do i = 1, ncomp
     dindex = dindex + 1
     sp10(i,jx,jy,jz) = data(dindex)
  end do

  ! sp (primary species log concs)
  do i = 1, ncomp
     dindex = dindex + 1
     sp(i,jx,jy,jz) = data(dindex)
     spold(i,jx,jy,jz) = data(dindex)
  end do

  ! primary activity coeff
  do i = 1, ncomp
     dindex = dindex + 1
     gam(i,jx,jy,jz) = data(dindex)
  end do

  ! aqueous complexes activity coeff
  do i = 1, nspec
     dindex = dindex + 1
     gam(i+ncomp,jx,jy,jz) = data(dindex)
  end do

  ! ion exchange
  do i = 1, nexchange
     dindex = dindex + 1
     spex(i,jx,jy,jz) = data(dindex)
     spexold(i,jx,jy,jz) = data(dindex)
  end do

  ! ion exchange
  do i = 1, nexchange
     dindex = dindex + 1
     spex10(i,jx,jy,jz) = data(dindex)
  end do

 ! ion exchange
  do i = 1, ncomp
     dindex = dindex + 1
     sexold(i,jx,jy,jz) = data(dindex)
  end do

! equilibrium surface complexation
  do i = 1, nsurf
     dindex = dindex + 1
     spsurf(i,jx,jy,jz) = data(dindex)
  end do

  ! equilibrium surface complexation
  do i = 1, nsurf
     dindex = dindex + 1
     LogTotalSurface(i,jx,jy,jz) = data(dindex)
  end do

    ! equilibrium surface complexation
  do i = 1, npot
     dindex = dindex + 1
     LogPotential(i,jx,jy,jz) = data(dindex)
  end do

  ! equilibrium surface complexation
  do i = 1, ncomp
     dindex = dindex + 1
     ssurfold(i,jx,jy,jz) = data(dindex)
  end do

  ! equilibrium surface complexation
  do i = 1, ncomp
     dindex = dindex + 1
     skdold(i,jx,jy,jz) = data(dindex)
  end do

! initial porosity 
  dindex = dindex + 1
  porin(jx,jy,jz) = data(dindex)

  return

end subroutine UnpackAlquimiaAuxiliaryData

! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PrintSizes(sizes)

  use AlquimiaContainers_module, only : AlquimiaSizes

  implicit none

  ! function parameters
  type (AlquimiaSizes), intent(in) :: sizes

  write (*, '(a)') "size : "
  write (*, '(a, i4)') "  num primary : ", sizes%num_primary
  write (*, '(a, i4)') "  num minerals : ", sizes%num_minerals
  write (*, '(a, i4)') "  num aqueous complexes : ", sizes%num_aqueous_complexes
  write (*, '(a, i4)') "  num surface sites : ", sizes%num_surface_sites
  write (*, '(a, i4)') "  num ion exchange sites : ", sizes%num_ion_exchange_sites
  write (*, '(a, i4)') "  num aux integers : ", sizes%num_aux_integers
  write (*, '(a, i4)') "  num aux doubles : ", sizes%num_aux_doubles
  return

end subroutine PrintSizes


! **************************************************************************** !
subroutine PrintState(state)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaState

  implicit none

  ! function parameters
  type (AlquimiaState), intent(in) :: state

  ! local variables
  integer :: i
  real (c_double), pointer :: conc(:)
  real (c_double), pointer :: immo(:)
  real (c_double), pointer :: vofx(:)
  real (c_double), pointer :: surf(:)
  real (c_double), pointer :: sites(:)
  real (c_double), pointer :: cec(:)
  
  write (*, '(a)') "state : "
  write (*, '(a, 1es13.6)') "  density water : ", state%water_density
  write (*, '(a, 1es13.6)') "  porosity : ", state%porosity
  write (*, '(a, 1es13.6)') "  temperature : ", state%temperature
  write (*, '(a, 1es13.6)') "  aqueous pressure : ", state%aqueous_pressure
  write (*, '(a, i4, a)') "  total mobile (", state%total_mobile%size, ") : "
  call c_f_pointer(state%total_mobile%data, conc, (/state%total_mobile%size/))
  do i=1, state%total_mobile%size
     write (*, '(1es13.6)') conc(i)
  end do
  write (*, '(a, i4, a)') "  total immobile (", state%total_immobile%size, ") : "
  call c_f_pointer(state%total_immobile%data, immo, (/state%total_immobile%size/))
  do i=1, state%total_immobile%size
     write (*, '(1es13.6)') immo(i)
  end do
  write (*, '(a, i4, a)') "  mineral volume fraction (", state%mineral_volume_fraction%size, ") : "
  call c_f_pointer(state%mineral_volume_fraction%data, vofx, (/state%mineral_volume_fraction%size/))
  do i=1, state%mineral_volume_fraction%size
     write (*, '(1es13.6)') vofx(i)
  end do
  write (*, '(a, i4, a)') "  mineral surface area (", state%mineral_specific_surface_area%size, ") : "
  call c_f_pointer(state%mineral_specific_surface_area%data, surf, (/state%mineral_specific_surface_area%size/))
  do i=1, state%mineral_specific_surface_area%size
     write (*, '(1es13.6)') surf(i)
  end do
  write (*, '(a, i4, a)') "  surface site density (", state%surface_site_density%size, ") : "
  call c_f_pointer(state%surface_site_density%data, sites, (/state%surface_site_density%size/))
  do i=1, state%surface_site_density%size
     write (*, '(1es13.6)') sites(i)
  end do
  write (*, '(a, i4, a)') "  cation exchange capacity (", state%cation_exchange_capacity%size, ") : "
  call c_f_pointer(state%cation_exchange_capacity%data, cec, (/state%cation_exchange_capacity%size/))
  do i=1, state%cation_exchange_capacity%size
     write (*, '(1es13.6)') cec(i)
  end do
end subroutine PrintState


! **************************************************************************** !
subroutine PrintEngineFunctionality(functionality)

  use AlquimiaContainers_module, only : AlquimiaEngineFunctionality

  implicit none

  ! function parameters
  type (AlquimiaEngineFunctionality), intent(in) :: functionality

  write (*, '(a)') "functionality : "
  write (*, '(a, L1)') "  thread safe : ", functionality%thread_safe
  write (*, '(a, L1)') "  temperature dependent : ", functionality%temperature_dependent
  write (*, '(a, L1)') "  pressure dependent : ", functionality%pressure_dependent
  write (*, '(a, L1)') "  porosity update : ", functionality%porosity_update
  write (*, '(a, i4)') "  index base : ", functionality%index_base
end subroutine PrintEngineFunctionality

! **************************************************************************** !
subroutine PrintProblemMetaData(meta_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer

  use c_f_interface_module, only : c_f_string_ptr

  use AlquimiaContainers_module, only : &
       AlquimiaProblemMetaData, kAlquimiaMaxStringLength

  implicit none

  ! function parameters
  type (AlquimiaProblemMetaData), intent(in) :: meta_data

  ! local variables
  type (c_ptr), pointer :: names(:)
  character(len=kAlquimiaMaxStringLength) :: name
  integer (c_int) :: i

  write (*, '(a)') "meta_data : "
  write (*, '(a)') "  primary names : "
  call c_f_pointer(meta_data%primary_names%data, names, (/meta_data%primary_names%size/))
  do i=1, meta_data%primary_names%size
     call c_f_string_ptr(names(i), name)
     write (*, '(a)') trim(name)
  end do
end subroutine PrintProblemMetaData


! **************************************************************************** !
subroutine PrintStatus(status)

  use AlquimiaContainers_module, only : AlquimiaEngineStatus

  implicit none

  ! function parameters
  type (AlquimiaEngineStatus), intent(in) :: status

  write (*, '(a)') "status : "
  write (*, '(a)') "  num rhs evaluation  : ", status%num_rhs_evaluations
  write (*, '(a)') "  num jacobian evaluations : ", status%num_jacobian_evaluations
  write (*, '(a)') "  num newton iterations : ", status%num_newton_iterations
  write (*, '(a)') "  converged  : ", status%converged

end subroutine PrintStatus

! Pflotran specific -- need to write those for CrunchFlow?
!!! **************************************************************************** !
!!subroutine PrintTranConstraint(tran_constraint)
!!
!!  use Constraint_module, only : tran_constraint_type
!!
!!  implicit none
!!
!!  ! function parameters
!!  type (tran_constraint_type), pointer :: tran_constraint
!!
!!  write (*, '(a)') "TranConstraint :"
!!  write (*, '(a, i4)') "    id : ", tran_constraint%id
!!  write (*, '(a, a)') "    name : ", tran_constraint%name
!!  write (*, '(a, L1)') "    requires equilibration : ", tran_constraint%requires_equilibration
!!  call PrintAqueousSpeciesConstraint(tran_constraint%aqueous_species)
!!  call PrintMineralConstraint(tran_constraint%minerals)
!!!    type(mineral_constraint_type), pointer :: minerals
!!!    type(srfcplx_constraint_type), pointer :: surface_complexes
!!!    type(colloid_constraint_type), pointer :: colloids
!!!    type(immobile_constraint_type), pointer :: immobile_species
!!
!!
!!end subroutine PrintTranConstraint
!!
!!subroutine PrintAqueousSpeciesConstraint(aqueous_species)
!!
!!  use Reaction_aux_module, only : aq_species_constraint_type
!!
!!  implicit none
!!
!!  type (aq_species_constraint_type), pointer :: aqueous_species 
!!
!!  write (*, '(a)') "    Aqueous species :"
!!  write (*, '(a)') aqueous_species%names
!!  write (*, '(a)') "    Constraint Conc :"
!!  write (*, '(f18.8)') aqueous_species%constraint_conc(:)
!!  write (*, '(a)') "    Constraint basis molarity :"
!!  write (*, '(f18.8)') aqueous_species%basis_molarity(:)
!!  write (*, '(a)') "    Constraint type :"
!!  write (*, '(i4)') aqueous_species%constraint_type(:)
!!  write (*, '(a)') "    Constraint spec id :"
!!  write (*, '(i4)') aqueous_species%constraint_spec_id(:)
!!  write (*, '(a)') "    Constraint aux string :"
!!  write (*, '(a)') aqueous_species%constraint_aux_string(:)
!!
!!end subroutine PrintAqueousSpeciesConstraint
!!
!!subroutine PrintMineralConstraint(minerals)
!!
!!  use Mineral_aux_module, only : mineral_constraint_type
!!
!!  implicit none
!!
!!  type (mineral_constraint_type), pointer :: minerals
!!
!!  write (*, '(a)') "    Mineral species :"
!!  write (*, '(a)') minerals%names
!!  write (*, '(a)') "    volume fraction :"
!!  write (*, '(f18.8)') minerals%constraint_vol_frac(:)
!!  write (*, '(a)') "    Constraint area :"
!!  write (*, '(f18.8)') minerals%constraint_area(:)
!!  write (*, '(a)') "    Constraint aux string :"
!!  write (*, '(a)') minerals%constraint_aux_string(:)
!!
!!end subroutine PrintMineralConstraint

subroutine PrintAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                              nexchange, nsurf, ndecay, npot, &
                              aux_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData

  use crunchtype
  use concentration, only : sp10, sp, gam, &
                            spex, spex10, &
                            spsurf, LogTotalSurface, &
                            sexold, ssurfold, &
                            skdold, &
                            jinit
  use mineral, only : LogPotential
  use medium, only : porin

  ! function parameters  
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  integer (c_int) :: num_ints, num_doubles
  integer(c_int), pointer :: idata(:)
  real (c_double), pointer :: data(:)
  integer :: i, dindex, idindex
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1

  ! local variables
  real (c_double), pointer :: conc(:)
  real (c_double), pointer :: immo(:)
  real (c_double), pointer :: vofx(:)
  real (c_double), pointer :: surf(:)
  real (c_double), pointer :: sites(:)
  real (c_double), pointer :: cec(:)
  

  call GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             num_ints, num_doubles)

  write (*, '(a)') "auxiliary data for crunchflow: "
  write (*, '(a)') "-- integers: "

! integers
  call c_f_pointer(aux_data%aux_ints%data, idata, (/num_ints/))
  idindex = 0

! jinit
  idindex = idindex + 1
  write (*, '(a, i4)') "  jinit : ",   idata(idindex)

  write (*, '(a)') "-- doubles: "
  
! doubles
  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0
  
  ! sp10 (primary species concs)
  write (*, '(a, i4, a)') "  primary species concs (", ncomp, ") : "
  do i = 1, ncomp
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do
  
  ! sp (primary species log concs)
  write (*, '(a, i4, a)') "  primary species log concs (", ncomp, ") : "
  do i = 1, ncomp
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! primary activity coeff
  write (*, '(a, i4, a)') "  primary activity coeff (", ncomp, ") : "
  do i = 1, ncomp
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! secondary aqueous complexes activity coeffs
  write (*, '(a, i4, a)') "  secondary aqueous complexes activity coeffs (", nspec, ") : "
  do i = 1, nspec
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! ion exchange
  write (*, '(a, i4, a)') "  ion exchange spex (", nexchange, ") : "
  do i = 1, nexchange
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do  

  ! ion exchange
  write (*, '(a, i4, a)') "  ion exchange spex10 (", nexchange, ") : "
  do i = 1, nexchange
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do  

  ! ion exchange
  write (*, '(a, i4, a)') "  ion exchange sexold (", ncomp, ") : "
  do i = 1, ncomp
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! surface primary
  write (*, '(a, i4, a)') "  surface primary spsurf (", nsurf, ") : "
  do i = 1, nsurf
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! total surface
  write (*, '(a, i4, a)') "  total surface LogTotalSurface (", nsurf, ") : "
  do i = 1, nsurf
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! surface potential
  write (*, '(a, i4, a)') "  surface potential LogPotential (", npot, ") : "
  do i = 1, npot
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! surface
  write (*, '(a, i4, a)') "  surface total old ssurfold (", ncomp, ") : "
  do i = 1, ncomp
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

  ! surface
  write (*, '(a, i4, a)') "  kd surface conc old skold (", ncomp, ") : "
  do i = 1, ncomp
     dindex = dindex + 1
     write (*, '(1es13.6)') data(dindex)
  end do

! initial porosity 
  dindex = dindex + 1
  write (*, '(a, 1es13.6)') "  initial porosity : ", data(dindex)

  return

end subroutine PrintAuxiliaryData

end module CrunchAlquimiaInterface_module

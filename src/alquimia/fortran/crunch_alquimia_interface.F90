
!
! Alquimia Copyright (c) 2013, The Regents of the University of California, 
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
! Author: Sergi Molins <smolins@lbl.gov>
! using plfotran_alquimia_interface by Ben Andre as template
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
!       ConvertAlquimiaConditionToPflotran, &
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
       PrintStatus

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
     integer(i4b)             :: nx
     integer(i4b)             :: ny
     integer(i4b)             :: nz
     integer(i4b)             :: igamma
     integer(i4b)             :: ikph
     integer(i4b)             :: jpor
     real(dp)                 :: corrmax
     real(dp)                 :: deltmin
     real(dp)                 :: time
  end type CrunchEngineState
  
contains


! **************************************************************************** !
!
! Public interface routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine Setup(input_filename, cf_engine_state, sizes, functionality, status)
!  NOTE: Function signature is dictated by the alquimia API.
!
!  NOTE: Assumes that MPI_Init() and / or PetscInitialize() have already
!    been called by the driver

  use, intrinsic :: iso_c_binding, only : c_char, c_ptr

  use c_f_interface_module, only : f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow global variables
  use crunchtype
  use params, only: mls
  use medium, only: isaturate

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
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
  CALL start98(ncomp,nspec,nkin,nrct,ngas,npot,nx,ny,nz,data1,ipath,igamma,  &
               ikmast,ikph,iko2,ltitle,tstep,delt,deltmin,ttol,jpor,ikin,nstop,       &
               corrmax,nseries,nexchange,nexch_sec,nsurf,nsurf_sec,ndecay,str_mon,    &
               str_day,str_hr,str_min,str_sec,NumInputFiles,InputFileCounter,alquimia)

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
                        sizes)

  ! Engine functionality
  ! note: no need to pass anything for now
  call SetEngineFunctionality(functionality)
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
subroutine ProcessCondition(cf_engine_state, condition, material_properties, &
     state, aux_data, status)
!  NOTE: Function signature is dictated by the alquimia API.

  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  use c_f_interface_module, only : c_f_string_ptr, f_c_string_ptr

  use AlquimiaContainers_module

  ! crunchflow
  use crunchtype
  
  use runtime, only: nchem
  use concentration, only: condlabel, stmp

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaGeochemicalCondition), intent(in) :: condition
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
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

  ! NOTE(bja): the data stored in alquimia's aux_data is uninitialized
  ! at this point, so don't want to copy it! (copy_auxdata = false)
  call CopyAlquimiaToAuxVars(copy_auxdata, state, aux_data, material_properties, &
                             ncomp, nspec, nkin, nrct, ngas, nexchange, nsurf, ndecay, npot)
  !
  ! process the condition
  !
  call  c_f_string_ptr(condition%name, name)
  !! note: need to replace with crunchflow print message call for rank 0!!
  !!engine_state%option%io_buffer = "processing : " // trim(name)
  !!call printMsg(engine_state%option) 

  if (condition%aqueous_constraints%size > 0) then
     ! the driver is supplying the constraint data, so we need to
     ! construct a geochemical condition in crunchflow's internal format.

     !! (smr) fixme: not supported yet
     status%error = kAlquimiaErrorUnsupportedFunctionality
     call f_c_string_ptr("ERROR: crunchflow interface does not support externally supplied geochemical conditions!", &
          status%message, kAlquimiaMaxStringLength)
     return
  
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

  if (allocated(stmp)) DEALLOCATE(stmp) !! to revisit ?
  allocate(stmp(ncomp))

! call to crunchflow subroutines are packaged in ProcessCrunchConstraint for clarity
!     
  call ProcessCrunchConstraint(engine_state, nco)

  if (allocated(stmp)) DEALLOCATE(stmp) !! to revisit ?
!
! repack the processed constraint data into the alquimia struct for the driver
!
  call CopyAuxVarsToAlquimia(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             state, aux_data )
  status%error = kAlquimiaNoError

end subroutine ProcessCondition


! **************************************************************************** !
subroutine ReactionStepOperatorSplit(cf_engine_state, &
     delta_t, material_properties, state, aux_data, status)
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
                           xgram, xgramOld
  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  real (c_double), intent(in) :: delta_t
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
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

  !call PrintState(engine_state%reaction, state)

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

  call CopyAlquimiaToAuxVars(copy_auxdata, state, aux_data, material_properties, &
                             ncomp, nspec, nkin, nrct, ngas, nexchange, nsurf, ndecay, npot)

  IF (nexchange > 0) THEN
    CALL UpdateExchanger(nx,ny,nz,nexchange)
  END IF

! set delt variable in CrunchFlow units
! IN= seconds, CF=years
  delt = delta_t / secyr

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

    ! --> error reduce time step
    write(*,*)'dead '
    stop

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
      ! update old concentrations to be used in next time step as initial guesses
      ! straight from CrunchFlow.f90
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            CALL oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)  
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
      
      ! send state variables back to Alquimia state with solution of this solve
      call CopyAuxVarsToAlquimia(ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot, &
                                 state, aux_data )
    
      ! write stats of solve into status
      status%num_newton_iterations = iterat
      status%error = kAlquimiaNoError

    end if
  
  end if   

end subroutine ReactionStepOperatorSplit


! **************************************************************************** !
subroutine GetAuxiliaryOutput( &
     cf_engine_state, &
     material_properties, &
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
  type (AlquimiaMaterialProperties), intent(in) :: material_properties
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
                           namexc
  use mineral, only: umin

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: cf_engine_state
  type (AlquimiaProblemMetaData), intent(out) :: meta_data
  type (AlquimiaEngineStatus), intent(out) :: status

  ! local variables
  type (c_ptr), pointer :: name_list(:)
  character (c_char), pointer :: name
  integer :: i, list_size
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

  if (meta_data%primary_names%size /= engine_state%ncomp) then
     write (*, '(a, i3, a, i3, a)') "meta_data%primary_names%size (", &
          meta_data%primary_names%size, ") != crunchflow%ncomp(", &
          engine_state%ncomp, ")"
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
  ! CF does not track individual Kd species, it is all primary or none
  ! (smr) TODO
  list_size = meta_data%isotherm_species_names%size
  call c_f_pointer(meta_data%isotherm_species_names%data, name_list, &
       (/list_size/))
  do i = 1, list_size
     call c_f_pointer(name_list(i), name)
     call f_c_string_chars( &
          trim(ulab(i)), &
          name, kAlquimiaMaxStringLength)
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
  
  ! set the pflotran input file name
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
subroutine SetEngineFunctionality(functionality)

  use AlquimiaContainers_module, only : AlquimiaEngineFunctionality
    
  implicit none

  ! function parameters
  type (AlquimiaEngineFunctionality), intent(out) :: functionality

  functionality%thread_safe = .false. ! (smr) crunch global variables do not allow for thread safe engine
  functionality%temperature_dependent = .false.
  functionality%pressure_dependent = .false.
  functionality%porosity_update = .false.
  functionality%operator_splitting = .true.
  functionality%global_implicit = .false.
  functionality%index_base = 1

end subroutine SetEngineFunctionality

! **************************************************************************** !
subroutine SetAlquimiaSizes(ncomp, nspec, nkin, nrct, ngas, &
                            nexchange, nsurf, ndecay, npot, &
                            sizes)

! note(smr): zeros denote unavailable processes for now, aug 8 2014

  use AlquimiaContainers_module, only : AlquimiaSizes

  use crunchtype

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
  type (AlquimiaSizes), intent(out) :: sizes

  sizes%num_primary = ncomp
  if (nsurf > 0 .or. nexchange > 0) then
     sizes%num_sorbed = ncomp ! as per alquimia convention
  else
     sizes%num_sorbed = 0
  end if
  sizes%num_kinetic_minerals = nkin ! or nrct -- need to check what each one is
  sizes%num_aqueous_complexes = nspec
  sizes%num_surface_sites = nsurf
  sizes%num_ion_exchange_sites = nexchange
  ! number of retardation species are not tracked explicitly in CF outside start98
  ! the array distrib (declared in concentration.f90) is allocated for ncomp
  ! with distrib(i) = kd if i retarded, and distrib(i) = 0 if not 
  ! set to zero for now and figure out things later 
  !! TODO (smr)
  sizes%num_isotherm_species = 0
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
                           xgram, xgramOld
  use mineral, only: LogPotential_tmp, iedl, &
                     volfx, volin, &
                     VolumeLastTimeStep, &
                     area, areain, &
                     mintype, LogPotential
  use medium, only: AqueousToBulkCond, SaturationCond, constantpor, &
                    porin, por, porOld, porcond, isaturate                   
  use temperature, only: tempcond, rocond, ro
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
 
  spgastmp10 = spgastmp10*denmol
  spgastmp = DLOG(spgastmp10)

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
        if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz)
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
  CALL oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)  
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

!! TO DO
!!!!! **************************************************************************** !
!!!!function ConvertAlquimiaConditionToCrunch(&
!!!!     option, reaction, alquimia_condition)
!!!!  use, intrinsic :: iso_c_binding
!!!!
!!!!  use c_f_interface_module, only : c_f_string_ptr
!!!!
!!!!  use AlquimiaContainers_module
!!!!
!!!!  ! pflotran
!!!!  use Option_module, only : option_type, printErrMsg, printMsg
!!!!  use Reaction_aux_module, only : reaction_type, aq_species_constraint_type, &
!!!!       AqueousSpeciesConstraintCreate
!!!!  use Mineral_aux_module, only : mineral_constraint_type, MineralConstraintCreate
!!!!  use String_module, only : StringCompareIgnoreCase
!!!!  use Constraint_module, only : tran_constraint_type, TranConstraintCreate, &
!!!!       CONSTRAINT_FREE, CONSTRAINT_TOTAL, CONSTRAINT_TOTAL_SORB, &
!!!!       CONSTRAINT_PH, CONSTRAINT_MINERAL, &
!!!!       CONSTRAINT_GAS, CONSTRAINT_CHARGE_BAL
!!!!
!!!!
!!!!  implicit none
!!!!
!!!!  ! function parameters
!!!!  type(option_type), pointer, intent(in) :: option
!!!!  type(reaction_type), pointer, intent(in) :: reaction
!!!!  type (AlquimiaGeochemicalCondition), intent(in) :: alquimia_condition
!!!!
!!!!  ! Return value
!!!!  type (tran_constraint_type), pointer :: ConvertAlquimiaConditionToPflotran
!!!!
!!!!  ! local variables
!!!!  integer :: i
!!!!  character (kAlquimiaMaxStringLength) :: name, constraint_type
!!!!  character (kAlquimiaMaxStringLength) :: associated_species
!!!!  type (tran_constraint_type), pointer :: tran_constraint
!!!!  type(aq_species_constraint_type), pointer :: pft_aq_species_constraint
!!!!  type(mineral_constraint_type), pointer :: pft_mineral_constraint
!!!!  type (AlquimiaAqueousConstraint), pointer :: alq_aqueous_constraints(:)
!!!!  type (AlquimiaMineralConstraint), pointer :: alq_mineral_constraints(:)
!!!!
!!!!
!!!!  call c_f_string_ptr(alquimia_condition%name, name)
!!!!  option%io_buffer = "building : " // trim(name)
!!!!  call printMsg(option)
!!!!  tran_constraint => TranConstraintCreate(option)
!!!!  tran_constraint%name = trim(name)
!!!!  ! NOTE(bja): requires_equilibration not used in pflotran?
!!!!  tran_constraint%requires_equilibration = PETSC_FALSE
!!!!
!!!!  !
!!!!  ! aqueous species
!!!!  !
!!!!  if (alquimia_condition%aqueous_constraints%size /= reaction%naqcomp) then
!!!!     option%io_buffer = 'Number of aqueous constraints ' // &
!!!!          'does not equal the number of primary chemical ' // &
!!!!          'components in constraint: ' // &
!!!!          trim(tran_constraint%name)
!!!!     call printErrMsg(option)
!!!!  end if
!!!!
!!!!  ! NOTE(bja) : this is the container for ALL aqueous constraints
!!!!  pft_aq_species_constraint => &
!!!!       AqueousSpeciesConstraintCreate(reaction, option)
!!!!
!!!!  call c_f_pointer(alquimia_condition%aqueous_constraints%data, &
!!!!       alq_aqueous_constraints, (/alquimia_condition%aqueous_constraints%size/))
!!!!
!!!!  do i = 1, alquimia_condition%aqueous_constraints%size
!!!!     call c_f_string_ptr(alq_aqueous_constraints(i)%primary_species_name, name)
!!!!     pft_aq_species_constraint%names(i) = trim(name)
!!!!
!!!!     pft_aq_species_constraint%constraint_conc(i) = alq_aqueous_constraints(i)%value
!!!!
!!!!     call c_f_string_ptr(alq_aqueous_constraints(i)%constraint_type, constraint_type)
!!!!
!!!!     call c_f_string_ptr(alq_aqueous_constraints(i)%associated_species, &
!!!!          associated_species)
!!!!
!!!!     if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringFree)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_FREE
!!!!
!!!!     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringTotalAqueous)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_TOTAL
!!!!
!!!!     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringTotalSorbed)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_TOTAL_SORB
!!!!     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringPH)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_PH
!!!!
!!!!     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringMineral)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_MINERAL
!!!!        pft_aq_species_constraint%constraint_aux_string(i) = &
!!!!             trim(associated_species)
!!!!
!!!!     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringGas)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_GAS
!!!!        pft_aq_species_constraint%constraint_aux_string(i) = &
!!!!             trim(associated_species)
!!!!
!!!!     else if (StringCompareIgnoreCase(constraint_type, kAlquimiaStringCharge)) then
!!!!        pft_aq_species_constraint%constraint_type(i) = CONSTRAINT_CHARGE_BAL
!!!!
!!!!     else
!!!!        option%io_buffer = 'Constraint type: ' // trim(constraint_type) // &
!!!!             ' not recognized in constraint,concentration'
!!!!        call printErrMsg(option)
!!!!     end if
!!!!  end do
!!!!  tran_constraint%aqueous_species => pft_aq_species_constraint
!!!!
!!!!  !
!!!!  ! minerals
!!!!  !
!!!!  ! FIXME(bja): are these checks the correct thing to do in all cases...?
!!!!  if (alquimia_condition%mineral_constraints%size > 0 .and. &
!!!!       alquimia_condition%mineral_constraints%size /= &
!!!!       reaction%mineral%nkinmnrl) then
!!!!     option%io_buffer = &
!!!!          'Number of mineral constraints is not equal to ' // &
!!!!          'number of kinetic minerals in condition: ' // &
!!!!          trim(tran_constraint%name)
!!!!     call printErrMsg(option)
!!!!  end if
!!!!
!!!!  pft_mineral_constraint => MineralConstraintCreate(reaction%mineral, option)
!!!!
!!!!  call c_f_pointer(alquimia_condition%mineral_constraints%data, &
!!!!       alq_mineral_constraints, (/alquimia_condition%mineral_constraints%size/))
!!!!  do i = 1, alquimia_condition%mineral_constraints%size
!!!!     call c_f_string_ptr(alq_mineral_constraints(i)%mineral_name, name)
!!!!     pft_mineral_constraint%names(i) = trim(name)
!!!!     pft_mineral_constraint%constraint_vol_frac(i) = &
!!!!          alq_mineral_constraints(i)%volume_fraction
!!!!     pft_mineral_constraint%constraint_area(i) = &
!!!!          alq_mineral_constraints(i)%specific_surface_area
!!!!  end do
!!!!  tran_constraint%minerals => pft_mineral_constraint
!!!!
!!!!  ConvertAlquimiaConditionToPflotran => tran_constraint
!!!!
!!!!end function ConvertAlquimiaConditionToCrunch


! **************************************************************************** !
subroutine CopyAlquimiaToAuxVars(copy_auxdata, state, aux_data, material_prop, &
                                 ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module

  ! crunchflow
  use crunchtype
  use temperature, only : t, ro
  use transport, only : satliq
  use medium, only : por
  use concentration, only : sn, &
                            ssurfold, &
                            sexold, &
                            exchangesites, &
                            ssurfn

  use mineral, only : volfx, area

  implicit none

  ! function parameters
  logical, intent(in) :: copy_auxdata
  type (AlquimiaState), intent(in) :: state
  type (AlquimiaAuxiliaryData), intent(in) :: aux_data
  type (AlquimiaMaterialProperties), intent(in) :: material_prop
  integer(i4b), intent(in) :: ncomp
  integer(i4b), intent(in) :: nspec
  integer(i4b), intent(in) :: nkin
  integer(i4b), intent(in) :: nrct
  integer(i4b), intent(in) :: ngas
  integer(i4b), intent(in) :: nexchange
  integer(i4b), intent(in) :: nsurf
  integer(i4b), intent(in) :: ndecay
  integer(i4b), intent(in) :: npot

  ! local variables
  real (c_double), pointer :: data(:)
  integer :: i
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1
  
  !write (*, '(a)') "Crunch_Alquimia_CopyAlquimiaToAuxVars() :"

  !
  ! copy state to crunchflow
  !
  ro(jx,jy,jz)            = state%water_density
  satliq(jx,jy,jz)        = material_prop%saturation
  t(jx,jy,jz)             = state%temperature
  !! not needed by crunch = state%aqueous_pressure
  por(jx,jy,jz)           = state%porosity  
  !! not needed by crunch = material_prop%volume
  !!write(*,*)'HERE dens, temp, por: ', state%water_density, state%temperature, state%porosity
  !
  ! primary aqueous
  !
  call c_f_pointer(state%total_mobile%data, data, (/ncomp/))
  do i = 1, ncomp
     sn(i,jx,jy,jz) = data(i) ! fixed point for nonlinear problem is in 'sn'
     !!write(*,*)'receive sn(i,jx,jy,jz): ',sn(i,jx,jy,jz),i
  end do

  ! sorbed primary
  if (nsurf > 0 .or. nexchange > 0) then
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
  ! isotherms (smr) leave out for now, see distrib
  !
  !call c_f_pointer(material_prop%isotherm_kd%data, data, &
  !     (/material_prop%isotherm_kd%size/))
  !do i = 1, material_prop%isotherm_kd%size
  !   reaction%eqkddistcoef(i) = data(i)
  !end do

  !call c_f_pointer(material_prop%langmuir_b%data, data, &
  !     (/material_prop%langmuir_b%size/))
  !do i = 1, material_prop%langmuir_b%size
  !   reaction%eqkdlangmuirb(i) = data(i)
  !end do

  !call c_f_pointer(material_prop%freundlich_n%data, data, &
  !     (/material_prop%freundlich_n%size/))
  !do i = 1, material_prop%freundlich_n%size
  !   reaction%eqkdfreundlichn(i) = data(i)
  !end do

  if (copy_auxdata) then
     call UnpackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                      nexchange, nsurf, ndecay, npot, &
                                      aux_data)
  end if

end subroutine CopyAlquimiaToAuxVars

! **************************************************************************** !
subroutine CopyAuxVarsToAlquimia(ncomp, nspec, nkin, nrct, ngas, &
                                 nexchange, nsurf, ndecay, npot, &
                                 state, aux_data)

  use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer

  use AlquimiaContainers_module

  ! crunch
  use crunchtype
  use temperature, only : t, ro
  !!use transport, only : satliq
  use medium, only : por !!, porin
  use concentration, only : sn, &
                             ssurfold, &
                             sexold, &
                             exchangesites, &
                             ssurfn
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
  type (AlquimiaState), intent(inout) :: state
  type (AlquimiaAuxiliaryData), intent(inout) :: aux_data

  ! local variables
  real (c_double), pointer :: data(:)
  integer :: i
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1
  
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
     data(i) = sn(i,jx,jy,jz)
     !!write(*,*)'sn(i,jx,jy,jz): ',sn(i,jx,jy,jz)
  end do  
  !
  ! sorbed
  !
  if (nsurf > 0 .or. nexchange > 0) then
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
  end do

  ! NOTE(bja): isotherms are material properties, and can't be changed
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

  num_ints = 0

  num_doubles = &
       3 * ncomp + &
       nspec + &
       nexchange + &
       2 * nsurf + &
       npot

end subroutine GetAuxiliaryDataSizes

! **************************************************************************** !
subroutine PackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                     nexchange, nsurf, ndecay, npot, &
                                     aux_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData

  use crunchtype
  use concentration, only : sp10, sp, gam, &
                            spex, spsurf, &
                            LogTotalSurface
  use mineral, only : LogPotential

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
  real (c_double), pointer :: data(:)
  integer :: i, dindex
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1

  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0
  
  call GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             num_ints, num_doubles)

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

  ! surface primary
  do i = 1, nsurf
     dindex = dindex + 1
     data(dindex) = spsurf(i,jx,jy,jz) 
  end do

  ! total surface
  do i = 1, nsurf
     dindex = dindex + 1
     data(dindex) = LogTotalSurface(i,jx,jy,jz) 
     !write(*,*)'pack i, LogTotalSurface(i,jx,jy,jz) ',i, LogTotalSurface(i,jx,jy,jz) 
  end do

  ! surface potential
  do i = 1, npot
     dindex = dindex + 1
     data(dindex) = LogPotential(i,jx,jy,jz) 
  end do

  return

end subroutine PackAlquimiaAuxiliaryData

! **************************************************************************** !
subroutine UnpackAlquimiaAuxiliaryData(ncomp, nspec, nkin, nrct, ngas, &
                                       nexchange, nsurf, ndecay, npot, &
                                       aux_data)

  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_f_pointer

  use AlquimiaContainers_module, only : AlquimiaAuxiliaryData
  
  use crunchtype
  use concentration, only : sp10, sp, gam, &
                            spex, spsurf, &
                            LogTotalSurface
  use mineral, only : LogPotential

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
  real (c_double), pointer :: data(:)
  integer :: i, dindex
  integer(i4b), parameter :: jx=1
  integer(i4b), parameter :: jy=1
  integer(i4b), parameter :: jz=1

  call GetAuxiliaryDataSizes(ncomp, nspec, nkin, nrct, ngas, &
                             nexchange, nsurf, ndecay, npot, &
                             num_ints, num_doubles)

  call c_f_pointer(aux_data%aux_doubles%data, data, (/num_doubles/))
  dindex = 0

  ! sp10 (primary species concs)
  do i = 1, ncomp
     dindex = dindex + 1
     sp10(i,jx,jy,jz) = data(dindex)
  end do

  ! sp10 (primary species log concs)
  do i = 1, ncomp
     dindex = dindex + 1
     sp(i,jx,jy,jz) = data(dindex)
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
  write (*, '(a, i4)') "  num kinetics minerals : ", sizes%num_kinetic_minerals
  write (*, '(a, i4)') "  num aqueous complexes : ", sizes%num_aqueous_complexes
  write (*, '(a, i4)') "  num surface sites : ", sizes%num_surface_sites
  write (*, '(a, i4)') "  num ion exchange sites : ", sizes%num_ion_exchange_sites
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

  write (*, '(a)') "state : "
  write (*, '(a, 1es13.6)') "  density water : ", state%water_density
  write (*, '(a, 1es13.6)') "  porosity : ", state%porosity
  write (*, '(a, 1es13.6)') "  temperature : ", state%temperature
  write (*, '(a, 1es13.6)') "  aqueous pressure : ", state%aqueous_pressure
  write (*, '(a, i4, a)') "  total primary (", state%total_mobile%size, ") : "
  call c_f_pointer(state%total_mobile%data, conc, (/state%total_mobile%size/))
  do i=1, state%total_mobile%size
     write (*, '(1es13.6)') conc(i)
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


end module CrunchAlquimiaInterface_module

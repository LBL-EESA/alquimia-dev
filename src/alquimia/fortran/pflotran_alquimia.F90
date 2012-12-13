! **************************************************************************** !
!
! Alquimia inteface
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PFloTranAlquimia_Setup(input_filename, metadata, sizes)

#include "pflotran_alquimia.h90"

  character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (alquimia_metadata_f), intent(inout) :: metadata
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

  sizes%num_primary = 1
  sizes%num_kinetic_minerals = 2
  sizes%num_aqueous_complexes = 3
  sizes%num_surface_sites = 4
  sizes%num_ion_exchange_sites = 5
  call PFloTranAlquimia_PrintSizes(sizes)

  metadata%thread_safe = .true.
  metadata%temperature_dependent = .false.
  metadata%pressure_dependent = .false.
  metadata%porosity_update = .true.
  call PFloTranAlquimia_PrintMetaData(metadata)

end subroutine PFloTranAlquimia_Setup


! **************************************************************************** !
subroutine PFloTranAlquimia_ProcessCondition()

#include "pflotran_alquimia.h90"

  print *, "Fortran process constraint."
end subroutine PFloTranAlquimia_ProcessCondition


! **************************************************************************** !
subroutine PFloTranAlquimia_ReactionStepOperatorSplit()

#include "pflotran_alquimia.h90"

  print *, "Fortran reaction step operator split."
end subroutine PFloTranAlquimia_ReactionStepOperatorSplit


! **************************************************************************** !
subroutine PFloTranAlquimia_GetAuxiliaryOutput()

#include "pflotran_alquimia.h90"

  print *, "Fortran get auxiliary output."
end subroutine PFloTranAlquimia_GetAuxiliaryOutput


! **************************************************************************** !
subroutine PFloTranAlquimia_GetEngineFunctionality(metadata)

#include "pflotran_alquimia.h90"

  type (alquimia_metadata_f), intent(out) :: metadata
  print *, "Fortran get engine functionality."
  metadata%thread_safe = .true.
  metadata%temperature_dependent = .true.
  metadata%pressure_dependent = .false.
  metadata%porosity_update = .false.

end subroutine PFloTranAlquimia_GetEngineFunctionality



! **************************************************************************** !
!
! Output helper routines
!
! **************************************************************************** !

! **************************************************************************** !
subroutine PFloTranAlquimia_PrintSizes(sizes)

#include "pflotran_alquimia.h90"

  type (alquimia_sizes_f), intent(inout) :: sizes
  print *, "size : "
  print *, "  num primary : ", sizes%num_primary
  print *, "  num kinetics minerals : ", sizes%num_kinetic_minerals
  print *, "  num aqueous complexes : ", sizes%num_aqueous_complexes
  print *, "  num surface sites : ", sizes%num_surface_sites
  print *, "  num ion exchange sites : ", sizes%num_ion_exchange_sites
end subroutine PFloTranAlquimia_PrintSizes


! **************************************************************************** !
subroutine PFloTranAlquimia_PrintMetadata(metadata)

#include "pflotran_alquimia.h90"

  type (alquimia_metadata_f), intent(inout) :: metadata
  print *, "metadata : "
  print *, "  thread safe : ", metadata%thread_safe
  print *, "  temperature dependent : ", metadata%temperature_dependent
  print *, "  pressure dependent : ", metadata%pressure_dependent
  print *, "  porosity update : ", metadata%porosity_update
end subroutine PFloTranAlquimia_PrintMetadata


! **************************************************************************** !
subroutine PFloTranAlquimia_PrintStatus(status)

#include "pflotran_alquimia.h90"

  type (alquimia_engine_status_f), intent(inout) :: status
  print *, "status : "
  print *, "  num rhs evaluation  : ", status%num_rhs_evaluations
  print *, "  num jacobian evaluations : ", status%num_jacobian_evaluations
  print *, "  num newton iterations : ", status%num_newton_iterations
  print *, "  converged  : ", status%converged

end subroutine PFloTranAlquimia_PrintStatus

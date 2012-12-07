#include "cc-f90-strings.h"

subroutine reverse_string(in_str, in_length, out_str)
  character(len=STRING_MAX_LENGTH), intent(in) :: in_str
  integer(4), intent(in) :: in_length
  character(len=STRING_MAX_LENGTH), intent(out) :: out_str
  integer(4) :: i

  do i = 1, in_length
     out_str(i:i) = in_str(in_length+1-i:in_length+1-i)
  end do
  ! NOTE: we need to do something with the null terminator...?
  !write(*, '(a)') out_str
end subroutine reverse_string

! NOTE: this fails for some reasaon, the length of the string is
! always MAX_LENGTH
!
!  integer :: length
!  character(len=STRING_MAX_LENGTH) :: temp
!  temp = trim(in_str)
!  length = len_trim(in_str)
!  write(*, *) length, temp
!  do i = 1, length
!    out_str(i:i) = in_str(length-i:length-i)
!  end do

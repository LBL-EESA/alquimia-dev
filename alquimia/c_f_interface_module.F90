
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

! Based on the c_interface_module API found at :
!   http://fortranwiki.org/fortran/show/c_interface_module
!
! NOTE(bja, 2013-07): Since that code isn't explicitly licensed, we
! can't use it.
!
module c_f_interface_module
  use, intrinsic :: iso_c_binding, only :

  implicit none

  private :: To_lower

  public :: &
       c_f_string_ptr, &
       c_f_string_chars, &
       f_c_string_ptr, &
       f_c_string_chars, &
       CaseInsensitiveStrcmp

contains

subroutine c_f_string_ptr(c_string_ptr, f_string)
  ! copy a null terminated c string into a fortran string.
  ! null pointer --> empty string
  ! empty fortran strings and padding at the end are expected to be spaces.
  use iso_c_binding, only : c_ptr, c_char, c_associated, c_f_pointer, c_null_char

  implicit none

  type(c_ptr), intent(in) :: c_string_ptr
  character(len=*), intent(out) :: f_string

  integer :: c
  character(kind=c_char, len=1), dimension(:), pointer :: c_chars

  f_string = ' '
  if (c_associated(c_string_ptr)) then
     call c_f_pointer(c_string_ptr, c_chars, [huge(0)])
     c = 1
     do while(c <= len(f_string) .and. c_chars(c) /= c_null_char)
        f_string(c:c) = c_chars(c)
        c = c + 1
     end do
  end if
end subroutine c_f_string_ptr

subroutine c_f_string_chars(c_string_chars, f_string)
  ! copy a c char array reference to fortran
  ! assumes c char array is properly null terminated to end the loop!
  use iso_c_binding, only : c_char, c_null_char

  implicit none
  character(kind=c_char, len=1), intent(in) :: c_string_chars(*)
  character(len=*), intent(out) :: f_string

  integer :: c

  f_string = ' '
  c = 1
  do while(c <= len(f_string) .and. c_string_chars(c) /= c_null_char)
     f_string(c:c) = c_string_chars(c)
     c = c + 1
  end do
end subroutine c_f_string_chars

subroutine f_c_string_ptr(f_string, c_string_ptr, c_string_len)
  ! copy a fortran string to a c string pointer
  use iso_c_binding, only : c_char, c_ptr, c_associated, c_f_pointer, c_null_char

  implicit none

  character(len=*), intent(in) :: f_string
  type(c_ptr), intent(in) :: c_string_ptr ! pointer can't be changed, target can
  integer, intent(in) :: c_string_len
 
  integer :: c, string_len
  character(kind=c_char, len=1), dimension(:), pointer :: c_chars

  string_len = min(len(f_string), c_string_len-1)
  if (c_associated(c_string_ptr)) then
     if (string_len >= 0) then
        call c_f_pointer(c_string_ptr, c_chars, [string_len+1])
        do c = 1, string_len
           c_chars(c) = f_string(c:c)
        end do
        c_chars(string_len + 1) = c_null_char
     end if
  end if

end subroutine f_c_string_ptr

subroutine f_c_string_chars(f_string, c_string_chars, c_string_len)
  ! copy a fortran string to a c string passed by char-array refrence
  ! c_string_len is the max string length, including null terminator
  use iso_c_binding, only : c_char, c_null_char

  implicit none

  character(len=*), intent(in) :: f_string
  character(kind=c_char, len=1), intent(out) :: c_string_chars(*)
  integer, intent(in) :: c_string_len 

  integer :: c, string_len

  string_len = min(len(f_string), c_string_len-1)

  do c = 1, string_len
     c_string_chars(c) = f_string(c:c)
  end do
  c_string_chars(string_len + 1) = c_null_char

end subroutine f_c_string_chars


function CaseInsensitiveStrcmp(str1, str2)
  ! this need to be here for alquimia to compile without pflotran
  implicit none

  character(len=*), intent(in) :: str1, str2
  character(len=512) :: low1, low2
  integer :: len1, len2, i
  logical :: CaseInsensitiveStrcmp

  len1 = len_trim(str1)
  len2 = len_trim(str2)

  if (len1 /= len2) then
    CaseInsensitiveStrcmp = .false.
    return
  endif

  low1 = str1
  low2 = str2 

  call To_lower(low1)
  call To_lower(low2)

  do i=1,len1
    if (low1(i:i) /= low2(i:i)) then
      CaseInsensitiveStrcmp = .false.
      return
    endif
  enddo

  CaseInsensitiveStrcmp = .true.
end function CaseInsensitiveStrcmp


subroutine To_lower(str)
  character(*), intent(in out) :: str
  integer :: i
  do i = 1, len(str)
    select case(str(i:i))
      case("A":"Z")
        str(i:i) = achar(iachar(str(i:i))+32)
    end select
  end do  
end subroutine To_Lower

end module c_f_interface_module

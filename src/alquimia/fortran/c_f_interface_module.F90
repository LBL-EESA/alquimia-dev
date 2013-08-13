! Based on the c_interface_module API found at :
!   http://fortranwiki.org/fortran/show/c_interface_module
!
! NOTE(bja, 2013-07): Since that code isn't explicitly licensed, we
! can't use it.
!
module c_f_interface_module
  use, intrinsic :: iso_c_binding, only :

  implicit none

  private

  public :: &
       c_f_string_ptr, &
       c_f_string_chars, &
       f_c_string_ptr, &
       f_c_string_chars

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

end module c_f_interface_module

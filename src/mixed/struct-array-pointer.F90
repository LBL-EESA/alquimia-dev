
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

#define MAX_STRING_LENGTH 256

! string processing functions from :
!   http://fortranwiki.org/fortran/show/c_interface_module
!
! From the fortranwiki.org copyright page:
! http://fortranwiki.org/fortran/show/Copyrights
!
! "All source code on the Fortran Wiki is in the Public domain unless
! otherwise noted."
!
! NOTE(bja): As of 2012-12-19, the c_interface_module page does not
! contain any additional copyright notice, so we can assume it is in
! the public domain.

! Copy a C string, passed by pointer, to a Fortran string.
! If the C pointer is NULL, the Fortran string is blanked.
! C_string must be NUL terminated, or at least as long as F_string.
! If C_string is longer, it is truncated. Otherwise, F_string is
! blank-padded at the end.
subroutine C_F_string_ptr(C_string, F_string)
  use, intrinsic :: iso_c_binding
  implicit none

  type(C_ptr), intent(in) :: C_string
  character(len=*), intent(out) :: F_string
  character(len=1,kind=C_char), dimension(:), pointer :: p_chars
  integer :: i
  if (.not. C_associated(C_string)) then
     F_string = ' '
  else
     call C_F_pointer(C_string,p_chars,[huge(0)])
     i=1
     do while(p_chars(i)/=c_null_char .and. i<=len(F_string))
        F_string(i:i) = p_chars(i)
        i=i+1
     end do
     if (i<len(F_string)) F_string(i:) = ' '
  end if
end subroutine C_F_string_ptr

! Copy a Fortran string to an allocated C string pointer.
! If the C pointer is NULL, no action is taken. (Maybe auto allocate via libc call?)
! If the length is not passed, the C string must be at least: len(F_string)+1
! If the length is passed and F_string is too long, it is truncated.
subroutine F_C_string_ptr(F_string, C_string, C_string_len)
  use, intrinsic :: iso_c_binding
  implicit none

  character(len=*), intent(in) :: F_string
  type(C_ptr), intent(in) :: C_string ! target = intent(out)
  integer, intent(in) :: C_string_len  ! Max string length,
  ! INCLUDING THE TERMINAL c_null_char
  character(len=1,kind=C_char), dimension(:), pointer :: p_chars
  integer :: i, strlen
  strlen = len(F_string)
  if (C_string_len <= 0) return
  strlen = min(strlen,C_string_len-1)
  if (.not. C_associated(C_string)) then
     return
  end if
  call C_F_pointer(C_string,p_chars,[strlen+1])
  forall (i=1:strlen)
     p_chars(i) = F_string(i:i)
  end forall
  p_chars(strlen+1) = c_null_char
end subroutine F_C_string_ptr


subroutine process_container(c)
  use, intrinsic :: iso_c_binding
  implicit none
  type, bind(c) :: container
     integer (c_int) :: num
     type (c_ptr) :: data
  end type container
  type, bind(c) :: data
     type (c_ptr) :: name
     type (c_ptr) :: outstr
     integer (c_int) :: xlength
     type (c_ptr) :: x
     integer (c_int) :: ylength
     type (c_ptr) :: y
  end type data
  type (container), intent(inout) :: c
  type (data), pointer :: local(:)
  character (len=MAX_STRING_LENGTH) :: name, outstr
  integer (c_int), pointer :: xlocal(:)
  real (c_double), pointer :: ylocal(:)
  integer :: i, j
  real :: ysum

  write (*, '(a)') "In fortran : "
  write (*, '(a i3)') "  c%num : ", c%num
  ! assign the array of data pointers to a fortran object
  call c_f_pointer(c%data, local, (/c%num/))

  ! loop through the local array of data objects
  do i = 1, c%num
     call C_F_string_ptr(local(i)%name, name)
     write (*, '(a i3 a a)') "   data(", i, ") : ", trim(name)
     ! associate xlocal with the c data
     call c_f_pointer(local(i)%x, xlocal, (/local(i)%xlength/))
     ! print the x data
     write (*, '(a i3)') "    f xlength : ", local(i)%xlength
     do j = 1, local(i)%xlength
        write (*, '     (i4)', advance='no') xlocal(j)
     enddo
     write(*, *)
     ! now we create the ylocal array, fill it, then associate the c pointer
     local(i)%ylength = 2 + local(i)%xlength
     ylocal => NULL()
     allocate(ylocal(local(i)%ylength))
     ! assign some data
     do j = 1, local(i)%ylength
        ylocal(j) = 0.5 + j
     end do
     ! set the c pointer to the correct location
     local(i)%y = c_loc(ylocal(1))
     ! print the y data
     write (*, '(a i3)') "    f ylength : ", local(i)%ylength
     do j = 1, local(i)%ylength
        write (*, '     (f6.2)', advance='no') ylocal(j)
     enddo
     write(*, *)
     ! now we write the output string
     ysum = sum(ylocal)
     write (outstr, '(a a f6.2)') trim(name), " : ", ysum
     write (*, '(a a)') "    outstr : ", trim(outstr)
     call F_C_string_ptr(trim(outstr), local(i)%outstr, MAX_STRING_LENGTH)
     write (*, '(a a a)') "    outstr : '", local(i)%outstr, "'"
  end do

end subroutine process_container

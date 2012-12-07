subroutine process_data(d)
  use, intrinsic :: iso_c_binding
  implicit none
  type, bind(c) :: data
     integer (c_int) :: xlength
     type (c_ptr) :: x
     integer (c_int) :: ylength
     type (c_ptr) :: y
  end type data
  type (data), intent(inout) :: d
  integer :: i
  real (c_double), allocatable, target, save :: ylocal(:)
  real (c_float), pointer :: x_array(:)

  ! Associate x_array with an array allocated in c
  call c_f_pointer (d%x, x_array, (/d%xlength/) ) 

  write (*, '(a i3)') "f x_length : ", d%xlength
  do i = 1, d%xlength
     write (*, '(f6.2)', advance='no') x_array(i)
  enddo
  write(*, *)

  ! Allocate a fortran array and make it available to C
  d%ylength = 2*d%xlength
  allocate (ylocal(d%ylength))
  do i = 1, d%ylength
     ylocal(i) = 0.5 * i
  enddo

  write (*, '(a i3)') "f y_length : ", d%ylength
  do i = 1, d%ylength
     write (*, '(f6.2)', advance='no') ylocal(i)
  enddo
  write(*, *)

  d%y = c_loc(ylocal)
end subroutine process_data

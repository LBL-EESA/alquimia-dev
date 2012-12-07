subroutine negate(p)
  type point
     integer(8) :: w
     real(4) :: x
     integer(4) :: y
     real(8) :: z
  end type point
  type (point) :: p
  p%w = -p%w
  p%x = -p%x
  p%y = -(p%y + 1)
  p%z = -p%z
end subroutine negate

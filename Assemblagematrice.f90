MODULE assemblagematrice 

    IMPLICIT NONE 

    CONTAINS 

    SUBROUTINE  GRADn (x1,y1, x2, y2, x3, y3, Dx, Dy)
       REAL(KIND=8), INTENT(in)  :: x1,y1, x2, y2, x3, y3 ! input parameters
       REAL(kind=8), INTENT(out) :: Dx, Dy ! output parameters
       REAL(kind=8)              :: D ! local variable
          D = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
         write(*,*)D
          Dx = (y2 - y3) / D
          Dy = (x3 - x2) / D
    END SUBROUTINE GRADn

   subroutine STIFF_DIAGn(n, x, y, sii )
integer, intent(in) :: n ! input parameters
real(kind=8), dimension(0:n),intent(in)::x,y ! input array
real(kind=8), intent(out) :: sii ! output parameters
real(kind=8) :: D, Dx, Dy ! local variables
integer :: j, j1 ! local variables
sii = 0.
do j=1,n
j1 = mod(j+1, n) + 1
call GRADn(x(0), y(0), x(j), y(j), x(j1), y(j1), Dx, Dy)
D = x(0)*(y(j)-y(j1) ) + x(j)*(y(j1)-y(0) ) + x(j1)*(y(0)-y(j) )
sii = sii + (Dx*Dx + Dy*Dy)*D/2 ! D/2 = area of triangle
end do
end subroutine STIFF_DIAGn

END MODULE assemblagematrice
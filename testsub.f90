program testsub
use assemblagematrice
    implicit none
    integer,parameter       :: n=10
    integer       ::i,j
    real(kind=8)  ::sii,d1,d2
    real(kind=8),dimension(:,:),allocatable::a
    real(kind=8),dimension(:),allocatable  ::x,y

    allocate(a(0:n,0:n),x(0:n),y(0:n))


    do i=0,n
        x(i)=i
        y(i)=i**2
    end do 
    a=0.
   write(*,*)x(1),y(1)
   call GRADn(x(0),y(0),x(1),y(1),x(2),y(2),d1,d2)
     call STIFF_DIAGn (n,x,y,sii)
    write(*,*)d1,d2




end program testsub
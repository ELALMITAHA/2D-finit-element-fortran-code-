  subroutine gradientconj(n,kmax,eps,b,u,A)
        implicit none
        real, dimension(1:n,1:n), intent(in) :: A
        real, dimension(1:n), intent(inout) :: u
        real, dimension(1:n), intent(in) :: b       
        real, intent(in) :: eps
        integer, intent(in) :: n,kmax
        real, dimension(1:n) :: d,r,w,r1
        real :: alpha, beta
        integer :: l    
 
        l=0
        r=matmul(A,u)-b
        d=r
        do while ( 0<=l .and. l<=kmax .and. sqrt(dot_product(r,r))>eps)
            w=(matmul(A,d))
            alpha=(dot_product(d,r))/(dot_product(d,w))
            u=u-alpha*d
            r1=r-alpha*w
            beta=(dot_product(r1,r1))/(dot_product(r,r))
            d=r1+beta*d
            l=l+1
            r=r1               
        end do
   end subroutine gradientconj
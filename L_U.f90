module L_u
    implicit none 

contains 

subroutine decente_remonter (A,n,x,b)
    
    integer,intent(in)::n 
    real(8),dimension(n,n),intent(in)::A
    real(8),dimension(n)  ,intent(in)::b
    real(8),dimension(n)  ,intent(out)::x
    real(8),dimension(n,n)            ::L,U
    real(8),dimension(n)              ::y,q,z
    integer                        ::i,k  
    real(8)                           ::s,t   

    call  LU(A,n,U,L)
    t=0.
    y(1)=b(1)/L(1,1)

    do i=2,n
        t=0.
	    do k=1,i-1
            t=t+L(i,k)*y(k)
	    end do 
    y(i)=(b(i)-t)/L(i,i)
    end do
    x(n)=y(n)/U(n,n)
    do i=n-1,1,-1
        s=0.
    	do k=i+1,n
            s=s+U(i,k)*x(k)
    	end do 
        x(i)=(y(i)-s)/U(i,i)
    end do 
    
end subroutine decente_remonter

subroutine LU(A,n,U,L)

  
   integer,intent(in)             ::n 
   real(8),dimension(n,n),intent(in) ::A
   real(8),dimension(n,n),intent(out)::L,U
   integer                        :: i ,k,j
   L=0.
    do i=1,n
        L(i,i)=1.
    end do 
    U=A
    do k=1,n-1
         if (U(k,k)==0)then 
            write(*,*)"la matrice n'admet pas de dÃ©composition LU"
        end if 
         do i=k+1,n
            L(i,k)=U(i,k)/U(k,k)
            do j=k+1,n
                U(i,j)=U(i,j)-L(i,k)*U(k,j)
            end do 
            U(i,k)=0
        end do 
    end do                                         
end subroutine LU

end module   L_U

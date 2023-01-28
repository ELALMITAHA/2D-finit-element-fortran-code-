program  stockage_matrices
!stockage crs d'une matrice A 
!--------------------------------------
implicit none 
 
integer,parameter::n=5
integer  ::i,j,z,k,c
integer,dimension(:),allocatable::aa,ja,b
integer,dimension(6)::ia
integer,dimension(5,5)::A
integer,dimension(5)::x,y

!  write(*,*)"entrez la valeur du vecteur x"
  
! do i =1,5

! read(*,*)x(i)
! end do 
A(1,1)=1;A(1,2)=0;A(1,3)=0;A(1,4)=2;A(1,5)=0
A(2,1)=3;A(2,2)=4;A(2,3)=0;A(2,4)=5;A(2,5)=0
A(3,1)=6;A(3,2)=0;A(3,3)=7;A(3,4)=8;A(3,5)=9
A(4,1)=0;A(4,2)=0;A(4,3)=10;A(4,4)=11;A(4,5)=0
A(5,1)=0;A(5,2)=0;A(5,3)=0;A(5,4)=0;A(5,5)=12
 
DO i=1,5
  write(*,*)


 z=0   
       do i=1,n
          do j=1,n
             if (A(i,j)/=0)then
             z=z+1
             end if       
       end do 
           end do 
allocate(b(1:n))
b=0       
 do i=1,n
          do j=1,n
             if (A(i,j)/=0)then
             b(i)=b(i)+1
             end if       
       end do 
           end do 
 allocate(aa(1:z),ja(1:z))   
       
k=1
                do i=1,n
                 do j=1,n
            if (A(i,j)/=0)then 
                aa(k)=A(i,j)
                    k=k+1
                    
            end if 
                 end do 
                        end do 


 ia(1)=1

          do i =2,n+1
      ia(i)=ia(i-1)+b(i-1)
          end do 
k=1                   
                    do i=1,n
                        do j=1,n
                                          
                                if( A(i,j)/=0)then
                            ja(k)=j
                                k=k+1
                             end if            
                         end do 

                    end do 
 
write(*,*)"les valuers non nul  de la matrice A est aa :",aa
! write(*,*)"la somme des élément non nul sur chaque ligne est ia :",ia
! write(*,*)"les indices de colone de chaque élément non nul sont ja :",ja
! call nn(2,1,5,A,c)
! call mat_vect(A,x,y)
! write(*,*)"le vecteur produit de a et x est:",y
! write(*,*)"la veuleur du i,j non nul est :",c

contains 
!subroutine qui renvoie une le couple (i,j) connue le coefficient A(i,j) si il est non nul
!--------------------------------------------------------------------------------------------
subroutine nn(i,j,n,A,c)


implicit none

integer,intent(in)::i,j,n
integer,dimension(n,n),intent(in)::A
integer,intent(out)::c
integer ::p
          do p=ia(i),ia(i+1)-1
       if(ja(p)==j)then
       c=aa(p)
        exit
     end if 
         end do 

end subroutine nn 
!subroutine qui calcule le produit d un vaceteur et une matrice sans tenir compte des 0
!-----------------------------------------------------------------------------------------
subroutine mat_vect(A,x,y)
integer,dimension(n,n),intent(in)::A
integer,dimension(n),intent(in)::x
integer,dimension(n),intent(out)::y
y=0
         do i=1,n
           do k=ia(i),ia(i+1)-1
            y(i)=y(i)+aa(k)*x(ja(k))
          end do 
        end do
end subroutine mat_vect                         

end program stockage_matrices

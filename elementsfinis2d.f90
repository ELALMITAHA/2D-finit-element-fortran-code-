program ELEMENTSFINIS2D
!===========================================================================!
!****************Chargement des module**************************************!
!===========================================================================!
use L_U
!===========================================================================!
!PROGRAMME RESOLUTION DU LAPLACIEN 2D AVEC CONDITION DE BORDS dIRICHLET     !
!NON HOMOGENE PAR EL ALMI TAHA DECEMENBRE 2017                              !
!===========================================================================! 
!Pour débeugué je resous le problème  $-\DELTA uex=fs$
!                                     $uex=0 sur le bord de [0,1]x[0,1]     ! 
!dans un premier temps je n'utilise pas de stockage SCR la matrice étant pas 
!tres grande 
 implicit none
!===========================================================================!
!***********************Défintion des variables*****************************!
!===========================================================================!
!***************Liste chainé pour stocker le sommet*************************!
!===========================================================================!
  TYPE cellule
    INTEGER                                   ::js ! Contenu de chaque cellule
    TYPE(cellule),POINTER::suiv
    END TYPE cellule
    TYPE(cellule),DIMENSION(:),POINTER        ::HashTab
    TYPE(cellule),POINTER                     ::Newcel,Ptcell,PtCellpred
 !========================================================================!
  !************Variables de défintion de maillage**************************!
  !========================================================================!
  real(8)                                     ::s
  INTEGER                                     :: i,j,Ns, itemp,temp1,temp2,temp3,temp4,Nt,tmp,nbrtrint
  REAL(8), DIMENSION(:,:), ALLOCATABLE        :: coord,A,aloc
  INTEGER,DIMENSION(:,:),ALLOCATABLE          ::Nu,Nubo,Trianglint
  INTEGER,DIMENSION(:),ALLOCATABLE            ::Nsupvois,Jposi,JvCell,IndPL
  REAL(8),DIMENSION(:),ALLOCATABLE            ::Tmat                                
  INTEGER                                     ::is,jj,ii,k,next,is1,is2,js,iv,jv,kv
  INTEGER                                     ::ismin,ismax,ntypn,Nsegement,jt,info 
  REAL(kind=8),DIMENSION(:),ALLOCATABLE       :: U,F,Floc,w1,w2,w3,w4  
  INTEGER                                     ::Ncoefmat
 !============================================================================!
 !********************Lecture du maillage générer par GMSH====================!  
 !Ns      : Le nombre des noeuds du maillage                                  !
 !coord   : Un tableau de dimension 2 qui contient cood(x,y) de chaque noeuds !
 !Nt      : Nombre des traingles du maillage  (meme ceux qui sont sur le bord )!
 !============================================================================! 
  i=1
  OPEN(12,file="maillage.msh") ! ouverture du fichier
  DO WHILE (i<53) 
     READ(12,*)
     i=i+1
   END DO 
  Ns=108
   ALLOCATE(coord(1:2,1:Ns))
  DO i=1,Ns
     read(12,*) itemp, coord(1,i), coord(2,i) 
  END DO 
CLOSE(12)

!========================================================================!
!*************Lecture de la table des triangle***************************!
!========================================================================!
 Nt=254! Nombre totale des triangles et nbrtrint est le nombre des trinagle qui sont a l interieur
OPEN(7,file="triangle.msh")
ALLOCATE(Nu(1:3,1:Nt))
nbrtrint=0
DO i=1,nt
    read(7,*)itemp,temp1,temp2,temp3,temp4,Nu(1,i),Nu(2,i),Nu(3,i) ! contient la liste des triangle 
END DO 
CLOSE(7)
j=0
do i=1,254
   ! write(*,*) Nu(1,i),Nu(2,i),Nu(3,i)
    if(Nu(1,i)>40 .and.Nu(2,i)>40 .and.Nu(3,i)>40 )then 
        nbrtrint=nbrtrint+1
    end if 
end do 

ALLOCATE(Trianglint(1:3,1:nbrtrint))
j=0
do i=1,254
if(Nu(1,i)>40 .and.Nu(2,i)>40 .and.Nu(3,i)>40 )then 
j=j+1
Trianglint(1:3,j) =nu(1:3,i)
end if 
end do 
 !=============================================================================!
! !***************Assemblage de la matrice**************************************!
! !=============================================================================!
ALLOCATE(A(1:ns,1:ns),aloc(1:3,1:3),u(1:Ns))

DO k=1,nbrtrint
    call matricelocal(aloc,k)
  do i=1,3
    do j=1,3
       ii=Trianglint(i,k)-40
       jj=Trianglint(j,k)-40
       a(ii,jj)=a(ii,jj)+aloc(i,j)
   end do 
end do 
end do 


!===========================================================================!
!********************Assemblage second momebre******************************!
!===========================================================================!
  ALLOCATE(F(1:Ns),Floc(1:3),w1(1:2),w2(1:2),w3(1:2),w4(1:2))
 ! w1(1)=1/3;w1(2)=1/3;w2(1)=0.2;w2(2)=0.6;w3(1)=0.2;w3(2)=0.2;w4(1)=0.6;w4(2)=0.2 !De Gauss
 !  DO k=1,nbrtrint     !Nt tous les triangles meme ceux qui touche le bord 
 !  call airtriangle(k,s)
 !    floc(1)=((-27./96.)*lambda1(w1,k)+(25./96.)*lambda1(w2,k)+(25./96.)*(lambda1(w3,k)+lambda1(w4,k)))*2*s
 !    Floc(2)=((-27./96.)*lambda2(w1,k)+(25./96.)*lambda2(w2,k)+(25./96.)*(lambda2(w3,k)+lambda2(w4,k)))*2*s
 !    Floc(3)=((-27./96.)*lambda3(w1,k)+(25./96.)*lambda3(w2,k)+(25./96.)*(lambda3(w3,k)+lambda3(w4,k)))*2*s
 !     Do i=1,3
 !        ii=Trianglint(i,k)-40
 !        F(ii)=F(ii)+floc(i)
 !   END DO 
 !  END DO 
f=0.

call decente_remonter(A,Ns,U,F)
open(9,file="fichier.vtk")
do i=1,Ns
    write(9,*) coord(1,i),coord(2,i),uex(coord(1,i),coord(1,2))
end do 




!=============================================================================!
 Contains 
!=============================================================================!

!=============================================================================!
!********************Fonction second Membre **********************************!
!=============================================================================!
REAL FUNCTION fs(x,y)
  REAL(8),INTENT(In):: x, y
  Fs=(y-y**2+x-x**2)*2
END FUNCTION
!=============================================================================!
!*********************Transformation affine***********************************!
!=============================================================================!
SUBROUTINE G(x,y,k)
  REAL(8),DIMENSION(1:2),INTENT(In)::x
  REAL(8),DIMENSION(1:2),INTENT(out)::y
  INTEGER                        ::k   
 y(1)=(coord(1,Trianglint(2,k)-40)-coord(1,Trianglint(1,k)-40))*x(1)+(coord(1,Trianglint(3,k)-40)&
    -coord(1,Trianglint(1,k))-40)*x(2)+coord(1,Trianglint(1,k)-40)+coord(1,Trianglint(1,k)-40)
 y(2)=(coord(2,Trianglint(2,k)-40)-coord(2,Trianglint(1,k)-40))*x(1)+(coord(2,Trianglint(3,k)-40)&
    -coord(2,Trianglint(1,k))-40)*x(2)+coord(2,Trianglint(1,k)-40)+coord(2,Trianglint(1,k)-40)

END SUBROUTINE G
!============================================================================!
!****************Fonctions de base element de référence**********************!
!============================================================================!
REAL FUNCTION lambda1(x,k)
  INTEGER,INTENT(in)::k
  REAL(8),DIMENSION(1:2),INTENT(IN)::x 
  REAL(8),DIMENSION(1:2)::vect
  call G(x,vect,k)
  lambda1=(1-x(1)-x(2))*Fs(vect(1),vect(2))
END FUNCTION lambda1
!===========================================================================!
REAL FUNCTION lambda2(x,k)
  INTEGER,INTENT(in)::k
  REAL(8),DIMENSION(1:2),INTENT(IN)::x
  REAL(8),DIMENSION(1:2)::vect
  call G(x,vect,k)
  lambda2=x(1)*Fs(vect(1),vect(2))
END FUNCTION lambda2
!===========================================================================!
REAL FUNCTION lambda3(x,k)
  INTEGER,INTENT(in)::k
  REAL(8),DIMENSION(1:2),INTENT(IN)::x
  REAL(8),DIMENSION(1:2)::vect
 call G(x,vect,k)
  lambda3=x(2)*Fs(vect(1),vect(2))

END FUNCTION lambda3
!===========================================================================
!*********************Solution exacte**************************************!
!===========================================================================
REAL FUNCTION uex(x,y)
  REAL(8),INTENT(IN)::x,y
  uex=x*y*(1-x)*(1-y)
END FUNCTION
!=============================================================================
!**************Caclule de l'aire du triagle k*********************************
!=============================================================================
SUBROUTINE airtriangle(k,s)
  INTEGER,INTENT(IN) ::k
  REAL(kind=8),INTENT(OUT)   ::s
  s=0.5*((coord(1,Trianglint(2,k)-40)-coord(1,Trianglint(1,k)-40))*(coord(2,Trianglint(3,k)-40)-coord(2,Trianglint(1,k)-40))-&
      (coord(2,Trianglint(2,k)-40)-coord(2,Trianglint(1,k)-40))*(coord(1,Trianglint(3,k)-40)-coord(1,Trianglint(1,k)-40)))
END SUBROUTINE airtriangle
!============================================================================
!************************Cacul norme vercteur********************************
!============================================================================
SUBROUTINE Norme(x1,x2,y1,y2,r)
  REAL(8),INTENT(in)::x1,x2,y1,y2
  REAL(8),INTENT(out)::R
  r=(x1-y1)**2+(x2-y2)**2
END SUBROUTINE
!============================================================================
!*********************produit scalaire***************************************
!============================================================================
SUBROUTINE prdscal(x1,x2,y1,y2,p)
  REAL(8),INTENT(IN)::x1,x2,y1,y2
  REAL(8),INTENT(out)::P
  P=x1*y1+x2*y2
END SUBROUTINE prdscal
!============================================================================
!********************Coordonnée vecteur**************************************
!============================================================================
SUBROUTINE corvect(x1,x2,y1,y2,c)
  REAL(8),INTENT(in)             ::x1,x2,y1,y2
  REAL(8),DIMENSION(1:2),INTENT(out)::c
  c(1)=y1-x1
  c(2)=y2-x2
END SUBROUTINE corvect
!============================================================================
!*************************Calcule de la matricelocal*************************
!============================================================================
!Calcul l'la matrice elementaire sur un triagle k pour k=1,...Nt
SUBROUTINE matricelocal(M,k)   !
  INTEGER,INTENT(In)::k
  REAL(8),DIMENSION(1:3,1:3),INTENT(out)::M
  REAL(8)                               ::s
  REAL (8)                              ::p1,p2,p3,n1,n2,n3
  Real(8),DIMENSION(1:2)                ::c1,c2,c3                                   
  call airtriangle(k,s)  
  call corvect(coord(1,Trianglint(1,k)-40),coord(2,Trianglint(1,k)-40),coord(1,Trianglint(2,k)-40),coord(2,Trianglint(2,k)-40),c1)!Coordonnée de AB
  call corvect(coord(1,Trianglint(2,k)-40),coord(2,Trianglint(2,k)-40),coord(1,Trianglint(3,k)-40),coord(2,Trianglint(3,k)-40),c2)!Coordonnée de BC
  call corvect(coord(1,Trianglint(1,k)-40),coord(2,Trianglint(1,k)-40),coord(1,Trianglint(3,k)-40),coord(2,Trianglint(3,k)-40),c3)!Coordonnée de AC
  call Norme(coord(1,Trianglint(2,k)-40),coord(2,Trianglint(2,k)-40),coord(1,Trianglint(3,k)-40),coord(2,Trianglint(3,k)-40),n1)!Norme du vecteur BC 
  call Norme(coord(1,Trianglint(1,k)-40),coord(2,Trianglint(1,k)-40),coord(1,Trianglint(3,k)-40),coord(2,Trianglint(3,k)-40),n2)!Norme du vecteur AC
  call Norme(coord(1,Trianglint(1,k)-40),coord(2,Trianglint(1,k)-40),coord(1,Trianglint(2,k)-40),coord(2,Trianglint(2,k)-40),n3)!Norme du vecteur AB
  Call prdscal(c3(1),c3(2),c2(1),c2(2),p1)!AC.BC 
  Call prdscal(c1(1),c1(2),c3(1),c3(2),p2)!BA.CA
  Call prdscal(c1(1),c1(2),c2(1),c2(2),p3)!AB.CB
  M=0.
  m(1,1)=n1;m(1,2)=-p1;m(1,3)=-p2
  m(2,1)=-p1;m(2,2)=n2;m(2,3)=p3
  m(3,1)=-p2;m(3,2)=p3 ;m(3,3)=n3
  m=m/(4.*s)
  END SUBROUTINE


END PROGRAM ELEMENTSFINIS2D
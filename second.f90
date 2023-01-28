
  ALLOCATE(F(1:Ns),Floc(1:3),w1(1:2),w2(1:2),w3(1:2),w4(1:2))
 w1(1)=1/3;w1(2)=1/3;w2(1)=0.2;w2(2)=0.6;w3(1)=0.2;w3(2)=0.2;w4(1)=0.6;w4(2)=0.2 !De Gauss
  DO k=1,Nt     !Nt tous les triangles meme ceux qui touche le bord 
  call airtriangle(k,s)
    floc(1)=((-27./96.)*lambda1(w1,k)+(25./96.)*lambda1(w2,k)+(25./96.)*(lambda1(w3,k)+lambda1(w4,k)))*2*s
    Floc(2)=((-27./96.)*lambda2(w1,k)+(25./96.)*lambda2(w2,k)+(25./96.)*(lambda2(w3,k)+lambda2(w4,k)))*2*s
    Floc(3)=((-27./96.)*lambda3(w1,k)+(25./96.)*lambda3(w2,k)+(25./96.)*(lambda3(w3,k)+lambda3(w4,k)))*2*s
     Do i=1,3
        ii=Nu(i,k)
       IF(II>40)THEN 
       F(ii)=F(ii)+floc(i)
   END IF 
     END DO 
  END DO 
REAL FUNCTION fs(x,y)
  REAL(8),INTENT(In):: x, y
  Fs=(y-y**2-x-x**2)*2
END FUNCTION
!=============================================================================!
!*********************Transformation affine***********************************!
!=============================================================================!
SUBROUTINE G(x,y,k)
  REAL(8),DIMENSION(1:2),INTENT(In)::x
  REAL(8),DIMENSION(1:2),INTENT(out)::y
  INTEGER                        ::k   
 y(1)=(coord(1,Nu(2,k))-coord(1,Nu(1,k)))*x(1)+(coord(1,Nu(3,k))-coord(1,Nu(1,k)))*x(2)+coord(1,Nu(1,k))+coord(1,Nu(1,k))
 y(2)=(coord(2,Nu(2,k))-coord(2,Nu(1,k)))*x(1)+(coord(2,Nu(3,k))-coord(2,Nu(1,k)))*x(2)+coord(2,Nu(1,k))+coord(2,Nu(1,k))

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
  s=0.5*((coord(1,Nu(2,k))-coord(1,Nu(1,k)))*(coord(2,Nu(3,k))-coord(2,Nu(1,k)))-&
      (coord(2,Nu(2,k))-coord(2,Nu(1,k)))*(coord(1,Nu(3,k))-coord(1,Nu(1,k))))
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
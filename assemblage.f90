Module assemblage 
 contains 

!========================================================================!
! ***************Construction de la liste des segement********************!
!========================================================================!
!cette partie est bien vérifié elle fonctionne bien                       !
!=========================================================================!
! Ns=SIZE(coord,DIM=2)
! ALLOCATE(HashTab(1:Ns))
! DO is =1,Ns
!   HashTab(is)%js=is
!   NULLIFY(HashTab(is)%suiv)
! END DO 
! ALLOCATE(Nsupvois(1:Ns))
! Nsupvois(1:Ns)=0
! Nsegement=0

!  NULLIFY(Ptcell,PtCellpred)
!  ntypn=3

! DO jt=1,Nt
!   do k=1,3
!     next=MOD(k,ntypn)+1
!     is1=Nu(k,jt)
!     is2=Nu(next,jt)
!     ismin=MIN(is1,is2);ismax=MAX(is1,is2)
!     Ptcell=>HashTab(ismin)%suiv
!     PtCellpred=>HashTab(ismin)
!     Do 
!       IF(.NOT.ASSOCIATED(Ptcell))EXIT
!       IF(Ptcell%js.GT.ismax)EXIT
!      PtCellpred=>PtCell
!       Ptcell=>PtCell%suiv
!     END DO 

!     IF (PtCellpred%js.LT.ismax)THEN
!       ALLOCATE(Newcel)
!       Newcel%js=ismax
!       Newcel%suiv=>Ptcell
!       PtCellpred%suiv=>Newcel
!       Nsupvois(ismin)=Nsupvois(ismin)+1
!       Nsegement=Nsegement+1
!     END IF
!   END DO 
! END DO  
! ALLOCATE(Nubo(2,Nsegement))
! k=0

! DO is=1,Ns
!   IF(Nsupvois(is).NE.0)THEN
!     PtCell=>HashTab(is)%suiv
!     DO iv=1,Nsupvois(is)
!       k=k+1
!       js=Ptcell%js
!       Nubo(1,k)=is;Nubo(2,k)=js
!       Ptcell=>Ptcell%suiv
!     ENd DO 
!   END IF 
! END DO 
! DEALLOCATE(HashTab)
! !=============================================================================!
! !************************ Stockage Morse*
! !=============================================================================!
!   !Calcul de Ncoefmat
!   Ncoefmat=Ns!    On est sur qu'au moins il ya Ns coef non dans la matrice car les teremes diagonaux sont non nul
!   DO i= 1,Nsegement
!     is=Nubo(1,i)
!     js=Nubo(2,i)
!     IF (is<=Ns.and.js<=Ns)THEN 
!       Ncoefmat=Ncoefmat+2
!     END IF 
!   END DO 
!   !Allocation des tableaux
!   ALLOCATE(Tmat(1:Ncoefmat))        ! contient le Nombre des éléments Non null de la matrice 
!   ALLOCATE(Jposi(1:Ns+1))           ! Indice dans le tableau Tmat du premier element non nul de la ligne i
!   ALLOCATE(JvCell(1:Ncoefmat))      !
!   !Construction de Jposi
!   DO i=1,Ns+1
!     Jposi(i)=i 
!   END DO  
! DO i=1,Nsegement
!   is=Nubo(1,i)
!   js=Nubo(2,i)
!   IF(is<=Ns.and.js<=Ns)then 
!     Jposi(is+1:Ns+1)=Jposi(is+1:Ns+1)+1
!     Jposi(js+1:Ns+1)=Jposi(js+1:Ns+1)+1
!   END IF 
! END DO 
! !Construction du tebaleau JvCell
! ALLOCATE(IndPL(1:Ns))
! IndPL(1:Ns)=Jposi(1:Ns)

! !Indices des colone des termes diagonaux

! DO i=1,Ns
!   JvCell(IndPL(i))=i
!   IndPL(i)=IndPL(i)+1
! END DO 
! ! Indices des termes extrats diagonaux
! DO i=1,Nsegement
!   is=Nubo(1,i)
!   js=Nubo(2,i)
!   if (is<=Ns .and. js<=Ns)THEN
!     JvCell(IndPL(is))=js
!     IndPL(is)=IndPL(is)+1
!     Jvcell(Indpl(js))=is
!       IndPL(js)=IndPL(js)+1
!     END IF 
!   END DO 

!   DO is=1,Ns
!     Do jv=Jposi(is+1)-1,Jposi(is),-1
!       do kv =Jposi(is)+1,jv
!         IF(Jvcell(kv-1)>Jvcell(kv))THEN 
!           tmp=Jvcell(kv-1)
!           Jvcell(kv-1)=Jvcell(kv)
!           Jvcell(kv)=tmp
!         END IF 
!       END DO 
!     END DO 
!   END DO 
 !============================================================================
!*********************Assemblage de matrice**********************************
!============================================================================
!SUBROUTINE Assemble(m,Aloc,Tmat,Jposi,JvCell)
!============================================================================
SUBROUTINE Ajout (ii,jj,Coef,TMat,Jposi,Jvcell)
  INTEGER ,INTENT(IN) ::ii,jj
  REAL(8) ,INTENT(IN) ::coef
  REAL(8) ,DIMENSION(1:Ncoefmat),INTENT(OUT)::Tmat
  INTEGER ,DIMENSION(1:Ns+1),    INTENT(OUT)::Jposi
  INTEGER ,DIMENSION(1:Ncoefmat),INTENT(OUT)::Jvcell
  INTEGER                                   ::i,j 
  LOGICAL                                   ::trouver 
  trouver=.False.
  do j=jposi(ii),Jposi(ii+1)-1
    IF (Jvcell(j)==jj)then
      Tmat(j)=Tmat(j)+Coef
      Trouver=.TRUE.
      EXIT
    END IF 
  END DO 
    IF(trouver.eqv..FALSE.)THEN
      write(*,*)"Probleme d'Assemblage de la matrice A"
      STOP 
    END IF
  END SUBROUTINE Ajout
  SUBROUTINE Assemble(m,Aloc,Tmat,jposi,Jvcell)
    INTEGER ,INTENT(IN)                        ::m
    REAL(8) ,DIMENSION(1:3,1:3),INTENT(IN)     ::Aloc    
    REAl(8) ,DIMENSION(1:Ncoefmat),INTENT(OUT)::Tmat
    INTEGER ,DIMENSION(1:Ns+1),    INTENT(OUT)::Jposi
    INTEGER ,DIMENSION(1:Ncoefmat),INTENT(OUT)::Jvcell
    INTEGER                                   ::i,j,k
    i=(1,m)
    j=Nu(2,m)
    k=Nu(3,m)
    IF (i<=Ns)THEN 
      CALL Ajout(i,i,aloc(1,1),Tmat,Jposi,Jvcell)
      IF (j<=Ns)THEN 
        CALL Ajout(i,j,aloc(1,2),TMat,Jposi,Jvcell)
      END IF 
        IF (k<=Ns)THEN 
        CALL Ajout(i,k,aloc(1,3),TMat,Jposi,Jvcell)
      END IF 
    END IF 
    IF (j<=Ns)THEN 
      CALL Ajout(j,j,aloc(2,2),Tmat,Jposi,Jvcell)
      IF (i<=Ns)THEN 
        CALL Ajout(j,i,aloc(2,1),TMat,Jposi,Jvcell)
      END IF 
        IF (k<=Ns)THEN 
        CALL Ajout(j,k,aloc(2,3),TMat,Jposi,Jvcell)
      END IF 
    END IF 
    IF (k<=Ns)THEN 
      CALL Ajout(k,k,aloc(3,3),Tmat,Jposi,Jvcell)
      IF (j<=Ns)THEN 
        CALL Ajout(k,j,aloc(3,2),TMat,Jposi,Jvcell)
      END IF 
        IF (i<=Ns)THEN 
        CALL Ajout(k,i,aloc(3,1),TMat,Jposi,Jvcell)
      END IF 
    END IF 
  END SUBROUTINE Assemble
!===================================================================!
!********************Assemblage second Membre***********************!
!*******************************************************************!

  end Module
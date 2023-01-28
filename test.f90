program lecture

  implicit none
  integer :: i, N, itemp
  real(8), dimension(:,:), allocatable :: coord
  
  i=1
  open(12,file="maillage.msh") ! ouverture du fichier

  ! on va passer les 12 premieres lignes car elles ne 
  ! nous interessent pas
  do while (i<13) 
     read(12,*)
     i=i+1
     print*,i
  end do

  ! on lit le nombre de noeuds du maillage
  read(12,*) N

  print*, N

  ! on peut alors allouer le tableau des coordonnees
  allocate(coord(1:2,N))

  ! on lit les coordonnees x et y de chaque noeud
  ! itemp sert juste a stocker le numero du noeud pour 
  ! pouvoir lire ses coordonnees dans les colonnes suivantes
  do i=1,N
     read(12,*) itemp, coord(1,i), coord(2,i) 
  end do

  ! Verification :
  do i=1,N
     print*, i, coord(1,i), coord(2,i)
  end do

  close(12)


end program lecture
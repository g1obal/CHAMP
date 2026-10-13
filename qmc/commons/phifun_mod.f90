module phifun_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: phin(:,:), dphin(:,:,:), d2phin(:,:)
 
 ! 2D HESSIAN ARRAYS
 double precision, allocatable :: d2phin_xx(:,:), d2phin_yy(:,:), d2phin_xy(:,:)

end module phifun_mod

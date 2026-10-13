module d2orb_mod

 use constants_mod
 implicit none
 save

 ! 2D HESSIAN ARRAYS
 double precision, allocatable :: d2orb_xx(:,:), d2orb_yy(:,:), d2orb_xy(:,:)

end module d2orb_mod

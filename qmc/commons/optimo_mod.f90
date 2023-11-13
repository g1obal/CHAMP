module optimo_mod

 use constants_mod
 implicit none
 save

 integer nparmot,notype,iconstrain_gauss_orbs
 integer, allocatable :: iwo(:,:), nparmo(:)
 integer, allocatable :: norb_constraints(:), orb_constraints(:,:,:)
 real*8 oparm3_max


end module optimo_mod

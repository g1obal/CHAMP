module pairden_mod

 implicit none
 save
 
 double precision, allocatable :: xx0probut(:,:,:),xx0probuu(:,:,:),xx0probud(:,:,:),xx0probdt(:,:,:)
 double precision, allocatable :: xx0probdu(:,:,:),xx0probdd(:,:,:),den2d_t(:,:),den2d_d(:,:),den2d_u(:,:)
 double precision, allocatable :: pot_ee2d_t(:,:), pot_ee2d_u(:,:), pot_ee2d_d(:,:)       
 double precision, allocatable :: dos(:)
 double precision, allocatable :: thetafix(:) !GO
 integer, allocatable          :: imeshfix1(:), imeshfix2(:) !GO
 double precision              :: dos_dele
 double precision              :: delxi(3),xmax,xfix(3)
 integer ithetafix !GO
 integer ifixe
 integer, parameter :: NAX = 50
 
! Variables for M-fold orientational order parameter
 integer :: M_pd, nrings_pd
 integer, dimension(20) :: conf_pd
 
! Variables for Inter-Ring Phase Correlation
 integer, parameter :: NIRBINS_pd = 360
 double precision, allocatable :: irphase_pd(:)

end module pairden_mod

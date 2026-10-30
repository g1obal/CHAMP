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
 double precision              :: pair_hits_u, pair_hits_d  ! Exact trigger trackers
 integer ithetafix !GO
 integer ifixe
 integer, parameter :: NAX = 50
 
! Select a radial shell by index; M_bff is its derived population, not an input.
 integer :: ring_id_bff, M_bff, nrings_bff
 integer, dimension(20) :: conf_bff
 double precision, allocatable :: bffden2d_t(:,:), bffden2d_u(:,:), bffden2d_d(:,:)
 double precision :: psi_M_acc_bff(20)
 
! Variables for Inter-Ring Phase Correlation
 integer :: M_ir_bff ! Common harmonic of the two outermost shells.
 integer, parameter :: NIRBINS_bff = 100
 double precision, allocatable :: irphase_bff(:)

end module pairden_mod

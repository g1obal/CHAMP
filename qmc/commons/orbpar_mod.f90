module orbpar_mod

 implicit none
 save

 double precision, allocatable :: oparm(:,:,:)
 double precision, allocatable :: oparm_sav(:,:)
 double precision, allocatable :: oparm_best(:,:)
 double precision :: Cdamp_pgs

 ! Fixed Gaussian setup, prepared by initialize_gaussian_basis.
 logical, allocatable :: gauss_central(:)
 double precision :: gauss_scale_r

end module orbpar_mod

module orbpar_mod

 implicit none
 save

 double precision, allocatable :: oparm(:,:,:)
 double precision, allocatable :: oparm_sav(:,:)
 double precision, allocatable :: oparm_best(:,:)
 double precision :: Cdamp_pgs

end module orbpar_mod

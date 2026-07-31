module coefs_mod

 use constants_mod

 integer               :: nbasis,norb
 real(dp), allocatable :: coef(:,:,:)
 logical               :: coef_is_diag
 integer               :: morb = 0 !JT

end module coefs_mod

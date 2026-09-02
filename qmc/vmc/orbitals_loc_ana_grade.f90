      subroutine orbitals_loc_ana_grade(iel,rvec_en,r_en,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! Calculate localized orbitals and derivatives for all or 1 electrons
      use control_mod
      use coefs_mod
      use dim_mod
      use wfsec_mod
      use contrl_per_mod
      use phifun_mod
      use const_mod
      use atom_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent) &
     &,orb(norb),dorb(3,norb),ddorb(norb)

! get basis functions
      if(ndim.eq.3) then
        call basis_fns(iel,rvec_en,r_en)
       elseif(ndim.eq.2) then
         if(ibasis.eq.1) then
           call basis_fns_2d(iel,rvec_en,r_en)
         elseif(ibasis.eq.4) then
           call basis_fns_2dgauss(iel,rvec_en,r_en)
         elseif(ibasis.eq.5) then
           call basis_fns_polargauss(iel,rvec_en,r_en)
         elseif(ibasis.eq.6) then
           call basis_fns_2dgauss_noncirc(iel, rvec_en,r_en)
         elseif(ibasis.eq.7) then
           if(iper_gaussian_type.eq.1) then
             call basis_fns_2dgauss_periodic(iel,rvec_en,r_en)
           else
             call basis_fns_2dgauss_periodic2(iel,rvec_en,r_en)
           endif
         elseif(ibasis.eq.8) then
           call basis_fns_polargauss_damped(iel,rvec_en,r_en)
         else
           stop 'orbitals_loc_ana: ibasis must be 1,4,5, 6, 7 or 8 for 2d systems'
         endif
      endif

! If ibasis.eq.4 or 5 and coef is a multiple of the unit matrix
      if((ibasis.eq.4 .or. ibasis.eq.5 .or. ibasis.eq.8) .and. coef_is_diag) then
        ! >>> MASSIVE O(N) SPEEDUP FOR DIAGONAL ORBITALS <<<
        do 24 iorb=1,norb
          c_val = coef(iorb,iorb,iwf)
          orb(iorb)=c_val*phin(iorb,iel)
          dorb(1,iorb)=c_val*dphin(1,iorb,iel)
          dorb(2,iorb)=c_val*dphin(2,iorb,iel)
          dorb(3,iorb)=c_val*dphin(3,iorb,iel)
   24     ddorb(iorb)=c_val*d2phin(iorb,iel)
      else
        ! >>> EXACT, CACHE-OPTIMIZED O(N^2) FALLBACK FOR LCAO <<<
        orb = 0.0d0
        dorb = 0.0d0
        ddorb = 0.0d0
        do iorb=1,norb
          do m=1,nbasis
            c_val = coef(m,iorb,iwf)
            orb(iorb)=orb(iorb)+c_val*phin(m,iel)
            dorb(1,iorb)=dorb(1,iorb)+c_val*dphin(1,m,iel)
            dorb(2,iorb)=dorb(2,iorb)+c_val*dphin(2,m,iel)
            dorb(3,iorb)=dorb(3,iorb)+c_val*dphin(3,m,iel)
            ddorb(iorb)=ddorb(iorb)+c_val*d2phin(m,iel)
          enddo
        enddo
      endif

      return
      end

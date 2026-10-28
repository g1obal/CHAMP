      subroutine orbitals_loc_ana(iel,rvec_en,r_en,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! Calculate localized orbitals and derivatives for all or 1 electrons
      use all_tools_mod
      use control_mod
      use coefs_mod
      use const_mod
      use dim_mod
      use wfsec_mod
      use contrl_per_mod
      use phifun_mod
      use atom_mod
      use d2orb_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent) &
     &,orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb)

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

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
           call basis_fns_2dgauss_noncirc(iel,rvec_en,r_en)
         elseif(ibasis.eq.7) then
           if(iper_gaussian_type.eq.1) then
             call basis_fns_2dgauss_periodic(iel,rvec_en,r_en)
           else
             call basis_fns_2dgauss_periodic2(iel,rvec_en,r_en)
           endif
         elseif(ibasis.eq.8) then
           call basis_fns_polargauss_damped(iel,rvec_en,r_en)
         else
           stop 'orbitals_loc_ana: ibasis must be 1,4,5,6,7 or 8 for 2d systems'
         endif
      endif

      if(ndim.eq.2) dorb(3,:,:)=0.d0

! If ibasis.eq.4 or 5 and coef is a multiple of the unit matrix
      if((ibasis.eq.4 .or. ibasis.eq.5 .or. ibasis.eq.8) .and. coef_is_diag) then
        ! Diagonal coefficients require only a scaled copy.
        do 24 iorb=1,norb
          c_val = coef(iorb,iorb,iwf)
          do 24 ie=nelec1,nelec2
            orb(ie,iorb) = c_val*phin(iorb,ie)
            dorb(1,ie,iorb) = c_val*dphin(1,iorb,ie)
            dorb(2,ie,iorb) = c_val*dphin(2,iorb,ie)
            if(ndim.eq.3) dorb(3,ie,iorb) = c_val*dphin(3,iorb,ie)
            ddorb(ie,iorb) = c_val*d2phin(iorb,ie)
            
            d2orb_xx(ie,iorb) = c_val*d2phin_xx(iorb,ie)
            d2orb_yy(ie,iorb) = c_val*d2phin_yy(iorb,ie)
   24       d2orb_xy(ie,iorb) = c_val*d2phin_xy(iorb,ie)
      else 
        ! Accumulate the dense coefficient expansion.
        orb(nelec1:nelec2,:) = 0.0d0
        dorb(:,nelec1:nelec2,:) = 0.0d0
        ddorb(nelec1:nelec2,:) = 0.0d0
        
        d2orb_xx(nelec1:nelec2,:) = 0.0d0
        d2orb_yy(nelec1:nelec2,:) = 0.0d0
        d2orb_xy(nelec1:nelec2,:) = 0.0d0
        
        do iorb=1,norb
          do m=1,nbasis
            c_val = coef(m,iorb,iwf)
            do ie=nelec1,nelec2
              orb(ie,iorb) = orb(ie,iorb) + c_val*phin(m,ie)
              dorb(1,ie,iorb) = dorb(1,ie,iorb) + c_val*dphin(1,m,ie)
              dorb(2,ie,iorb) = dorb(2,ie,iorb) + c_val*dphin(2,m,ie)
              if(ndim.eq.3) dorb(3,ie,iorb) = dorb(3,ie,iorb) + c_val*dphin(3,m,ie)
              ddorb(ie,iorb) = ddorb(ie,iorb) + c_val*d2phin(m,ie)
              
              d2orb_xx(ie,iorb) = d2orb_xx(ie,iorb) + c_val*d2phin_xx(m,ie)
              d2orb_yy(ie,iorb) = d2orb_yy(ie,iorb) + c_val*d2phin_yy(m,ie)
              d2orb_xy(ie,iorb) = d2orb_xy(ie,iorb) + c_val*d2phin_xy(m,ie)
            enddo
          enddo
        enddo
      endif

!     do  iorb=1,norb
!       do  ie=nelec1,nelec2
!         write(6,*) 'orb(ie,iorb),ie,iorb=',orb(ie,iorb),ie,iorb
!       enddo
!     enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine orbitals_loc_anae(iel,rvec_en,r_en,orb)
! Written by Cyrus Umrigar
! Calculate localized orbitals for electron iel
      use coefs_mod
      use dim_mod
      use wfsec_mod
      use phifun_mod
      use atom_mod
      use const_mod
      use contrl_per_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent) &
     &,orb(norb)

! get basis functions
      if(ndim.eq.3) then
        call basis_fnse2(iel,rvec_en,r_en)
       elseif(ndim.eq.2) then
        if(ibasis.eq.4) then
          call basis_fns_2dgauss(iel,rvec_en,r_en)
        elseif(ibasis.eq.5) then
          call basis_fns_polargauss(iel,rvec_en,r_en)
        elseif(ibasis.eq.8) then
          call basis_fns_polargauss_damped(iel,rvec_en,r_en)
        else
          call basis_fns_2de2(iel,rvec_en,r_en)
        endif
      endif

      if((ibasis.eq.4 .or. ibasis.eq.5 .or. ibasis.eq.8) .and. coef_is_diag) then !GO
        ! DIAGONAL SHORTCUT
        do iorb=1,norb
          orb(iorb)=coef(iorb,iorb,iwf)*phin(iorb,iel)
        enddo
      else
        orb = 0.0d0
        do iorb=1,norb
          do m=1,nbasis
            orb(iorb)=orb(iorb)+coef(m,iorb,iwf)*phin(m,iel)
          enddo
        enddo
      endif
      
      return
      end

!---------------------------------------------------------------------------

      subroutine deriv_orbitals(rvec_en,r_en,orb,dorb,ddorb)
! Written by A.D.Guclu (Apr 2005) starting from orbitals_loc_ana.f
! Modified by Gokhan Oztarhan (Oct 2023)
! Calculate localized orbitals, coo. and parameter derivatives for all electrons
      use control_mod
      use coefs_mod
      use const_mod
      use dim_mod
      use wfsec_mod
      use contrl_per_mod
      use phifun_mod
      use optimo_mod
      use atom_mod
      use deriv_phifun_mod
      use d2orb_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent)
      dimension orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb)

      if(ndim.eq.2) dorb(3,:,:)=0.d0

      nelec1=1
      nelec2=nelec

! get basis functions
      if(ndim.ne.2) stop 'deriv_orbitals: ndim must be 2'
      if(ibasis.eq.4) then
        call deriv_2dgauss(rvec_en,r_en)
      elseif(ibasis.eq.5) then
        call deriv_polargauss(rvec_en,r_en)
      elseif(ibasis.eq.6) then
        call deriv_2dgauss_noncirc(rvec_en,r_en)
      elseif(ibasis.eq.7) then
        if (iper_gaussian_type.eq.1) then
          call deriv_2dgauss_periodic(rvec_en,r_en)
        else
          call deriv_2dgauss_periodic2(rvec_en,r_en)
        endif
      elseif(ibasis.eq.8) then
        call deriv_polargauss_damped(rvec_en,r_en)
      else
        stop 'deriv_orbitals: ibasis must be 4, 5, 6, 7 or 8'
      endif
      
      if ((ibasis.eq.4 .or. ibasis.eq.5 .or. ibasis.eq.8) .and. coef_is_diag) then
        ! Diagonal coefficient expansion.
        ! Every output is assigned once.
        do iorb=1,norb
          c_val = coef(iorb,iorb,iwf)
          do ie=nelec1,nelec2
            orb(ie,iorb) = c_val*phin(iorb,ie)
            do idim=1,ndim
              dorb(idim,ie,iorb) = c_val*dphin(idim,iorb,ie)
            enddo
            ddorb(ie,iorb) = c_val*d2phin(iorb,ie)
            
            d2orb_xx(ie,iorb) = c_val*d2phin_xx(iorb,ie)
            d2orb_yy(ie,iorb) = c_val*d2phin_yy(iorb,ie)
            d2orb_xy(ie,iorb) = c_val*d2phin_xy(iorb,ie)
          enddo
        enddo
        
      else
        ! Start the dense expansion from zero.
        orb = 0.0d0
        dorb = 0.0d0
        ddorb = 0.0d0
      
        d2orb_xx = 0.0d0
        d2orb_yy = 0.0d0
        d2orb_xy = 0.0d0
        
        do iorb=1,norb
          do m=1,nbasis
            c_val = coef(m,iorb,iwf)
            do ie=nelec1,nelec2
              orb(ie,iorb)=orb(ie,iorb)+c_val*phin(m,ie)
              do idim=1,ndim
                dorb(idim,ie,iorb)=dorb(idim,ie,iorb)+c_val*dphin(idim,m,ie)
              enddo
              ddorb(ie,iorb)=ddorb(ie,iorb)+c_val*d2phin(m,ie)
              
              d2orb_xx(ie,iorb)=d2orb_xx(ie,iorb)+c_val*d2phin_xx(m,ie)
              d2orb_yy(ie,iorb)=d2orb_yy(ie,iorb)+c_val*d2phin_yy(m,ie)
              d2orb_xy(ie,iorb)=d2orb_xy(ie,iorb)+c_val*d2phin_xy(m,ie)

            enddo
          enddo
        enddo
      endif

      return
      end

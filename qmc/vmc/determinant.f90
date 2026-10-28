      subroutine determinant(x,rvec_en,r_en,ddet_det,d2lndet,div_vd,determ)
! Written by Cyrus Umrigar starting from Kevin Schmidt's routine
      use all_tools_mod
      use constants_mod
      use control_mod
      use deriv_csf_mod
      use orb_mod
      use dorb_mod
      use coefs_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contrl_opt_mod
      use contr3_mod
      use optimo_mod
      use kinet_mod
      use d2orb_mod
      implicit real*8(a-h,o-z)

! Routine to calculate the value, gradient and Laplacian of the
! determinantal part of the wavefunction.
! Also derivatives wrt csf_coefs for optimizing them.

! Naming conventions are:
! "u" at end means up-spin
! "d" at end means dn-spin
! "d" at beginning of name indicates gradient
! "d2" at beginning of name indicates Laplacian
! "i" at end of portion of name indicates derivative wrt. csf_coef(i).
! "e" at end of portion of name indicates that it is not summed over all electrons
! "ln" indicates natural logarithm
! k=cartesian component, i=electron index
! slmui(i,idet)        = Inverse of transposed Slater matrix for idet stored as a vector
! fpu(k,i,idet)        = Gradient of elements of transposed Slater matrix for idet stored as a vector
! fppu(i,idet)         = Laplacian of elements of transposed Slater matrix for idet stored as a vector
! ddeti_deti(1,i,idet) = gradient(ln(det)) for electron i and determinant idet
! d2edeti_deti(i,idet) = Laplacian(det)/det for electron i and determinant idet
! determ               = sum of products of determinants times their coeffs.
! ddet_det(k,i)        = gradient(ln(determ))
! d2lndet              = Laplacian(ln(determ)) = Laplacian(determ)/determ - (gradient(determ)/determ)**2
! div_vd(i)            = Laplacian(ln(determ)) = Laplacian(determ)/determ - (gradient(determ)/determ)**2 for each electron
! ekinen(i)            = -0.5 Laplacian(determ)/determ

! Note that the first dimension of the slater matrices is MMAT_DIM = (MELEC/2)**2.
! The first dimension of the Slater matrices must be at least max(nup**2,ndn**2)
! So, we check in read_input that nup and ndn are each <= MELEC/2.

      common /dojasderiv/ ijasderiv
      common /s2l2/ el2_det, s2_det, s2o, el2o, s2sum, el2sum, s2cum, el2cum, s2cm2, el2cm2
      common /s2l2_switches/ icalc_L2, icalc_S2

      dimension x(3,*),rvec_en(3,nelec,*),r_en(nelec,*),ddet_det(3,*),div_vd(nelec)

! initialize the derivative arrays to zero
      do 10 i=1,nelec
        ekinen(i)=0
        ddet_det(1,i)=0
        ddet_det(2,i)=0
   10   ddet_det(3,i)=0
      d2lndet=0

! initialize the determinant arrays to one
      do 20 idet=1,ndetup
   20   detu(idet)=one
      do 22 idet=1,ndetdn
   22   detd(idet)=one
!     determ=0

! get orbitals and derivatives for all electrons
      if(iperiodic.eq.0 .or. iperiodic.eq.1) then

        if(inum_orb.eq.0) then
          if(nparmot.eq.0 .or. igradhess.eq.0) then
            call orbitals_loc_ana(0,rvec_en,r_en,orb,dorb,ddorb)
           else
            call deriv_orbitals(rvec_en,r_en,orb,dorb,ddorb)
          endif
         else
          call orbitals_loc_num(0,x,orb,dorb,ddorb)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pw(x,orb,dorb,ddorb)
!         write(6,'(''x(1,1),orb(1,1) from pw'',9f12.8)') x(1,1),orb(1,1),dorb(1,1,1),ddorb(1,1)
         else
          call orbitals_period_num(x,orb,dorb,ddorb)
!         write(6,'(''x(1,1),orb(1,1) from nu'',9f12.8)') x(1,1),orb(1,1),dorb(1,1,1),ddorb(1,1)
        endif

      endif

      call object_modified_by_index (orb_index)   !JT
      call object_modified_by_index (dorb_index)   !JT
      call object_modified_by_index (ddorb_index)   !JT

      if(ipr.ge.4) then
        do 26 iorb=1,norb
          write(6,'(''iorb,orb='',i3,(30f9.5))') iorb,(orb(i,iorb),i=1,nelec)
          write(6,'(''iorb,dorb1='',i3,(30f9.5))') iorb,(dorb(1,i,iorb),i=1,nelec)
          write(6,'(''iorb,dorb2='',i3,(30f9.5))') iorb,(dorb(2,i,iorb),i=1,nelec)
          write(6,'(''iorb,dorb3='',i3,(30f9.5))') iorb,(dorb(3,i,iorb),i=1,nelec)
   26     write(6,'(''iorb,ddorb='',i3,(30f9.5))') iorb,(ddorb(i,iorb),i=1,nelec)
      endif

! !fp
! Not only are we looking at the transpose of the Slater,
! but also, we fold the the matrix into a single column
! example: the following slater matrix
! ( phi_1(r1) phi_1(r2) )
! |                     |
! ( phi_2(r1) phi_2(r2) )
! is represented as
! ( phi_1(r1) )
! | phi_1(r2) |
! | phi_2(r1) |
! ( phi_2(r2) )
!
! The 3 nested loops are over
! 1) up electrons,
! 2) determinants
! 3) basis states setting up transpose of the Slater
! matrix in slmui to get inverse transpose.
! Also put derivatives in fpu and fppu.
! iworbdup(i,j) says which is the ith orbital of the jth up-determinant

! fpu example: (for idet=1)
! ( dphi_1(r1) / dx  dphi_2(r1) / dx  dphi_1(r2) / dx  dphi_2(r2) / dx )
! | dphi_1(r1) / dy  dphi_2(r1) / dy  dphi_1(r2) / dy  dphi_2(r2) / dy |
! ( dphi_1(r1) / dz  dphi_2(r1) / dz  dphi_1(r2) / dz  dphi_2(r2) / dz )

! fppu example: (for idet=1)
! ( d^2 phi_1(r1) / dx^2 + d^2 phi_1(r1) / dy^2 + d^2 phi_1(r1) / dz^2 )
! ( d^2 phi_2(r1) / dx^2 + d^2 phi_2(r1) / dy^2 + d^2 phi_2(r1) / dz^2 )
! ( d^2 phi_1(r2) / dx^2 + d^2 phi_1(r2) / dy^2 + d^2 phi_1(r2) / dz^2 )
! ( d^2 phi_2(r2) / dx^2 + d^2 phi_2(r2) / dy^2 + d^2 phi_2(r2) / dz^2 )
! !fp
      ik=-nup
      do 30 i=1,nup
        ik=ik+nup
        do 30 idet=1,ndetup
          jk=-nup
          do 30 j=1,nup
            jk=jk+nup
            slmui(i+jk,idet)=orb(i,iworbdup(j,idet))
            fpu(1,j+ik,idet)=dorb(1,i,iworbdup(j,idet))
            fpu(2,j+ik,idet)=dorb(2,i,iworbdup(j,idet))
            fpu(3,j+ik,idet)=dorb(3,i,iworbdup(j,idet))
   30       fppu(j+ik,idet)=ddorb(i,iworbdup(j,idet))

! loop through number of determinants calculating the inverse
! transpose matrices and their determinants
      do 40 idet=1,ndetup
!       write(6,'(''before matinv 40'',i2,9d14.6)') idet,slmui(1,idet),slmui(2,idet),slmui(3,idet),slmui(4,idet)
!         write(6,'(''before matinv 40'')')
!         do irow=0,(nup-1)
!            write(6,'(i3,(100f8.4))') irow,(slmui(irow*nup+k,idet),k=1,nup)
!         enddo
!         write(6,'(''x-posns of electrons: '',50f8.2)') (rvec_en(1,k,1), k=1,nelec)
!         write(6,'(''y-posns of electrons: '',50f8.2)') (rvec_en(2,k,1), k=1,nelec)
   40   call matinv(slmui(1,idet),nup,detu(idet))
!  40   write(6,'(''after matinv 40'',i2,9d14.6)') idet,slmui(1,idet),slmui(2,idet),slmui(3,idet),slmui(4,idet)

      call object_modified_by_index (slmui_index) !JT
      call object_modified_by_index (detu_index)  !JT

! repeat above for down spins
      ik=-ndn
      do 50 i=1,ndn
        ik=ik+ndn
        do 50 idet=1,ndetdn
          jk=-ndn
          do 50 j=1,ndn
            jk=jk+ndn
            slmdi(i+jk,idet)=orb(i+nup,iworbddn(j,idet))
            fpd(1,j+ik,idet)=dorb(1,i+nup,iworbddn(j,idet))
            fpd(2,j+ik,idet)=dorb(2,i+nup,iworbddn(j,idet))
            fpd(3,j+ik,idet)=dorb(3,i+nup,iworbddn(j,idet))
   50       fppd(j+ik,idet)=ddorb(i+nup,iworbddn(j,idet))

      do 60 idet=1,ndetdn
   60   call matinv(slmdi(1,idet),ndn,detd(idet))

      call object_modified_by_index (slmdi_index)  !JT
      call object_modified_by_index (detd_index)  !JT

      if(ipr.ge.4) write(6,'(''detu,detd'',9d12.5)') detu(1),detd(1)

! set up sum of slater determinants along with their
! coefficients that were included in input data
      do 70 idet=1,max(ndetup,ndetdn)
!       determ=determ+detu(idet)*detd(idet)*cdet(idet,iwf)

! zero out temporary derivative arrays
        do 70 i=1,nelec
          ddeti_deti(1,i,idet)=0
          ddeti_deti(2,i,idet)=0
          ddeti_deti(3,i,idet)=0
  70      d2edeti_deti(i,idet)=0

! loop through up spin electrons
! take inner product of transpose inverse with derivative
! vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2
      do 80 idet=1,ndetup
        ik=-nup
        do 80 i=1,nup
          ik=ik+nup
          do 80 j=1,nup
            ddeti_deti(1,i,idet)=ddeti_deti(1,i,idet)+slmui(j+ik,idet)*fpu(1,j+ik,idet)
            ddeti_deti(2,i,idet)=ddeti_deti(2,i,idet)+slmui(j+ik,idet)*fpu(2,j+ik,idet)
            ddeti_deti(3,i,idet)=ddeti_deti(3,i,idet)+slmui(j+ik,idet)*fpu(3,j+ik,idet)
   80       d2edeti_deti(i,idet)=d2edeti_deti(i,idet)+slmui(j+ik,idet)*fppu(j+ik,idet)

! repeat above for down spins
      do 90 idet=1,ndetdn
        ik=-ndn
        do 90 i=nup+1,nelec
          ik=ik+ndn
          do 90 j=1,ndn
            ddeti_deti(1,i,idet)=ddeti_deti(1,i,idet)+slmdi(j+ik,idet)*fpd(1,j+ik,idet)
            ddeti_deti(2,i,idet)=ddeti_deti(2,i,idet)+slmdi(j+ik,idet)*fpd(2,j+ik,idet)
            ddeti_deti(3,i,idet)=ddeti_deti(3,i,idet)+slmdi(j+ik,idet)*fpd(3,j+ik,idet)
   90       d2edeti_deti(i,idet)=d2edeti_deti(i,idet)+slmdi(j+ik,idet)*fppd(j+ik,idet)

! combine results for up and down spins to get d(det)/dx
! and d2(det)/dx in ddet_det and d2lndet respectively
      determ=0
      do 117 icsf=1,ncsf
        do 117 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
!         term=detu(idet)*detd(idet)*cdet(idet,iwf)
          if(ndn.ge.1) then
            term=detu(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            term=detu(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
          determ=determ+term
          do 115 i=1,nup
            iwdet=iwdetup(idet)
            ddet_det(1,i)=ddet_det(1,i)+ddeti_deti(1,i,iwdet)*term
            ddet_det(2,i)=ddet_det(2,i)+ddeti_deti(2,i,iwdet)*term
            ddet_det(3,i)=ddet_det(3,i)+ddeti_deti(3,i,iwdet)*term
            ekinen(i)=ekinen(i)+d2edeti_deti(i,iwdet)*term
  115       d2lndet=d2lndet+d2edeti_deti(i,iwdet)*term
          do 117 i=nup+1,nelec
            iwdet=iwdetdn(idet)
            ddet_det(1,i)=ddet_det(1,i)+ddeti_deti(1,i,iwdet)*term
            ddet_det(2,i)=ddet_det(2,i)+ddeti_deti(2,i,iwdet)*term
            ddet_det(3,i)=ddet_det(3,i)+ddeti_deti(3,i,iwdet)*term
            ekinen(i)=ekinen(i)+d2edeti_deti(i,iwdet)*term
  117       d2lndet=d2lndet+d2edeti_deti(i,iwdet)*term

! Inverse of full sum of determinants
      detinv=one/determ

! If detinv is infinite set it to largest number
!     if(isinf(detinv)) then ! should also work
      if(0.1*detinv==detinv) then
        write(6,'(''determinant: determ, detinv1='',9es12.4)') determ, detinv
        determ=tiny(determ)
        detinv=1/determ
        determ=0
        write(6,'(''determinant: determ, detinv2='',9es12.4)') determ, detinv
      endif

! multiply through to set up logarithmic first and second derivatives
      d2lndet=d2lndet*detinv
      do 120 i=1,nelec
        div_vd(i)=ekinen(i)*detinv
        ekinen(i)=-half*ekinen(i)*detinv
        do 120 k=1,ndim
          ddet_det(k,i)=ddet_det(k,i)*detinv
          div_vd(i)=div_vd(i)-ddet_det(k,i)**2
  120     d2lndet=d2lndet-ddet_det(k,i)**2

! EXACT L_Z^2 AND S^2 OBSERVABLE ESTIMATORS
      if (icalc_L2 .eq. 1) call calc_l2_exact(x, detinv, el2_det)
      if (icalc_S2 .eq. 1) call calc_s2_exact(detinv, s2_det)

! Derivatives wrt to csf_coefs for optimizing them
! Note that the arrays that are needed for vmc and dmc are over ndet but
! those that are needed for optimization only are over nparmcsf.
      if(index(mode,'fit').ne.0 .or. igradhess.gt.0 .or. l_opt_csf) then
        d2det_det=0
        do 125 i=1,nelec
  125     d2det_det=d2det_det-2*ekinen(i)
        do 150 iparm=1,nparmcsf
          icsf=iwcsf(iparm)
          d2deti_det(iparm)=0
          deti_det(iparm)=0
          do 130 i=1,nelec
            do 130 k=1,ndim
  130         ddeti_det(k,i,iparm)=0
          do 150 idet_in_csf=1,ndet_in_csf(icsf)
            idet=iwdet_in_csf(idet_in_csf,icsf)
            if(ndn.ge.1) then
              term=detu(iwdetup(idet))*detd(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)*detinv
             else
              term=detu(iwdetup(idet))*cdet_in_csf(idet_in_csf,icsf)*detinv
            endif
            deti_det(iparm)=deti_det(iparm)+term
            do 140 i=1,nup
              iwdet=iwdetup(idet)
              d2deti_det(iparm)=d2deti_det(iparm)+d2edeti_deti(i,iwdet)*term
              do 140 k=1,ndim
  140           ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+ddeti_deti(k,i,iwdet)*term
            do 150 i=nup+1,nelec
              iwdet=iwdetdn(idet)
              d2deti_det(iparm)=d2deti_det(iparm)+d2edeti_deti(i,iwdet)*term
              do 150 k=1,ndim
  150           ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+ddeti_deti(k,i,iwdet)*term
        if(ipr.ge.4) write(6,'(''deti_det(iparm) in determinant'',40d12.4)') (deti_det(iparm),iparm=1,nparmcsf)


! Derivatives with respect to orbital parameters (not orbital coefficients!).
        if(mod(iopt,10).eq.2) then
          do iparm=1,nparmcsf+nparmot
            do jparm=1,nparmcsf+nparmot
              detij_det(iparm,jparm)=0
            enddo
          enddo
        endif

! JT beg : deti_det for first csf (needed for unitary parametrization)
!          icsf=1
!          det1_det=0
!          do idet_in_csf=1,ndet_in_csf(icsf)
!            idet=iwdet_in_csf(idet_in_csf,icsf)
!            term=detu(idet)*detd(idet)*cdet_in_csf(idet_in_csf,icsf)*detinv
!            det1_det=det1_det+term
!          enddo
!          call object_modified_by_index (det1_det_index)
! JT end


        if(nparmot.gt.0) then
          call deriv_det_orb(orb,dorb,ddorb,detinv)
        endif


      endif

      return
      end
      
      
!--------------------------------------------------------------------------------------

subroutine deriv_det_orb(orb,dorb,ddorb,detinv)
  ! Exact derivatives using the independent-parameter map and spin trace contractions.
  ! Real, localized Gaussian basis-parameter derivatives, ibasis=4,5,8.
  ! Compact basis derivatives are evaluated by deriv_orbitals.
  ! common_allocations_deriv_det prepares the constraint map and persistent workspace.
  ! Parameter configuration changes require that preparation to be repeated.
  ! Orbital coefficients are fixed with respect to the basis parameters.
  ! Slater inverses and basis derivatives must already be valid and finite.
  ! Trace identities: Filippi, Assaraf and Moroni, JCP 144, 194105 (2016).
  use constants_mod
  use const_mod
  use dim_mod
  use contrl_opt_mod
  use contrl_opt2_mod
  use coefs_mod
  use wfsec_mod
  use optim_mod
  use optimo_mod
  use dorb_mod
  use dets_mod
  use slater_mod
  use deriv_phifun_mod
  use deriv_det_mod
  implicit none

  real(dp), intent(in) :: orb(nelec,norb), dorb(3,nelec,norb), ddorb(nelec,norb)
  real(dp), intent(in) :: detinv

  integer :: ip, p, q, np, ntot
  integer :: icsf, idc, idet, iu, idn, ie, k, jcsf, nterm, ga, gb
  real(dp) :: base, weight, value, lap_up, lap_dn
  logical :: newton

  if (nparmot == 0) return
  np = nparmot
  ntot = nparmcsf + np
  newton = mod(iopt,10) == 2
  deti_det(nparmcsf+1:ntot) = 0.0_dp
  ddeti_det(:,:,nparmcsf+1:ntot) = 0.0_dp
  d2deti_det(nparmcsf+1:ntot) = 0.0_dp
  ! Linear CSF coefficients have zero CSF-CSF second derivatives.
  if (newton) detij_det(1:ntot,1:ntot) = 0.0_dp
  if (igradhess == 0) return

  ! Collapse CSF weights before evaluating each unique spin determinant.
  orb_wup = 0.0_dp
  orb_wdn = 0.0_dp
  orb_oup = 0.0_dp
  orb_odn = 0.0_dp
  do icsf = 1, ncsf
    do idc = 1, ndet_in_csf(icsf)
      idet = iwdet_in_csf(idc,icsf)
      iu = iwdetup(idet)
      base = detu(iu)*cdet_in_csf(idc,icsf)*detinv
      if (ndn > 0) then
        idn = iwdetdn(idet)
        base = base*detd(idn)
      endif
      weight = csf_coef(icsf,iwf)*base
      orb_wup(iu) = orb_wup(iu) + weight
      if (ndn > 0) orb_wdn(idn) = orb_wdn(idn) + weight
    enddo
  enddo

  if (nup > 0) call spin_part(nup,0,ndetup,iworbdup,slmui,orb_wup,det_orb_scratch%use_rank_up,orb_oup)
  if (ndn > 0) call spin_part(ndn,nup,ndetdn,iworbddn,slmdi,orb_wdn,det_orb_scratch%use_rank_dn,orb_odn)

  ! Complete spin-product derivatives and CSF-orbital cross derivatives.
  ! The single-spin normalized derivatives have already been accumulated.
  do icsf = 1, ncsf
    do idc = 1, ndet_in_csf(icsf)
      idet = iwdet_in_csf(idc,icsf)
      iu = iwdetup(idet)
      base = detu(iu)*cdet_in_csf(idc,icsf)*detinv
      orb_osum = orb_oup(:,iu)
      if (ndn > 0) then
        idn = iwdetdn(idet)
        base = base*detd(idn)
        orb_osum = orb_osum + orb_odn(:,idn)
      endif
      weight = csf_coef(icsf,iwf)*base
      if (ndn > 0 .and. weight /= 0.0_dp) then
        ! Cross terms need only parameters with a nonzero spin first derivative.
        nterm = 0
        do p = 1, np
          if (orb_oup(p,iu) == 0.0_dp .and. orb_odn(p,idn) == 0.0_dp) cycle
          nterm = nterm + 1
          orb_term_params(nterm) = p
        enddo
        lap_up = sum(d2edeti_deti(1:nup,iu))
        lap_dn = sum(d2edeti_deti(nup+1:nelec,idn))
        do ga = 1, nterm
          p = orb_term_params(ga)
          ip = nparmcsf + p
          do ie = 1, nup
            do k = 1, ndim
              ddeti_det(k,ie,ip) = ddeti_det(k,ie,ip) &
                + weight*ddeti_deti(k,ie,iu)*orb_odn(p,idn)
            enddo
          enddo
          do ie = nup+1, nelec
            do k = 1, ndim
              ddeti_det(k,ie,ip) = ddeti_det(k,ie,ip) &
                + weight*ddeti_deti(k,ie,idn)*orb_oup(p,iu)
            enddo
          enddo
          d2deti_det(ip) = d2deti_det(ip) + weight &
            * (lap_up*orb_odn(p,idn) + lap_dn*orb_oup(p,iu))
          if (.not. newton) cycle
          do gb = 1, ga
            q = orb_term_params(gb)
            value = weight*(orb_oup(p,iu)*orb_odn(q,idn) + orb_odn(p,idn)*orb_oup(q,iu))
            call add_second(p,q,value)
          enddo
        enddo
      endif
      if (.not. newton) cycle
      do jcsf = 1, nparmcsf
        if (iwcsf(jcsf) /= icsf) cycle
        do p = 1, np
          if (orb_osum(p) == 0.0_dp) cycle
          ip = nparmcsf + p
          value = base*orb_osum(p)
          detij_det(ip,jcsf) = detij_det(ip,jcsf) + value
          detij_det(jcsf,ip) = detij_det(jcsf,ip) + value
        enddo
      enddo
    enddo
  enddo

contains

  subroutine add_second(pa,pb,x)
    integer, intent(in) :: pa, pb
    real(dp), intent(in) :: x
    integer :: ia, ja
    ia = nparmcsf + pa
    ja = nparmcsf + pb
    detij_det(ia,ja) = detij_det(ia,ja) + x
    if (ia /= ja) detij_det(ja,ia) = detij_det(ja,ia) + x
  end subroutine add_second

  subroutine spin_part(n,offset,nd,occ,invflat,weights,use_rank,olog)
    integer, intent(in) :: n, offset, nd, occ(:,:)
    real(dp), intent(in) :: invflat(:,:), weights(:)
    real(dp), intent(out) :: olog(:,:)
    logical, intent(in) :: use_rank(:)

    integer :: nb, nm, ng, d, b, t, a, j, e, dim, row
    integer :: mu, nu, pa, pb, ra, rb, ga, gb, ta, tb, l
    real(dp) :: w, tr, x, lap0, scale

    olog = 0.0_dp

    do d = 1, nd
      det_orb_scratch%s(1:n,1:n) = reshape(invflat(1:n*n,d),[n,n])
      nb = 0
      if (coef_is_diag) then
        ! Only occupied physical bases contribute. A slave can contribute
        ! even when its master basis is absent from this determinant.
        do j = 1, n
          b = occ(j,d)
          if (.not. any(orb_pmap(:,b) > 0)) cycle
          nb = nb + 1
          det_orb_scratch%bases(nb) = b
          det_orb_scratch%r(nb,1:n) = coef(b,b,iwf)*det_orb_scratch%s(j,1:n)
        enddo
      else
        do l = 1, orb_nactive
          b = orb_active_basis(l)
          if (.not. any(coef(b,occ(1:n,d),iwf) /= 0.0_dp)) cycle
          nb = nb + 1
          det_orb_scratch%bases(nb) = b
          do j = 1, n
            det_orb_scratch%r(nb,j) = coef(b,occ(j,d),iwf)
          enddo
        enddo
        if (nb > 0) det_orb_scratch%r(1:nb,1:n) = matmul(det_orb_scratch%r(1:nb,1:n),det_orb_scratch%s(1:n,1:n))
      endif
      if (nb == 0) cycle

      ! For physical parameter mu=(t,b): A_mu = u_mu*c_b^T,
      ! r_b = c_b^T*A^{-1}, and O_mu = r_b*u_mu.
      nm = 0
      ng = 0
      det_orb_scratch%where = 0
      do row = 1, nb
        b = det_orb_scratch%bases(row)
        do t = 1, notype
          a = orb_pmap(t,b)
          if (a == 0) cycle
          nm = nm + 1
          det_orb_scratch%mb(nm) = row
          det_orb_scratch%mt(nm) = t
          det_orb_scratch%mp(nm) = a
          det_orb_scratch%ms(nm) = orb_psign(t,b)
          det_orb_scratch%u(1:n,nm) = dparam(t,b,offset+1:offset+n)
          det_orb_scratch%omu(nm) = dot_product(det_orb_scratch%r(row,1:n),det_orb_scratch%u(1:n,nm))
          olog(a,d) = olog(a,d) + det_orb_scratch%ms(nm)*det_orb_scratch%omu(nm)
          if (det_orb_scratch%where(a) /= 0) cycle
          ng = ng + 1
          det_orb_scratch%groups(ng) = a
          det_orb_scratch%where(a) = ng
        enddo
      enddo

      w = weights(d)
      ! Olog is needed for spin cross terms even when this total weight cancels.
      if (w == 0.0_dp) cycle
      do ga = 1, ng
        a = det_orb_scratch%groups(ga)
        deti_det(nparmcsf+a) = deti_det(nparmcsf+a) + w*olog(a,d)
      enddo

      ! One-body trace derivative: tr(S*B_mu)-tr(S*B*S*A_mu).
      ! Compute B*S once per operator, then apply it to all physical u_mu.
      do dim = 1, ndim
        do j = 1, n
          det_orb_scratch%op(1:n,j) = dorb(dim,offset+1:offset+n,occ(j,d))
        enddo
        det_orb_scratch%op(1:n,1:n) = matmul(det_orb_scratch%op(1:n,1:n),det_orb_scratch%s(1:n,1:n))
        det_orb_scratch%work(1:n,1:nm) = matmul(det_orb_scratch%op(1:n,1:n),det_orb_scratch%u(1:n,1:nm))
        do mu = 1, nm
          a = nparmcsf + det_orb_scratch%mp(mu)
          row = det_orb_scratch%mb(mu)
          b = det_orb_scratch%bases(row)
          t = det_orb_scratch%mt(mu)
          scale = w*det_orb_scratch%ms(mu)
          do e = 1, n
            x = det_orb_scratch%r(row,e)*(ddparam(dim,t,b,offset+e)-det_orb_scratch%work(e,mu)) &
              + ddeti_deti(dim,offset+e,d)*det_orb_scratch%omu(mu)
            ddeti_det(dim,offset+e,a) = ddeti_det(dim,offset+e,a) + scale*x
          enddo
        enddo
      enddo
      do j = 1, n
        det_orb_scratch%op(1:n,j) = ddorb(offset+1:offset+n,occ(j,d))
      enddo
      det_orb_scratch%op(1:n,1:n) = matmul(det_orb_scratch%op(1:n,1:n),det_orb_scratch%s(1:n,1:n))
      det_orb_scratch%work(1:n,1:nm) = matmul(det_orb_scratch%op(1:n,1:n),det_orb_scratch%u(1:n,1:nm))
      lap0 = sum(d2edeti_deti(offset+1:offset+n,d))
      do mu = 1, nm
        row = det_orb_scratch%mb(mu)
        b = det_orb_scratch%bases(row)
        t = det_orb_scratch%mt(mu)
        a = nparmcsf + det_orb_scratch%mp(mu)
        x = dot_product(det_orb_scratch%r(row,1:n),d2dparam(t,b,offset+1:offset+n)) &
          - dot_product(det_orb_scratch%r(row,1:n),det_orb_scratch%work(1:n,mu)) + lap0*det_orb_scratch%omu(mu)
        d2deti_det(a) = d2deti_det(a) + w*det_orb_scratch%ms(mu)*x
      enddo
      if (.not. newton) cycle

      ! Use the exact contraction selected and sized in common_allocations.
      if (use_rank(d)) then
        det_orb_scratch%proj(1:nb,1:nm) = matmul(det_orb_scratch%r(1:nb,1:n),det_orb_scratch%u(1:n,1:nm))
        do mu = 1, nm
          pa = det_orb_scratch%mp(mu)
          ra = det_orb_scratch%mb(mu)
          do nu = 1, mu
            pb = det_orb_scratch%mp(nu)
            rb = det_orb_scratch%mb(nu)
            ! For one physical basis, the trace product cancels O_mu*O_nu
            ! exactly. Only its direct curvature (below) remains.
            if (ra == rb) cycle
            x = w*det_orb_scratch%ms(mu)*det_orb_scratch%ms(nu) &
              * (det_orb_scratch%omu(mu)*det_orb_scratch%omu(nu)-det_orb_scratch%proj(ra,nu)*det_orb_scratch%proj(rb,mu))
            ! Both ordered physical pairs contribute to a shared diagonal.
            if (pa == pb .and. mu /= nu) x = 2.0_dp*x
            call add_second(pa,pb,x)
          enddo
        enddo
      else
        det_orb_scratch%y(1:n,1:n,1:ng) = 0.0_dp
        do mu = 1, nm
          ga = det_orb_scratch%where(det_orb_scratch%mp(mu))
          do j = 1, n
            det_orb_scratch%y(1:n,j,ga) = det_orb_scratch%y(1:n,j,ga) + &
            det_orb_scratch%ms(mu)*det_orb_scratch%r(det_orb_scratch%mb(mu),j)*det_orb_scratch%u(1:n,mu)
          enddo
        enddo
        do ga = 1, ng
          do gb = 1, ga
            tr = 0.0_dp
            do j = 1, n
              tr = tr + dot_product(det_orb_scratch%y(1:n,j,ga),det_orb_scratch%y(j,1:n,gb))
            enddo
            call add_second(det_orb_scratch%groups(ga),det_orb_scratch%groups(gb),-w*tr)
          enddo
        enddo
        ! Group derivatives can have rank greater than one: retain the
        ! full logarithmic trace term and normalization product.
        do ga = 1, ng
          pa = det_orb_scratch%groups(ga)
          do gb = 1, ga
            pb = det_orb_scratch%groups(gb)
            call add_second(pa,pb,w*olog(pa,d)*olog(pb,d))
          enddo
        enddo
      endif

      ! Direct curvature: only two physical parameters on the SAME basis
      ! have a nonzero A_mu,nu. Map BOTH indices, independently of masters.
      do row = 1, nb
        b = det_orb_scratch%bases(row)
        do ta = 1, notype
          pa = orb_pmap(ta,b)
          if (pa == 0) cycle
          do tb = 1, ta
            pb = orb_pmap(tb,b)
            if (pb == 0) cycle
            x = w*orb_psign(ta,b)*orb_psign(tb,b) &
              * dot_product(det_orb_scratch%r(row,1:n),d2param(ta,tb,b,offset+1:offset+n))
            call add_second(pa,pb,x)
          enddo
        enddo
      enddo
    enddo
  end subroutine spin_part
end subroutine deriv_det_orb
      
!--------------------------------------------------------------------------------------

      subroutine calc_l2_exact(x, detinv, el2_det_out)
! Evaluates exact analytical L_z^2 direct local estimator
! GO, Oct 2026
      use all_tools_mod
      use constants_mod
      use control_mod
      use deriv_csf_mod
      use orb_mod
      use dorb_mod
      use d2orb_mod
      use coefs_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contrl_opt_mod
      use contr3_mod
      use optimo_mod
      use kinet_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*)
      real*8 :: detinv, el2_det_out
      
      real*8, allocatable, save :: lz_up(:), lz2_up(:), lz_dn(:), lz2_dn(:)
      real*8, allocatable, save :: T_x_u(:,:), T_y_u(:,:), M_u(:,:)
      real*8, allocatable, save :: T_x_d(:,:), T_y_d(:,:), M_d(:,:)
      logical, save :: is_l2_allocated = .false.
      
      real*8 :: H_xx, H_yy, H_xy, tr_M, tr_M2, sum_H, sum_D, tx, ty
      real*8 :: term_l2, term_csf
      
      if (.not. is_l2_allocated) then
        allocate(lz_up(ndetup), lz2_up(ndetup))
        if (ndn .gt. 0) allocate(lz_dn(ndetdn), lz2_dn(ndetdn))
        allocate(T_x_u(nup,nup), T_y_u(nup,nup), M_u(nup,nup))
        if (ndn .gt. 0) allocate(T_x_d(ndn,ndn), T_y_d(ndn,ndn), M_d(ndn,ndn))
        is_l2_allocated = .true.
      endif

      el2_det_out = 0.d0

      ! Calculate L_z and L_z^2 components for Up determinants
      do idet = 1, ndetup
        do i = 1, nup
          do k = 1, nup
            tx = 0.d0
            ty = 0.d0
            do m = 1, nup
               tx = tx + dorb(1, i, iworbdup(m, idet)) * slmui(m + (k-1)*nup, idet)
               ty = ty + dorb(2, i, iworbdup(m, idet)) * slmui(m + (k-1)*nup, idet)
            enddo
            T_x_u(i,k) = tx
            T_y_u(i,k) = ty
          enddo
        enddo

        tr_M = 0.d0
        tr_M2 = 0.d0
        sum_H = 0.d0
        sum_D = 0.d0

        do i = 1, nup
          do k = 1, nup
            M_u(i,k) = x(1,i) * T_y_u(i,k) - x(2,i) * T_x_u(i,k)
          enddo
        enddo

        do i = 1, nup
          tr_M = tr_M + M_u(i,i)
          do k = 1, nup
            tr_M2 = tr_M2 + M_u(i,k) * M_u(k,i)
          enddo
          
          H_xx = 0.d0
          H_yy = 0.d0
          H_xy = 0.d0
          do m = 1, nup
             H_xx = H_xx + d2orb_xx(i, iworbdup(m, idet)) * slmui(m + (i-1)*nup, idet)
             H_yy = H_yy + d2orb_yy(i, iworbdup(m, idet)) * slmui(m + (i-1)*nup, idet)
             H_xy = H_xy + d2orb_xy(i, iworbdup(m, idet)) * slmui(m + (i-1)*nup, idet)
          enddo
          
          sum_H = sum_H + x(1,i)**2 * H_yy + x(2,i)**2 * H_xx - 2.d0 * x(1,i) * x(2,i) * H_xy
          sum_D = sum_D + x(1,i) * T_x_u(i,i) + x(2,i) * T_y_u(i,i)
        enddo

        lz_up(idet) = tr_M
        ! Exact Eigenstate Lz2 estimator: -hbar^2 [ (Tr M)^2 - Tr(M^2) + sum_H - sum_D ]
        lz2_up(idet) = - ( tr_M**2 - tr_M2 + sum_H - sum_D )
      enddo

      ! Calculate L_z and L_z^2 components for Down determinants
      if (ndn .ge. 1) then
        do idet = 1, ndetdn
          do i = 1, ndn
            ie = i + nup
            do k = 1, ndn
              tx = 0.d0
              ty = 0.d0
              do m = 1, ndn
                 tx = tx + dorb(1, ie, iworbddn(m, idet)) * slmdi(m + (k-1)*ndn, idet)
                 ty = ty + dorb(2, ie, iworbddn(m, idet)) * slmdi(m + (k-1)*ndn, idet)
              enddo
              T_x_d(i,k) = tx
              T_y_d(i,k) = ty
            enddo
          enddo

          tr_M = 0.d0
          tr_M2 = 0.d0
          sum_H = 0.d0
          sum_D = 0.d0

          do i = 1, ndn
            ie = i + nup
            do k = 1, ndn
              M_d(i,k) = x(1,ie) * T_y_d(i,k) - x(2,ie) * T_x_d(i,k)
            enddo
          enddo

          do i = 1, ndn
            ie = i + nup
            tr_M = tr_M + M_d(i,i)
            do k = 1, ndn
              tr_M2 = tr_M2 + M_d(i,k) * M_d(k,i)
            enddo
            
            H_xx = 0.d0
            H_yy = 0.d0
            H_xy = 0.d0
            do m = 1, ndn
               H_xx = H_xx + d2orb_xx(ie, iworbddn(m, idet)) * slmdi(m + (i-1)*ndn, idet)
               H_yy = H_yy + d2orb_yy(ie, iworbddn(m, idet)) * slmdi(m + (i-1)*ndn, idet)
               H_xy = H_xy + d2orb_xy(ie, iworbddn(m, idet)) * slmdi(m + (i-1)*ndn, idet)
            enddo
            
            sum_H = sum_H + x(1,ie)**2 * H_yy + x(2,ie)**2 * H_xx - 2.d0 * x(1,ie) * x(2,ie) * H_xy
            sum_D = sum_D + x(1,ie) * T_x_d(i,i) + x(2,ie) * T_y_d(i,i)
          enddo

          lz_dn(idet) = tr_M
          lz2_dn(idet) = - ( tr_M**2 - tr_M2 + sum_H - sum_D )
        enddo
      endif

      ! Combine components across all CSFs to form exact total L_z^2 expectation value
      do icsf = 1, ncsf
        do idet_in_csf = 1, ndet_in_csf(icsf)
          idet = iwdet_in_csf(idet_in_csf, icsf)
          if (ndn .ge. 1) then
            term_csf = detu(iwdetup(idet)) * detd(iwdetdn(idet)) * csf_coef(icsf,iwf) * cdet_in_csf(idet_in_csf,icsf)
            ! Total operator cross-terms: L_z^2 = L_{z,up}^2 + L_{z,dn}^2 + 2 L_{z,up} L_{z,dn}
            ! The -2.d0 corrects for the imaginary phase in the cross product since L_z_up is Tr(M_up).
            term_l2 = lz2_up(iwdetup(idet)) + lz2_dn(iwdetdn(idet)) - 2.d0 * lz_up(iwdetup(idet)) * lz_dn(iwdetdn(idet))
          else
            term_csf = detu(iwdetup(idet)) * csf_coef(icsf,iwf) * cdet_in_csf(idet_in_csf,icsf)
            term_l2 = lz2_up(iwdetup(idet))
          endif
          el2_det_out = el2_det_out + term_l2 * term_csf * detinv
        enddo
      enddo

      return
      end

!--------------------------------------------------------------------------------------

      subroutine calc_s2_exact(detinv, s2_det_out)
! Calculates S^2 using Sherman-Morrison rank-1 swap updates
! GO, Oct 2026
      use all_tools_mod
      use constants_mod
      use control_mod
      use deriv_csf_mod
      use orb_mod
      use dorb_mod
      use d2orb_mod
      use coefs_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contrl_opt_mod
      use contr3_mod
      use optimo_mod
      use kinet_mod
      implicit real*8(a-h,o-z)

      real*8 :: detinv, s2_det_out
      real*8 :: ratio_sum, rat_up, rat_dn, term

      s2_det_out = 0.d0
      if (ndn .ge. 1 .and. nup .ge. 1) then
        do i = 1, nup
          do j = nup+1, nelec
            ratio_sum = 0.d0
            do icsf = 1, ncsf
              do idet_in_csf = 1, ndet_in_csf(icsf)
                idet = iwdet_in_csf(idet_in_csf, icsf)
                rat_up = 0.d0
                ik = nup*(i-1)
                do k = 1, nup
                  rat_up = rat_up + orb(j, iworbdup(k, iwdetup(idet))) * slmui(k+ik, iwdetup(idet))
                enddo
                rat_dn = 0.d0
                ik = ndn*(j-nup-1)
                do m = 1, ndn
                  rat_dn = rat_dn + orb(i, iworbddn(m, iwdetdn(idet))) * slmdi(m+ik, iwdetdn(idet))
                enddo
                if (ndn .ge. 1) then
                  term = detu(iwdetup(idet)) * detd(iwdetdn(idet)) * csf_coef(icsf,iwf) * cdet_in_csf(idet_in_csf,icsf)
                else
                  term = detu(iwdetup(idet)) * csf_coef(icsf,iwf) * cdet_in_csf(idet_in_csf,icsf)
                endif
                ratio_sum = ratio_sum + term * rat_up * rat_dn
              enddo
            enddo
            s2_det_out = s2_det_out + ratio_sum * detinv
          enddo
        enddo
      endif

      return
      end
      


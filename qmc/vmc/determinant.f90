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

      dimension x(3,*),rvec_en(3,nelec,*),r_en(nelec,*),ddet_det(3,*),div_vd(nelec)
      dimension dporb(notype,nbasis,nelec,norb),d2porb(notype,notype,nbasis,nelec,norb)
      dimension ddporb(3,notype,nbasis,nelec,norb),d2dporb(notype,nbasis,nelec,norb)

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
            call deriv_orbitals(rvec_en,r_en,orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb)
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
        if(iopt.eq.2) then
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
          if(iconstrain_gauss_orbs.eq.1) then
            call constrained_deriv_det_orb(dporb,d2porb,ddporb,d2dporb)
          endif
          call deriv_det_orb(orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb,detinv)
        endif


      endif

      return
      end
      
      
!--------------------------------------------------------------------------------------
      subroutine deriv_det_orb(orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb,detinv)
! Written by A.D.Guclu, Apr 2006
! Heavily modified by Gokhan Oztarhan, Oct 2023

! Calculates derivatives wrt orbital parameters for optimization

! explanations on local variables:

! detui(iparm,idet)  = derivative wrt the parameter iparm in up determ. idet
! detuij(jparm,idet) = derivative wrt the parameter jparm, of the current detui, in up determ. idet
!                       in other words, second derivatives wrt iparm&jparm, but iparm is not stored.
! ddetui(1:3,ie,idet)= velocity of the current detui, iparm is not stored.
! d2detui(idet)      = laplacian of the current detui, iparm is not stored.

! slmfui(nup,nup) = Inverse of up Slater matrix (not transposed), stored in 2d array


! NOTES AFTER COMPLETE REVISION, OCT 2023

! Optimization of basis function parameters is enabled for the case 
! in which orbitals are a linear combination of basis functions. 

! Derivatives with respect to basis functions need traces of the form
! tr(A^-1 * dA * A^-1 * B)
! where A is the Slater matrix, dA and B are the matrices containing derivatives.

! Inverse of transposed Slater matrix, is stored in slmui and slmdi arrays.
! Trace formulas need the inverse of Slater matrix (not transposed).
! Transpose of slmui and slmdi are stored in slmfui and slmfdi 2d arrays
! in order to obtain inverse of Slater matrix.

! Note that the form of the Slater matrix is the following.
! ( phi_1(r1) phi_2(r1) )
! ( phi_1(r2) phi_2(r2) )

! Useful properties in order to reduce the number of operations:

! If A and B are NxN matrices, and B has only one non-zero column or row,
! tr(A*B)^2 = tr((A*B)^2).

! For the traces of the form tr(A*B) and tr(A*B*A*B), use cyclic property
! so that the resulting matrices of which trace is calculated have a single
! non-zero column or row. Therefore, only one element of the matrix products
! has to be calculated in order to obtain the trace.
! If B has only one non-zero column, calculate tr(A*B) and tr(A*B*A*B).
! If B has only one non-zero row, calculate tr(B*A) and tr(B*A*B*A).


      use dorb_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use wfsec_mod
      use contrl_opt_mod
      use optimo_mod
      use coefs_mod
      implicit real*8(a-h,o-z)

! arguments:
      dimension orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb)
      dimension dporb(notype,nbasis,nelec,norb),d2porb(notype,notype,nbasis,nelec,norb)
      dimension ddporb(3,notype,nbasis,nelec,norb),d2dporb(notype,nbasis,nelec,norb)

! local arrays:
      dimension detui(nparmd,ndetup),detdi(nparmd,ndetdn)
      dimension detuij(nparmd,ndetup),detdij(nparmd,ndetdn)
      dimension ddetui(3,nelec,ndetup),ddetdi(3,nelec,ndetdn)
      dimension d2detui(ndetup),d2detdi(ndetdn)
      dimension slmfui(nup,nup,ndetup),slmfdi(ndn,ndn,ndetdn)
      dimension amatu(nup,nup),amatd(ndn,ndn)
      dimension amattempu(nup,nup),amattempd(ndn,ndn)
      

! initializations of local and output variables:
! complete initialization of detij_det is done in determinant()
      do iparm0=1,nparmot
        iparm=iparm0+nparmcsf
        d2deti_det(iparm)=0
        deti_det(iparm)=0
        do idet=1,ndetup
          detui(iparm0,idet)=0
        enddo
        do idet=1,ndetdn
          detdi(iparm0,idet)=0
        enddo
        do i=1,nelec
          do k=1,ndim
            ddeti_det(k,i,iparm)=0
          enddo
        enddo
      enddo
      
      ! Save inverse of Slater matrix (not transposed) to 2d array 
      ! by calculating the transpose of transposed inverse of Slater matrix (slmui)
      slmfui = 0
      do idet = 1, ndetup
        do i = 1, nup
          jk = -nup
          do j = 1, nup
            jk = jk + nup
            slmfui(i,j,idet) = slmui(i+jk,idet)
          enddo
        enddo
      enddo
      slmfdi = 0
      do idet = 1, ndetdn
        do i = 1, ndn
          jk = -ndn
          do j = 1, ndn
            jk = jk + ndn
            slmfdi(i,j,idet) = slmdi(i+jk,idet)
          enddo
        enddo
      enddo
      
      iparm0 = 0
      do it = 1, notype
        do ip = 1, iabs(nparmo(it))
          iparm0 = iparm0 + 1
          
          ! Following are array assignments
          ddetui = 0
          d2detui = 0
          detuij = 0
          ddetdi = 0
          d2detdi = 0
          detdij = 0
        
          ! Up electrons
          do idet = 1, ndetup
          
            ! Constract A^-1 * dA * A^-1
            amattempu = 0
            do k = 1, nup
              do j = 1, nup
                do i = 1, nup
                  amattempu(i,j) = amattempu(i,j) + slmfui(i,j,idet) * dporb(it,ip,j,iworbdup(k,idet))
                enddo
              enddo
            enddo
            amatu = 0
            do k = 1, nup
              do j = 1, nup
                do i = 1, nup
                  amatu(i,j) = amatu(i,j) + amattempu(i,j) * slmfui(j,k,idet)
                enddo
              enddo
            enddo
            
            ! First parameter derivatives
            do i = 1, nup
              do j = 1, nup
                detui(iparm0,idet) = detui(iparm0,idet) &
     &             + slmfui(j,i,idet) * dporb(it,ip,i,iworbdup(j,idet))
              enddo
            enddo

            ! Velocity of first parameter derivatives
            do ie = 1, nup
              do j = 1, nup
                do idim=1,ndim
                  ddetui(idim,ie,idet) = ddetui(idim,ie,idet) &
     &               + ddporb(idim,it,ip,ie,iworbdup(j,idet)) * slmfui(j,ie,idet) &
     &               - dorb(idim,ie,iworbdup(j,idet)) * amatu(j,ie)
                enddo
              enddo
            enddo
            
            ! Laplacian of first parameter derivatives
            do ie = 1, nup
              do j = 1, nup
                d2detui(idet) = d2detui(idet) &
     &             + d2dporb(it,ip,ie,iworbdup(j,idet)) * slmfui(j,ie,idet) &
     &             - ddorb(ie,iworbdup(j,idet)) * amatu(j,ie)
              enddo
            enddo
            
            ! Second parameter derivatives if doing newton
            if (iopt .eq. 2) then
              jparm0 = 0
              do jt = 1, notype
                do jp = 1, iabs(nparmo(jt))
                  jparm0 = jparm0 + 1
                  if (jparm0 .le. iparm0) then ! symmetry of Hessian
                    if (jp .eq. ip) then
                      do i = 1, nup
                        do j = 1, nup
                          detuij(jparm0,idet) = detuij(jparm0,idet) &
     &                       + slmfui(j,i,idet) * d2porb(it,jt,jp,i,iworbdup(j,idet)) &
     &                       - amatu(j,i) * dporb(jt,jp,i,iworbdup(j,idet))
                        enddo
                      enddo
                    else
                      do i = 1, nup
                        do j = 1, nup
                          detuij(jparm0,idet) = detuij(jparm0,idet) &
     &                       - amatu(j,i) * dporb(jt,jp,i,iworbdup(j,idet))
                        enddo
                      enddo
                    endif ! jp .eq. ip
                  endif ! jparm0 .le. iparm0
                enddo ! jt
              enddo ! jp
            endif ! iopt

          enddo ! idet (up electrons)
          
          
          ! Down electrons
          do idet = 1, ndetdn
          
            ! Constract A^-1 * dA * A^-1
            amattempd = 0
            do k = 1, ndn
              do j = 1, ndn
                do i = 1, ndn
                  amattempd(i,j) = amattempd(i,j) + slmfdi(i,j,idet) * dporb(it,ip,j+nup,iworbddn(k,idet))
                enddo
              enddo
            enddo
            amatd = 0
            do k = 1, ndn
              do j = 1, ndn
                do i = 1, ndn
                  amatd(i,j) = amatd(i,j) + amattempd(i,j) * slmfdi(j,k,idet)
                enddo
              enddo
            enddo
            
            ! First parameter derivatives
            do i = 1, ndn
              do j = 1, ndn
                detdi(iparm0,idet) = detdi(iparm0,idet) &
     &             + slmfdi(j,i,idet) * dporb(it,ip,i+nup,iworbddn(j,idet))
              enddo
            enddo

            ! Velocity of first parameter derivatives
            do ie = 1, ndn
              do j = 1, ndn
                do idim=1,ndim
                  ddetdi(idim,ie,idet) = ddetdi(idim,ie,idet) &
     &               + ddporb(idim,it,ip,ie+nup,iworbddn(j,idet)) * slmfdi(j,ie,idet) &
     &               - dorb(idim,ie+nup,iworbddn(j,idet)) * amatd(j,ie)
                enddo
              enddo
            enddo
            
            ! Laplacian of first parameter derivatives
            do ie = 1, ndn
              do j = 1, ndn
                d2detdi(idet) = d2detdi(idet) &
     &             + d2dporb(it,ip,ie+nup,iworbddn(j,idet)) * slmfdi(j,ie,idet) &
     &             - ddorb(ie+nup,iworbddn(j,idet)) * amatd(j,ie)
              enddo
            enddo
            
            ! Second parameter derivatives if doing newton
            if (iopt .eq. 2) then
              jparm0 = 0
              do jt = 1, notype
                do jp = 1, iabs(nparmo(jt))
                  jparm0 = jparm0 + 1
                  if (jparm0 .le. iparm0) then ! symmetry of Hessian
                    if (jp .eq. ip) then
                      do i = 1, ndn
                        do j = 1, ndn
                          detdij(jparm0,idet) = detdij(jparm0,idet) &
     &                       + slmfdi(i,j,idet) * d2porb(it,jt,jp,i+nup,iworbddn(j,idet)) &
     &                       - amatd(j,i) * dporb(jt,jp,i+nup,iworbddn(j,idet))
                        enddo
                      enddo
                    else
                      do i = 1, ndn
                        do j = 1, ndn
                          detdij(jparm0,idet) = detdij(jparm0,idet) &
     &                       - amatd(j,i) * dporb(jt,jp,i+nup,iworbddn(j,idet))
                        enddo
                      enddo
                    endif ! jp .eq. ip
                  endif ! jparm0 .le. iparm0
                enddo ! jt
              enddo ! jp
            endif ! iopt
          
          enddo ! idet (down electrons)
          
          
! now put determinants together to get the final results:
          iparm=iparm0+nparmcsf
          do icsf=1,ncsf
            do idet_in_csf=1,ndet_in_csf(icsf)

              idet=iwdet_in_csf(idet_in_csf,icsf)

              if(ndn.ge.1) then
                term=detu(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
              else
                term=detu(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
              endif
              term=term*detinv
              
              deti_det(iparm)=deti_det(iparm)+ &
     &                 (detui(iparm0,iwdetup(idet))+detdi(iparm0,iwdetdn(idet)))*term
     
              d2deti_det(iparm) = d2deti_det(iparm) + term * (d2detui(iwdetup(idet)) + d2detdi(iwdetdn(idet)))

              do i = 1, nup
                iwdet = iwdetup(idet)
                do k = 1, ndim
                  ddeti_det(k,i,iparm) = ddeti_det(k,i,iparm) + term &
     &              * (ddetui(k,i,iwdet) + ddeti_deti(k,i,iwdet) &
     &              * (detui(iparm0,iwdetup(idet)) + detdi(iparm0,iwdetdn(idet)))) 
                enddo
                d2deti_det(iparm) = d2deti_det(iparm) &
     &            + term * (d2edeti_deti(i,iwdet) &
     &            * (detui(iparm0,iwdetup(idet)) + detdi(iparm0,iwdetdn(idet))))
              enddo
              do i = nup + 1, nelec
                iwdet = iwdetdn(idet)
                do k = 1, ndim
                  ddeti_det(k,i,iparm) = ddeti_det(k,i,iparm) + term &
     &              * (ddetdi(k,i-nup,iwdet) + ddeti_deti(k,i,iwdet) &
     &              * (detui(iparm0,iwdetup(idet)) + detdi(iparm0,iwdetdn(idet))))
                enddo
                d2deti_det(iparm) = d2deti_det(iparm) &
     &            + term * (d2edeti_deti(i,iwdet) &
     &            * (detui(iparm0,iwdetup(idet)) + detdi(iparm0,iwdetdn(idet))))
              enddo

              if (iopt .eq. 2) then
                do jparm0=1,iparm0
                  jparm=jparm0+nparmcsf
                  detij_det(iparm,jparm) = detij_det(iparm,jparm) + term &
     &              * (detuij(jparm0,iwdetup(idet)) + detdij(jparm0,iwdetdn(idet)) &
     &              + (detui(iparm0,iwdetup(idet)) + detdi(iparm0,iwdetdn(idet))) &
     &              * (detui(jparm0,iwdetup(idet)) + detdi(jparm0,iwdetdn(idet))))
                  detij_det(jparm,iparm)=detij_det(iparm,jparm)
                enddo
              endif
              
            enddo ! idet_in_csf
          enddo ! icsf
        
        enddo ! ip
      enddo ! it


      return
      end
      

!--------------------------------------------------------------------------------------

      subroutine constrained_deriv_det_orb(dporb,d2porb,ddporb,d2dporb)
! Created: Gokhan Oztarhan, Oct 2023

! Calculates the summation of basis function derivatives 
! if they are constrained to be equal (or equal with minus sign).
! Simply sum basis function derivatives together into the orbital derivatives 
! array before calculating the derivatives of the Slater determinant.
! Therefore, the value of the parameter, which will be optimized, 
! in dporb, d2porb, ddporb and d2dporb arrays also contain the basis function
! derivatives of the constrained parameters.

! This subroutine is created for revised deriv_det_orb subroutine since
! old constrained_deriv_det_orb (created by A. Mehta) is not compatible with
! the revised deriv_det_orb. However, the usage of norb_constraints and 
! orb_constraints arrays are not changed.
      
      use dorb_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use coefs_mod
      use wfsec_mod
      use optimo_mod
      use basic_tools_mod
      use objects_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)
      
      ! arguments
      dimension dporb(notype,nbasis,nelec,norb),d2porb(notype,notype,nbasis,nelec,norb)
      dimension ddporb(3,notype,nbasis,nelec,norb),d2dporb(notype,nbasis,nelec,norb)
      
      do it = 1, notype
        do ip = 1, norb_constraints(it)
        
          consgn = real(sign(1,orb_constraints(it,ip,2)))
          
          ipcon2 = iabs(orb_constraints(it,ip,2))
        
          do iorb = 1, norb
            do ie = 1, nelec
              dporb(it,ipcon2,ie,iorb) = &
     &          dporb(it,ipcon2,ie,iorb) &
     &          + consgn * dporb(it,orb_constraints(it,ip,1),ie,iorb)
     
              do idim = 1, ndim
                ddporb(idim,it,ipcon2,ie,iorb) = &
     &            ddporb(idim,it,ipcon2,ie,iorb) &
     &            + consgn * ddporb(idim,it,orb_constraints(it,ip,1),ie,iorb)
              enddo
              
              d2dporb(it,ipcon2,ie,iorb) = &
     &          d2dporb(it,ipcon2,ie,iorb) &
     &          + consgn * d2dporb(it,orb_constraints(it,ip,1),ie,iorb)
     
              do jt = 1, notype
                d2porb(it,jt,ipcon2,ie,iorb) = &
     &            d2porb(it,jt,ipcon2,ie,iorb) &
     &            + consgn * d2porb(it,jt,orb_constraints(it,ip,1),ie,iorb)
              enddo
            enddo ! ie
          enddo ! iorb
          
        enddo ! ip
      enddo ! it
      
      return
      end
      

!--------------------------------------------------------------------------------------

      subroutine constrained_deriv_det_orb_deprecated(orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb,detinv)
!  Written by Abhijit Mehta, May 2010
!    Extensively modified January 2011
!  Calculates derivatives wrt orbital parameters for optimization if some of
!    the parameters are constrained to be equal (eg, we may want all gaussians
!      to have the same width, so we can optimize just 1 parameter rather than
!       nelec individual parameters.)
!  This routine calls deriv_det_orb to calculate derivatives as if each
!     individual orbital parameter were independent, then just sums them up
!     (using the chain rule) to get the derivative with respect to the single
!       parameter.

      use dorb_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use coefs_mod
      use wfsec_mod
      use optimo_mod
      use basic_tools_mod
      use objects_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)
! arguments
      dimension orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb)
      dimension dporb(notype,nbasis,nelec,norb),d2porb(notype,notype,nbasis,nelec,norb)
      dimension ddporb(3,notype,nbasis,nelec,norb),d2dporb(notype,nbasis,nelec,norb)
! temporary variables, used to save values
      dimension nparmo_temp(notype)
      dimension deti_det_temp(nparmd), ddeti_det_temp(3,nelec,nparmd)
      dimension d2deti_det_temp(nparmd), detij_det_temp(nparmd, nparmd)
!      dimension detij_partial_temp(nparmd,nparmcsf+4*nbasis)
      dimension iwo_temp(nbasis,notype) !GO
! ACM: For optimized constraints, we shall run deriv_det_orb as if we were
!  varying each orbital parameter separately, so we must change the values of
!  nparmo(it), nparmot, nparmd, and iwo(ip,it).  We calculate overall derivatives using
!   chain rule, and restore values after calling deriv_det_orb
!  This is a quick and dirty solution, but it should work.
! Jan 2011: for now, calculate derivatives wrt all orbital params when there's a constraint
!    In the future, we should fix this to only compute the necessary parameters.
      nparmo_temp = nparmo  ! this is an array operation nparmo_temp(:) = nparmo(:)
      nparmot_temp = nparmot
      nparmd_temp = nparmd
      iwo_temp = iwo
      do it=1,notype
        if (nparmo(it).lt.0) then
          nparmot = nparmot - iabs(nparmo(it)) + nbasis
          nparmo(it) = nbasis
          do ib = 1,nbasis  ! we'll calculate derivatives wrt all orbital params
            iwo(ib, it) = ib
          enddo
        endif
      enddo
      nparmd = nparmcsf+nparmot

!  Housekeeping to make sure arrays are the right size - hope I'm doing this right
      call object_modified('nparmd')
      call object_modified('iwo')
      call alloc('deti_det', deti_det, nparmd)
      call alloc('ddeti_det', ddeti_det, 3,nelec, nparmd)
      call alloc('d2deti_det', d2deti_det, nparmd)
      call alloc('detij_det', detij_det, nparmd, nparmd)

! Derivatives with respect to orbital parameters (not orbital coefficients!).
!   This is a cut and paste from determinant()
      if(iopt.eq.2) then
        do iparm=1,nparmcsf+nparmot
          do jparm=1,nparmcsf+nparmot
            detij_det(iparm,jparm)=0
          enddo
        enddo
      endif

      call deriv_det_orb(orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb,detinv)

!      write(6,'(''After deriv_det_orb, detij_det ='')')  ! ACM debug
!      do i=1,nparmd
!        write(6,'(20g12.4)') (detij_det(i,j), j=1,nparmd)
!      enddo

!  Do chain rule to calculate derivative with respect to constrained parameter

        iparm0 = nparmcsf ! which parameter we're on for final result
        iparm1 = nparmcsf ! which parameter we're on from deriv_det_orb
        do it = 1,notype
          if (nparmo_temp(it).lt.0) then
!           Do chain rule by summing up derivatives:
            do icon=1,norb_constraints(it)
              iparm_sum_index = iparm1+iabs(orb_constraints(it,icon,2)) ! where to keep sum
              consgn = real(sign(1,orb_constraints(it,icon,2))) !whether this constraint is a mirror constraint
              iparm_summand_index = iparm1+orb_constraints(it,icon,1) ! constrained param
              deti_det(iparm_sum_index) = deti_det(iparm_sum_index) + consgn*deti_det(iparm_summand_index)
              ddeti_det(:,:,iparm_sum_index) = ddeti_det(:,:,iparm_sum_index) + consgn*ddeti_det(:,:,iparm_summand_index)
              d2deti_det(iparm_sum_index) = d2deti_det(iparm_sum_index) + consgn*d2deti_det(iparm_summand_index)
              do icon2=1,norb_constraints(it) ! sum to get partials -detij_det
                iparm_sum_index2 = iparm1+iabs(orb_constraints(it,icon2,2))
                consgn2 = real(sign(1,orb_constraints(it,icon2,2)))
                iparm_summand_index2 = iparm1 + orb_constraints(it,icon2,1)
                detij_det(iparm_sum_index,iparm_sum_index2) = detij_det(iparm_sum_index,iparm_sum_index2) &
     &             + consgn*consgn2*detij_det(iparm_summand_index,iparm_summand_index2)
              enddo
            enddo
!           Now put the summed values in output variables (currently labeled _temp), but we switch them later
            do ip=1,iabs(nparmo_temp(it))  ! set sums to output variables
              iparm_final_index = iparm0+ip
              iparm_sum_index = iparm1+iwo_temp(ip,it)
              deti_det_temp(iparm_final_index) = deti_det(iparm_sum_index)
              ddeti_det_temp(:,:,iparm_final_index) = ddeti_det(:,:,iparm_sum_index)
              d2deti_det_temp(iparm_final_index) = d2deti_det(iparm_sum_index)
              do ip2=1,iabs(nparmo_temp(it)) ! second sum for detij_det
                iparm_final_index2 = iparm0+ip2
                iparm_sum_index2 = iparm1+iwo_temp(ip2,it)
                detij_det_temp(iparm_final_index,iparm_final_index2) = detij_det(iparm_sum_index,iparm_sum_index2)
!                write(6,'(6i5)') ip, iwo_temp(ip,it), iparm_final_index, iparm_final_index2, iparm_sum_index, iparm_sum_index2
!                write(6,'(3g12.4)')  iparm_final_index, iparm_final_index2, detij_det_temp(iparm_final_index,iparm_final_index2)
!                write(6,'(3g12.4)') iparm_sum_index,iparm_sum_index2, detij_det(iparm_sum_index,iparm_sum_index2)
              enddo
            enddo
!    This code was from when nparmo(it)=-1 just meant constrain ALL params
!            iparm0 = iparm0 + 1
!            deti_det_temp(iparm0) = sum(deti_det(iparm1+1:iparm1+nbasis))
!            ddeti_det_temp(:,:,iparm0) = sum(ddeti_det(:,:,iparm1+1:iparm1+nbasis), DIM=3)
!            d2deti_det_temp(iparm0) = sum(d2deti_det(iparm1+1:iparm1+nbasis))
!            detij_partial_temp(iparm0,1:nparmd) = sum(detij_det(iparm1+1:iparm1+nbasis,:), DIM=1)
            iparm0 = iparm0 + iabs(nparmo_temp(it))
            iparm1 = iparm1 + nbasis
          else
            iparm0 = iparm0 + nparmo(it)
            iparm1 = iparm1 + nparmo(it)
          endif
        enddo
!      write(6,'(''After applying constraints, detij_det ='')') ! ACM debug
!      do i=1,nparmd
!        write(6,'(20g12.4)') (detij_det(i,j), j=1,nparmd)
!      enddo
!
!      write(6,'(''After applying constraints, detij_det_temp ='')')
!      do i=1,nparmd_temp
!        write(6,'(20g12.4)') (detij_det_temp(i,j), j=1,nparmd_temp)
!      enddo

!  Finish up doing sum over second index to get detij_det(iparm,jparm) - NO LONGER NEEDED
!        jparm0 = nparmcsf
!        jparm1 = nparmcsf
!        do it = 1,notype
!         if (nparmo_temp(it).eq.-1) then
!            jparm0 = jparm0 + 1
!            detij_det_temp(:,jparm0) = sum(detij_partial_temp(:,jparm1+1:jparm1+nbasis), DIM=2)
!            jparm1 = jparm1 + nbasis
!          else
!            jparm0 = jparm0 + nparmo(it)
!            jparm1 = jparm1 + nparmo(it)
!         endif
!        enddo

!   Restore variables to their original values, update sizes
!     of derivative arrays and give them correct values
        nparmo = nparmo_temp
        nparmot = nparmot_temp
        nparmd = nparmd_temp
        iwo = iwo_temp
        call object_modified('nparmd')
        call object_modified('iwo')
        call alloc('deti_det', deti_det, nparmd)
        call alloc('ddeti_det', ddeti_det, 3,nelec, nparmd)
        call alloc('d2deti_det', d2deti_det, nparmd)
        call alloc('detij_det', detij_det, nparmd, nparmd)
        deti_det = deti_det_temp
        ddeti_det = ddeti_det_temp
        d2deti_det = d2deti_det_temp
        detij_det = detij_det_temp


!      write(6,'(''After constrained_deriv_det_orb, detij_det ='')') ! ACM debug
!      do i=1,nparmd
!        write(6,'(20g12.4)') (detij_det(i,j), j=1,nparmd)
!      enddo


        return
        end

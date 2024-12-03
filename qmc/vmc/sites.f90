      subroutine sites(x,nelec,nsite)
! Written by Cyrus Umrigar
      use constants_mod
      use atom_mod
      use dim_mod
      use pseudo_mod
      use jel_sph2_mod
      use contrl_per_mod
      use orbpar_mod
      use wfsec_mod
      use dorb_mod
      use coefs_mod !GO
      use periodic_mod, only: rlatt_sim

      implicit real*8(a-h,o-z)

! Routine to put electrons down around centers for a VERY crude initial
! configuration if nothing else is available.  It is better to put them
! too close than to put them too far away because they equilibrate faster
! when they are too close.

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /cyldot/ cyldot_v, cyldot_s, cyldot_rho !GO
      common /gndot/ gndot_v0, gndot_rho, gndot_s, gndot_k !GO
      common /wire/ wire_w,wire_length,wire_length2,wire_radius2, wire_potential_cutoff,wire_prefactor,wire_root1

      dimension x(3,*),nsite(*)
      dimension indcoefavailup(nbasis), indcoefavaildn(nbasis) !GO
      dimension nbasup(nbasis), nbasdn(nbasis) !GO

! Loop over spins and centers. If odd number of electrons on all
! atoms then the up-spins have an additional electron.
! So assumption is that system is not strongly polarized.

!     gauss()=dcos(two*pi*rannyu(0))*dsqrt(-two*dlog(rannyu(0)))
      pi=4*datan(1.d0)
      
      if (nloc.eq.-6 .or. nloc.eq.-7) then !GO
        call basis_coef_avail(indcoefavailup, indcoefavaildn) 
        nbasup = 0 
        nbasdn = 0
      endif

      if(iperiodic.eq.3) then ! periodic solids
! Generate points uniformly in lattice coordinates and convert back to cartesian coordinates
        do ielec=1,nelec
          do k=1,ndim
            x(k,ielec)=0
            do i=1,ndim
              r_basis=rannyu(0)-0.5
              x(k,ielec)=x(k,ielec)+rlatt_sim(k,i)*r_basis
            enddo
          enddo
        enddo
        goto 20
      endif

      if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
        if(we.eq.0.d0) stop 'we should not be 0 in sites for quantum dots (nloc=-1)'
       elseif(nloc.eq.-3) then ! jellium RM
        if(zconst.eq.0.d0) stop 'zconst should not be 0 in sites for atoms in jellium (nloc=-3)'
      endif

      ielec=0
      do 10 ispin=1,2
        do 10 i=1,ncent
          if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
            znucc=dsqrt(we)
           elseif(nloc.eq.-3) then ! jellium RM
            znucc=zconst
           elseif(nloc.eq.-4) then ! quantum wire
            znucc=dsqrt(wire_w)
           elseif(nloc.eq.-6) then ! cylindrical quantum dot !GO
            znucc = cyldot_rho
           elseif(nloc.eq.-7) then ! gaussian quantum dot !GO
            znucc = gndot_rho 
           else ! atoms and molecules
            if(znuc(iwctype(i)).eq.0.d0) stop 'znuc should not be 0 in sites for atoms and molecules'
            znucc=znuc(iwctype(i))
          endif
          ju=(nsite(i)+2-ispin)/2
          do 10 j=1,ju
            ielec=ielec+1
            if(ielec.gt.nelec) return
            if(nloc.eq.-1 .or. nloc.eq.-5 .or. nloc.eq.-4) then
              sitsca=1/znucc
             elseif(nloc.eq.-6 .or. nloc.eq.-7) then
              sitsca=znucc/2 !GO
             elseif(j.eq.1) then
              sitsca=1/max(znucc,1.d0)
             elseif(j.le.5) then
              sitsca=2/max(znucc-2,1.d0)
             elseif(j.le.9) then
              sitsca=3/max(znucc-10,1.d0)
             elseif(j.le.18) then
              sitsca=4/max(znucc-18,1.d0)
             else
              sitsca=5/max(znucc-36,1.d0)
            endif


! sample position from exponentials or gaussian around center
! A.D.Guclu 5/2008: need circular coo. for ring shaped quantum dots
            if((nloc.eq.-1 .or. nloc.eq.-5) .and. rring.gt.0.d0) then
              if(ibasis.eq.5) then
                site = (0.5d0 - rannyu(0))/dsqrt(we*oparm(3, iworbd(ielec,1), iwf))
                angle = (0.5d0 - rannyu(0))/dsqrt(oparm(4, iworbd(ielec,1), iwf))
                site = site + oparm(1, iworbd(ielec,1), iwf)
                angle = angle + oparm(2, iworbd(ielec,1), iwf)
!  Make sure electron is near the center of some gaussian - might not work
!     if there's more than 1 slater determinant
                x(1,ielec)=site*dcos(angle)
                x(2,ielec)=site*dsin(angle)
              else
!               This code sampled from a gaussian:
!                site=-dlog(rannyu(0))
!                site=dsqrt(site)
!                site=sign(site,(rannyu(0)-half))
!               This code samples from a smaller, uniform region:
!                site = 2.0d0*(0.5d0 - rannyu(0))
                site = (0.5d0 - rannyu(0))/dsqrt(we)
                angle=2.0d0*pi*(dble(ielec) - rannyu(0))/dble(nelec)
                x(1,ielec)=(site+rring)*dcos(angle)
                x(2,ielec)=(site+rring)*dsin(angle)
!               x(1,ielec)=(sitsca*site+rring)*dcos(angle)
!               x(2,ielec)=(sitsca*site+rring)*dsin(angle)
              endif
! sample position gaussian around center   
! G.Oztarhan 06/2021 (modified: 09/2023):
!      for cylindrical and gaussian quantum dots, make sure that electrons are
!      located around centers of gaussian basis functions within an effective radius.
!      Be careful if there is more than 1 slater determinant
            elseif(nloc.eq.-6 .or. nloc.eq.-7) then
              call find_basis_center(ielec, indcoefavailup, indcoefavaildn, nbasup, nbasdn, iind)
              
              site = dsqrt(-dlog(rannyu(0)))
              if(ibasis.eq.4) then
                 site = site*min(sitsca,1.d0/dsqrt(oparm(3, iind, iwf)))
              else
                 site = site*sitsca
              endif
              angle = 2.0d0*pi*rannyu(0)

              x(1,ielec) = site*dcos(angle) + oparm(1, iind, iwf)
              x(2,ielec) = site*dsin(angle) + oparm(2, iind, iwf)
            else
               do 5 k=1,ndim
! sample position from exponentials or gaussian around center
! a.d.guclu: for wires distribute electrons linearly in y direction
! a.c.mehta: unless floating gaussians, then make sure electrons
!             are close to centers of gaussians
!  Warning:  this might not work if we have multiple slater determinants
                 site=-dlog(rannyu(0))
                 if(nloc.eq.-1 .or. nloc.eq.-4 .or. nloc.eq.-5) site=dsqrt(site)
                 site=sign(site,(rannyu(0)-half))

                 if(nloc.eq.-4) then
                   if (ibasis.eq.6 .or. ibasis.eq.7) then
                     site = (0.5d0 - rannyu(0))/dsqrt(we*oparm(k+2, iworbd(ielec,1), iwf))
!  Make sure electron is near the center of some gaussian - might not work
!     if there's more than 1 slater determinant
                     x(k,ielec) = site + oparm(k, iworbd(ielec,1), iwf)
                   else
                     if(k.eq.2) then
                       x(k,ielec)=sitsca*site+cent(k,i)
                     elseif(iperiodic.eq.0) then
                       x(k,ielec)=wire_length*(0.5d0-rannyu(0))
                     else
                       x(k,ielec)=wire_length*rannyu(0)
                     endif
                   endif
                 else                      ! molecules (mostly at least)
                   x(k,ielec)=sitsca*site+cent(k,i)
                 endif ! nloc.eq.-4
   5           enddo
            endif

   10     continue
      if(ielec.lt.nelec) stop 'bad input to sites'
      
      if (nloc.eq.-6 .or. nloc.eq.-7) then !GO
        nbassum = 0
        do iibas = 1, nbasis
          nbassum = nbassum + nbasup(iibas)
          nbassum = nbassum + nbasdn(iibas)
          if (nbasup(iibas) .gt. 1) stop 'bad input to sites (nbasup(iibas) .gt. 1)'
          if (nbasdn(iibas) .gt. 1) stop 'bad input to sites (nbasdn(iibas) .gt. 1)'
        enddo
        if (nbassum .ne. nelec) stop 'bad input to sites (nbassum .ne. nelec)'
      endif

!     write(6,*)
   20 write(6,'(a,i3,a)') '1 configuration for',nelec,' electrons has been generated by routine sites.'
      write(6,'(''sites:'',1000f12.6)') ((x(k,i),k=1,ndim),i=1,nelec)

      return
      end
      
!-------------------------------------------------------------------------------

      subroutine find_basis_center(ielec, indcoefavailup, indcoefavaildn, nbasup, nbasdn, iind)
      ! Find basis function centers that are available.
      ! 
      ! Try to distribute the electrons uniformly (and randomly) to the location
      ! of the basis functions (e.g. one electron at a site for nelec=ncent=nbasis).
      ! 
      ! Created: Gokhan Oztarhan, 09/2023
      
      use dorb_mod
      use coefs_mod
      use dets_mod
      
      implicit real*8(a-h,o-z)
      
      dimension nbasup(nbasis), nbasdn(nbasis)
      dimension indcoefavailup(nbasis), indcoefavaildn(nbasis)
      dimension indavail(nbasis)
      
      ! Find available basis indices of up and down at the same time
      indcoefavailsame = 0
      do i = 1, nbasis
        if (indcoefavailup(i) .eq. 1 .and. indcoefavaildn(i) .eq. 1) then
          indcoefavailsame = indcoefavailsame + 1
        endif
      enddo
      if (indcoefavailsame .eq. nbasis) then ! if all electrons can be placed around all basis centers
        indcoefavailsame = 0
      endif
      
      ! Find unoccupied basis indices
      indavail = 0
      i = 1
      if (ielec .le. nbasis .and. ielec .le. nup + ndn - indcoefavailsame) then
        do j = 1, nbasis
          if (nbasup(j) .eq. 0 .and. nbasdn(j) .eq. 0) then
            indavail(i) = j
            i = i + 1
          endif
        enddo
      else
        if (ielec .le. nup) then
          do j = 1, nbasis
            if (nbasup(j) .eq. 0) then
              indavail(i) = j
              i = i + 1
            endif
          enddo
        else
          do j = 1, nbasis
            if (nbasdn(j) .eq. 0) then
              indavail(i) = j
              i = i + 1
            endif
          enddo
        endif
      endif

      ! Compare the unoccupied indices to the available basis coef.
      ! If a basis index not available by the coefficients, remove it from
      ! unoccupied basis indices.
      if (ielec .le. nup) then
        do j = 1, i - 1
          if (indcoefavailup(indavail(j)) .eq. 0) then
            indavail(j) = 0
          endif
        enddo
      else
        do j = 1, i - 1
          if (indcoefavaildn(indavail(j)) .eq. 0) then
            indavail(j) = 0
          endif
        enddo
      endif
      
      ! Reorder the unoccupied basis indices
      do i = 1, nbasis
        do j = i + 1, nbasis
          if (indavail(i) .eq. 0 .and. indavail(j) .ne. 0) then
            indtemp = indavail(j)
            indavail(j) = indavail(i)
            indavail(i) = indtemp
          endif
        enddo
      enddo
      
      ! Find non-zero indices
      do i = 1, nbasis
        if (indavail(i) .eq. 0) exit
      enddo
      
      if ((i - 1) .eq. 0) then
        stop 'find_basis_center: bad input to sites (i - 1 = 0)'
      endif

      ! Randomly choose from unoccupied indices
      iind = indavail(int(rannyu(0) * (i - 1)) + 1)
      
      if (ielec .le. nup) then
        nbasup(iind) = nbasup(iind) + 1
      else
        nbasdn(iind) = nbasdn(iind) + 1
      endif

      return
      end
      
!-------------------------------------------------------------------------------

      subroutine basis_coef_avail(indcoefavailup, indcoefavaildn)
      ! Check availability of basis functions by analyzing orbital coefficients.
      ! 
      ! Note: Orbitals are linear combinations of basis functions. The orbital 
      ! coefficients being analyzed are the coefficients of those linear 
      ! combinations.
      ! 
      ! A summation of the orbital coefficients square of a particular basis 
      ! function is over orbital index for a particular spin for all determinants.
      ! If the sum is zero, it may be better not to initialize an electron at the
      ! location of the chosen basis function with that particular spin.
      ! If other basis functions have not any overlap (or a value) at that 
      ! location, Slater determinant may become zero.
      ! 
      ! Created: Gokhan Oztarhan, 09/2023
      
      use dorb_mod
      use coefs_mod
      use dets_mod

      implicit real*8(a-h,o-z)
      
      dimension coefsum(nbasis), coefsum_ratio(nbasis)
      dimension indcoefavailup(nbasis), indcoefavaildn(nbasis)

      coefsum = 0
      do i = 1, nbasis
        do j = 1, ndet
          do k = 1, nup
            coefsum(i) = coefsum(i) + coef(i,iworbd(k,j),1)**2
          enddo
        enddo
      enddo
      
      coefsum_max = maxval(coefsum)
      do i = 1, nbasis
        coefsum_ratio(i) = coefsum(i) / coefsum_max
      enddo
      
      coefsum_threshold = 1d-1
      indcoefavailsum = 0
      do while (indcoefavailsum .lt. nup)
        indcoefavailup = 1
        do i = 1, nbasis
          if (coefsum_ratio(i) .lt. coefsum_threshold) then
            indcoefavailup(i) = 0
          endif
        enddo
        coefsum_threshold = coefsum_threshold * 1d-1
        indcoefavailsum = sum(indcoefavailup)
      enddo

      coefsum = 0
      do i = 1, nbasis
        do j = 1, ndet
          do k = nup + 1, nup + ndn
            coefsum(i) = coefsum(i) + coef(i,iworbd(k,j),1)**2
          enddo
        enddo
      enddo
      
      coefsum_max = maxval(coefsum)
      do i = 1, nbasis
        coefsum_ratio(i) = coefsum(i) / coefsum_max
      enddo
      
      coefsum_threshold = 1d-1
      indcoefavailsum = 0
      do while (indcoefavailsum .lt. ndn)
        indcoefavaildn = 1
        do i = 1, nbasis
          if (coefsum_ratio(i) .lt. coefsum_threshold) then
            indcoefavaildn(i) = 0
          endif
        enddo
        coefsum_threshold = coefsum_threshold * 1d-1
        indcoefavailsum = sum(indcoefavaildn)
      enddo

      return
      end



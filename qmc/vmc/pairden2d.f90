      subroutine pairden2d(p,q,xold,xnew)

! Written by A.D.Guclu jun2005.
! Heavily edited by Gokhan Oztarhan Feb 2022.
! M-fold orientational order parameter implementation added.
      use dets_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      dimension xold(3,nelec),xnew(3,nelec)
      
      logical insideo, insiden
      
      ! Local variables for M-fold orientational order parameter procedure
      real*8 :: Xcm_o, Ycm_o, Xcm_n, Ycm_n
      real*8 :: xprime_o(nelec), yprime_o(nelec)
      real*8 :: xprime_n(nelec), yprime_n(nelec)
      real*8 :: rad_o(nelec), rad_n(nelec)
      integer :: idx_o(nelec), idx_n(nelec)
      real*8 :: C_o, S_o, C_n, S_n
      real*8 :: Theta_body_o, Theta_body_n, theta_j
      real*8 :: xbody_o(nelec), ybody_o(nelec)
      real*8 :: xbody_n(nelec), ybody_n(nelec)
      integer :: n_skip_pd, i, j, k, itemp
      
      ! Local variables for Inter-Ring Phase Correlation
      real*8 :: Cin_o_pd, Sin_o_pd, Cout_o_pd, Sout_o_pd
      real*8 :: Cin_n_pd, Sin_n_pd, Cout_n_pd, Sout_n_pd
      real*8 :: Thetain_o_pd, Thetaout_o_pd, dTheta_o_pd
      real*8 :: Thetain_n_pd, Thetaout_n_pd, dTheta_n_pd
      integer :: Min_pd, Mout_pd, ibin_o_pd, ibin_n_pd
      integer :: n_skip_in_pd, n_skip_out_pd
      real*8 :: pi_pd

      if (M_pd .eq. 0) then
        ! ---------------------------------------------------------
        ! ORIGINAL PAIR DENSITY CODE (Executes if M_pd = 0)
        ! ---------------------------------------------------------
        do ier=1,nelec ! reference electron     
          ! Check if electron in the old config is inside any of the fixed mesh point   
          ix1o = nint(delxi(1) * xold(1,ier)) 
          ix2o = nint(delxi(2) * xold(2,ier)) 
          insideo = .false. 
          do itheta = 1, ithetafix
            if (ix1o .eq. imeshfix1(itheta) .and. ix2o .eq. imeshfix2(itheta))  then
              insideo = .true.
              thetao = -thetafix(itheta)
              exit
            end if
          end do

          ! Check if electron in the new config is inside any of the fixed mesh point  
          ix1n = nint(delxi(1) * xnew(1,ier)) 
          ix2n = nint(delxi(2) * xnew(2,ier)) 
          insiden = .false. 
          do itheta = 1, ithetafix
            if (ix1n .eq. imeshfix1(itheta) .and. ix2n .eq. imeshfix2(itheta))  then
              insiden = .true.
              thetan = -thetafix(itheta)
              exit
            end if
          end do

          ! old config
          if (insideo) then
            do ie2=1,nelec ! electron relative to the reference electron
              if(ie2.ne.ier) then
                call rotate(thetao, xold(1,ie2), xold(2,ie2), x1roto, x2roto)
                ix1roto = nint(delxi(1) * x1roto) 
                ix2roto = nint(delxi(2) * x2roto) 
                
                if (ix1roto .lt. -NAX .or. ix1roto .gt. NAX .or. ix2roto .lt. -NAX .or. ix2roto .gt. NAX) cycle
                if(ier.le.nup) then
                  xx0probut(0,ix1roto,ix2roto)=xx0probut(0,ix1roto,ix2roto)+q
                  if(ie2.le.nup) then
                    xx0probuu(0,ix1roto,ix2roto)=xx0probuu(0,ix1roto,ix2roto)+q
                  else
                    xx0probud(0,ix1roto,ix2roto)=xx0probud(0,ix1roto,ix2roto)+q
                  endif
                else
                  xx0probdt(0,ix1roto,ix2roto)=xx0probdt(0,ix1roto,ix2roto)+q
                  if(ie2.le.nup) then
                    xx0probdu(0,ix1roto,ix2roto)=xx0probdu(0,ix1roto,ix2roto)+q
                  else
                    xx0probdd(0,ix1roto,ix2roto)=xx0probdd(0,ix1roto,ix2roto)+q
                  endif
                endif
              end if
            enddo
          end if
          
          ! new config
          if (insiden) then
            do ie2=1,nelec ! electron relative to the reference electron
              if(ie2.ne.ier) then
                call rotate(thetan, xnew(1,ie2), xnew(2,ie2), x1rotn, x2rotn)
                ix1rotn = nint(delxi(1) * x1rotn) 
                ix2rotn = nint(delxi(2) * x2rotn) 

                if (ix1rotn .lt. -NAX .or. ix1rotn .gt. NAX .or. ix2rotn .lt. -NAX .or. ix2rotn .gt. NAX) cycle
                if(ier.le.nup) then
                  xx0probut(0,ix1rotn,ix2rotn)=xx0probut(0,ix1rotn,ix2rotn)+p
                  if(ie2.le.nup) then
                    xx0probuu(0,ix1rotn,ix2rotn)=xx0probuu(0,ix1rotn,ix2rotn)+p
                  else
                    xx0probud(0,ix1rotn,ix2rotn)=xx0probud(0,ix1rotn,ix2rotn)+p
                  endif
                else
                  xx0probdt(0,ix1rotn,ix2rotn)=xx0probdt(0,ix1rotn,ix2rotn)+p
                  if(ie2.le.nup) then
                    xx0probdu(0,ix1rotn,ix2rotn)=xx0probdu(0,ix1rotn,ix2rotn)+p
                  else
                    xx0probdd(0,ix1rotn,ix2rotn)=xx0probdd(0,ix1rotn,ix2rotn)+p
                  endif
                endif
              end if
            enddo
          end if
          
        enddo ! reference electron
      
      else
        ! ---------------------------------------------------------
        ! M-FOLD BODY-FIXED FRAME CODE (Executes if M_pd /= 0)
        ! ---------------------------------------------------------
        
        ! 1. Translation to Center of Mass
        Xcm_o = 0.d0
        Ycm_o = 0.d0
        Xcm_n = 0.d0
        Ycm_n = 0.d0
        do i = 1, nelec
           Xcm_o = Xcm_o + xold(1,i)
           Ycm_o = Ycm_o + xold(2,i)
           Xcm_n = Xcm_n + xnew(1,i)
           Ycm_n = Ycm_n + xnew(2,i)
        end do
        Xcm_o = Xcm_o / nelec
        Ycm_o = Ycm_o / nelec
        Xcm_n = Xcm_n / nelec
        Ycm_n = Ycm_n / nelec

        do i = 1, nelec
           xprime_o(i) = xold(1,i) - Xcm_o
           yprime_o(i) = xold(2,i) - Ycm_o
           rad_o(i) = dsqrt(xprime_o(i)**2 + yprime_o(i)**2)
           idx_o(i) = i

           xprime_n(i) = xnew(1,i) - Xcm_n
           yprime_n(i) = xnew(2,i) - Ycm_n
           rad_n(i) = dsqrt(xprime_n(i)**2 + yprime_n(i)**2)
           idx_n(i) = i
        end do

        ! 2. Radial Sorting (Bubble Sort)
        do i = 1, nelec - 1
           do j = i + 1, nelec
              if (rad_o(idx_o(i)) .gt. rad_o(idx_o(j))) then
                 itemp = idx_o(i)
                 idx_o(i) = idx_o(j)
                 idx_o(j) = itemp
              end if
              if (rad_n(idx_n(i)) .gt. rad_n(idx_n(j))) then
                 itemp = idx_n(i)
                 idx_n(i) = idx_n(j)
                 idx_n(j) = itemp
              end if
           end do
        end do
        
        ! ---------------------------------------------------------
        ! Inter-Ring Phase Correlation (Executes if 2 or more rings)
        ! ---------------------------------------------------------
        if (nrings_pd .ge. 2) then
           ! Automatically target the two outermost rings
           Mout_pd = conf_pd(nrings_pd)
           Min_pd = conf_pd(nrings_pd - 1)
           
           n_skip_in_pd = 0
           do i = 1, nrings_pd - 2
              n_skip_in_pd = n_skip_in_pd + conf_pd(i)
           end do
           n_skip_out_pd = n_skip_in_pd + Min_pd

           Cin_o_pd = 0.d0; Sin_o_pd = 0.d0
           Cout_o_pd = 0.d0; Sout_o_pd = 0.d0
           Cin_n_pd = 0.d0; Sin_n_pd = 0.d0
           Cout_n_pd = 0.d0; Sout_n_pd = 0.d0

           ! Inner ring
           do k = 1, Min_pd
              j = idx_o(n_skip_in_pd + k)
              theta_j = datan2(yprime_o(j), xprime_o(j))
              Cin_o_pd = Cin_o_pd + dcos(Min_pd * theta_j)
              Sin_o_pd = Sin_o_pd + dsin(Min_pd * theta_j)

              j = idx_n(n_skip_in_pd + k)
              theta_j = datan2(yprime_n(j), xprime_n(j))
              Cin_n_pd = Cin_n_pd + dcos(Min_pd * theta_j)
              Sin_n_pd = Sin_n_pd + dsin(Min_pd * theta_j)
           end do
           Thetain_o_pd = (1.d0 / Min_pd) * datan2(Sin_o_pd, Cin_o_pd)
           Thetain_n_pd = (1.d0 / Min_pd) * datan2(Sin_n_pd, Cin_n_pd)

           ! Outer ring
           do k = 1, Mout_pd
              j = idx_o(n_skip_out_pd + k)
              theta_j = datan2(yprime_o(j), xprime_o(j))
              Cout_o_pd = Cout_o_pd + dcos(Mout_pd * theta_j)
              Sout_o_pd = Sout_o_pd + dsin(Mout_pd * theta_j)

              j = idx_n(n_skip_out_pd + k)
              theta_j = datan2(yprime_n(j), xprime_n(j))
              Cout_n_pd = Cout_n_pd + dcos(Mout_pd * theta_j)
              Sout_n_pd = Sout_n_pd + dsin(Mout_pd * theta_j)
           end do
           Thetaout_o_pd = (1.d0 / Mout_pd) * datan2(Sout_o_pd, Cout_o_pd)
           Thetaout_n_pd = (1.d0 / Mout_pd) * datan2(Sout_n_pd, Cout_n_pd)

           ! Phase difference
           dTheta_o_pd = Thetaout_o_pd - Thetain_o_pd
           dTheta_n_pd = Thetaout_n_pd - Thetain_n_pd

           ! Wrap to [-pi, pi] rigidly
           pi_pd = 4.d0 * datan(1.d0)
           do while (dTheta_o_pd .gt. pi_pd)
              dTheta_o_pd = dTheta_o_pd - 2.d0 * pi_pd
           end do
           do while (dTheta_o_pd .lt. -pi_pd)
              dTheta_o_pd = dTheta_o_pd + 2.d0 * pi_pd
           end do

           do while (dTheta_n_pd .gt. pi_pd)
              dTheta_n_pd = dTheta_n_pd - 2.d0 * pi_pd
           end do
           do while (dTheta_n_pd .lt. -pi_pd)
              dTheta_n_pd = dTheta_n_pd + 2.d0 * pi_pd
           end do

           ! Binning
           ibin_o_pd = int((dTheta_o_pd + pi_pd) / (2.d0 * pi_pd) * NIRBINS_pd) + 1
           if (ibin_o_pd .lt. 1) ibin_o_pd = 1
           if (ibin_o_pd .gt. NIRBINS_pd) ibin_o_pd = NIRBINS_pd
           irphase_pd(ibin_o_pd) = irphase_pd(ibin_o_pd) + q

           ibin_n_pd = int((dTheta_n_pd + pi_pd) / (2.d0 * pi_pd) * NIRBINS_pd) + 1
           if (ibin_n_pd .lt. 1) ibin_n_pd = 1
           if (ibin_n_pd .gt. NIRBINS_pd) ibin_n_pd = NIRBINS_pd
           irphase_pd(ibin_n_pd) = irphase_pd(ibin_n_pd) + p
        end if

        ! Determine core electrons to skip based on conf_pd
        n_skip_pd = 0
        do i = 1, nrings_pd
           if (conf_pd(i) .eq. M_pd) exit
           n_skip_pd = n_skip_pd + conf_pd(i)
        end do
        ! Fallback safety in case of bad input
        if (n_skip_pd + M_pd .gt. nelec) n_skip_pd = 0 

        ! 3. M-fold Orientational Order Parameter & Collective Angle
        C_o = 0.d0
        S_o = 0.d0
        C_n = 0.d0
        S_n = 0.d0
        
        do k = 1, M_pd
           j = idx_o(n_skip_pd + k)
           theta_j = datan2(yprime_o(j), xprime_o(j))
           C_o = C_o + dcos(M_pd * theta_j)
           S_o = S_o + dsin(M_pd * theta_j)

           j = idx_n(n_skip_pd + k)
           theta_j = datan2(yprime_n(j), xprime_n(j))
           C_n = C_n + dcos(M_pd * theta_j)
           S_n = S_n + dsin(M_pd * theta_j)
        end do

        Theta_body_o = (1.d0 / M_pd) * datan2(S_o, C_o)
        Theta_body_n = (1.d0 / M_pd) * datan2(S_n, C_n)

        ! 4. Transform configurations into the Body-Fixed Frame
        do i = 1, nelec
           xbody_o(i) = xprime_o(i) * dcos(-Theta_body_o) - yprime_o(i) * dsin(-Theta_body_o)
           ybody_o(i) = xprime_o(i) * dsin(-Theta_body_o) + yprime_o(i) * dcos(-Theta_body_o)

           xbody_n(i) = xprime_n(i) * dcos(-Theta_body_n) - yprime_n(i) * dsin(-Theta_body_n)
           ybody_n(i) = xprime_n(i) * dsin(-Theta_body_n) + yprime_n(i) * dcos(-Theta_body_n)
        end do

        ! 5. UNCONDITIONAL Binning in the Body-Fixed Frame
        ! No reference particle, no Cartesian box filter. Bin every step!
        do ie2 = 1, nelec
           
           ! Old configuration binning
           ix1roto = nint(delxi(1) * xbody_o(ie2))
           ix2roto = nint(delxi(2) * ybody_o(ie2))

           if (ix1roto .ge. -NAX .and. ix1roto .le. NAX .and. ix2roto .ge. -NAX .and. ix2roto .le. NAX) then
              if (ie2 .le. nup) then
                 xx0probut(0,ix1roto,ix2roto) = xx0probut(0,ix1roto,ix2roto) + q
                 xx0probuu(0,ix1roto,ix2roto) = xx0probuu(0,ix1roto,ix2roto) + q
              else
                 xx0probdt(0,ix1roto,ix2roto) = xx0probdt(0,ix1roto,ix2roto) + q
                 xx0probdd(0,ix1roto,ix2roto) = xx0probdd(0,ix1roto,ix2roto) + q
              endif
           end if

           ! New configuration binning
           ix1rotn = nint(delxi(1) * xbody_n(ie2))
           ix2rotn = nint(delxi(2) * ybody_n(ie2))

           if (ix1rotn .ge. -NAX .and. ix1rotn .le. NAX .and. ix2rotn .ge. -NAX .and. ix2rotn .le. NAX) then
              if (ie2 .le. nup) then
                 xx0probut(0,ix1rotn,ix2rotn) = xx0probut(0,ix1rotn,ix2rotn) + p
                 xx0probuu(0,ix1rotn,ix2rotn) = xx0probuu(0,ix1rotn,ix2rotn) + p
              else
                 xx0probdt(0,ix1rotn,ix2rotn) = xx0probdt(0,ix1rotn,ix2rotn) + p
                 xx0probdd(0,ix1rotn,ix2rotn) = xx0probdd(0,ix1rotn,ix2rotn) + p
              endif
           end if
           
        end do
      end if ! M_pd check

      return
      end

!------------------------------------------------------------------------------------

      subroutine rotate(theta,x1,x2,xrot1,xrot2)

! rotates (x1,x2) by theta. Result is (xrot1,xrot2)

      implicit real*8(a-h,o-z)

      thetarot=datan2(x2,x1)-theta
      r=dsqrt(x1*x1+x2*x2)
      xrot1=r*dcos(thetarot)
      xrot2=r*dsin(thetarot)

      return
      end

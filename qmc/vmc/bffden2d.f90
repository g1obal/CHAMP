      subroutine bffden2d(p,q,xold,xnew)

! Body-Fixed Frame Density utilizing the M-fold Orientational Order Parameter.
! Uses strict Laboratory Frame coordinates (no Center of Mass translation).
! Written by Gokhan Oztarhan, Sep 2026.

      use dets_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      dimension xold(3,nelec),xnew(3,nelec)
      
      real*8 :: p, q
      real*8 :: rad_o(nelec), rad_n(nelec)
      integer :: idx_o(nelec), idx_n(nelec)
      real*8 :: C_o, S_o, C_n, S_n
      real*8 :: Theta_body_o, Theta_body_n, theta_j
      real*8 :: xbody_o(nelec), ybody_o(nelec)
      real*8 :: xbody_n(nelec), ybody_n(nelec)
      integer :: n_skip_bff, i, j, k, itemp, iring, m_shell
      real*8 :: mag_o, mag_n
      
      ! Local variables for Inter-Ring Phase Correlation
      real*8 :: Cin_o_bff, Sin_o_bff, Cout_o_bff, Sout_o_bff
      real*8 :: Cin_n_bff, Sin_n_bff, Cout_n_bff, Sout_n_bff
      real*8 :: Thetain_o_bff, Thetaout_o_bff, dTheta_o_bff
      real*8 :: Thetain_n_bff, Thetaout_n_bff, dTheta_n_bff
      integer :: Min_bff, Mout_bff, ibin_o_bff, ibin_n_bff
      integer :: n_skip_in_bff, n_skip_out_bff
      real*8 :: pi_bff
      
      integer :: ix1roto, ix2roto, ix1rotn, ix2rotn

      ! 1. Lab frame radii (No Center of Mass translation)
      do i = 1, nelec
         rad_o(i) = dsqrt(xold(1,i)**2 + xold(2,i)**2)
         idx_o(i) = i

         rad_n(i) = dsqrt(xnew(1,i)**2 + xnew(2,i)**2)
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
      
      ! 3. Accumulate Order Parameters for ALL valid shells in conf_bff
      n_skip_bff = 0
      do iring = 1, nrings_bff
         m_shell = conf_bff(iring)
         if (m_shell .gt. 1) then
            C_o = 0.d0; S_o = 0.d0
            C_n = 0.d0; S_n = 0.d0
            do k = 1, m_shell
               j = idx_o(n_skip_bff + k)
               theta_j = datan2(xold(2,j), xold(1,j))
               C_o = C_o + dcos(m_shell * theta_j)
               S_o = S_o + dsin(m_shell * theta_j)

               j = idx_n(n_skip_bff + k)
               theta_j = datan2(xnew(2,j), xnew(1,j))
               C_n = C_n + dcos(m_shell * theta_j)
               S_n = S_n + dsin(m_shell * theta_j)
            end do
            mag_o = dsqrt(C_o**2 + S_o**2) / m_shell
            mag_n = dsqrt(C_n**2 + S_n**2) / m_shell
            psi_M_acc_bff(iring) = psi_M_acc_bff(iring) + q * mag_o + p * mag_n
         end if
         n_skip_bff = n_skip_bff + m_shell
      end do

      ! ---------------------------------------------------------
      ! Inter-Ring Phase Correlation (Executes if 2 or more rings)
      ! ---------------------------------------------------------
      if (nrings_bff .ge. 2) then
         ! Automatically target the two outermost rings
         Mout_bff = conf_bff(nrings_bff)
         Min_bff = conf_bff(nrings_bff - 1)
         
         n_skip_in_bff = 0
         do i = 1, nrings_bff - 2
            n_skip_in_bff = n_skip_in_bff + conf_bff(i)
         end do
         n_skip_out_bff = n_skip_in_bff + Min_bff

         Cin_o_bff = 0.d0; Sin_o_bff = 0.d0
         Cout_o_bff = 0.d0; Sout_o_bff = 0.d0
         Cin_n_bff = 0.d0; Sin_n_bff = 0.d0
         Cout_n_bff = 0.d0; Sout_n_bff = 0.d0

         ! Inner ring
         do k = 1, Min_bff
            j = idx_o(n_skip_in_bff + k)
            theta_j = datan2(xold(2,j), xold(1,j))
            Cin_o_bff = Cin_o_bff + dcos(Min_bff * theta_j)
            Sin_o_bff = Sin_o_bff + dsin(Min_bff * theta_j)

            j = idx_n(n_skip_in_bff + k)
            theta_j = datan2(xnew(2,j), xnew(1,j))
            Cin_n_bff = Cin_n_bff + dcos(Min_bff * theta_j)
            Sin_n_bff = Sin_n_bff + dsin(Min_bff * theta_j)
         end do
         Thetain_o_bff = (1.d0 / Min_bff) * datan2(Sin_o_bff, Cin_o_bff)
         Thetain_n_bff = (1.d0 / Min_bff) * datan2(Sin_n_bff, Cin_n_bff)

         ! Outer ring
         do k = 1, Mout_bff
            j = idx_o(n_skip_out_bff + k)
            theta_j = datan2(xold(2,j), xold(1,j))
            Cout_o_bff = Cout_o_bff + dcos(Mout_bff * theta_j)
            Sout_o_bff = Sout_o_bff + dsin(Mout_bff * theta_j)

            j = idx_n(n_skip_out_bff + k)
            theta_j = datan2(xnew(2,j), xnew(1,j))
            Cout_n_bff = Cout_n_bff + dcos(Mout_bff * theta_j)
            Sout_n_bff = Sout_n_bff + dsin(Mout_bff * theta_j)
         end do
         Thetaout_o_bff = (1.d0 / Mout_bff) * datan2(Sout_o_bff, Cout_o_bff)
         Thetaout_n_bff = (1.d0 / Mout_bff) * datan2(Sout_n_bff, Cout_n_bff)

         ! Phase difference
         dTheta_o_bff = Thetaout_o_bff - Thetain_o_bff
         dTheta_n_bff = Thetaout_n_bff - Thetain_n_bff

         ! Wrap to [-pi, pi] rigidly
         pi_bff = 4.d0 * datan(1.d0)
         do while (dTheta_o_bff .gt. pi_bff)
            dTheta_o_bff = dTheta_o_bff - 2.d0 * pi_bff
         end do
         do while (dTheta_o_bff .lt. -pi_bff)
            dTheta_o_bff = dTheta_o_bff + 2.d0 * pi_bff
         end do

         do while (dTheta_n_bff .gt. pi_bff)
            dTheta_n_bff = dTheta_n_bff - 2.d0 * pi_bff
         end do
         do while (dTheta_n_bff .lt. -pi_bff)
            dTheta_n_bff = dTheta_n_bff + 2.d0 * pi_bff
         end do

         ! Binning
         ibin_o_bff = int((dTheta_o_bff + pi_bff) / (2.d0 * pi_bff) * NIRBINS_bff) + 1
         if (ibin_o_bff .lt. 1) ibin_o_bff = 1
         if (ibin_o_bff .gt. NIRBINS_bff) ibin_o_bff = NIRBINS_bff
         irphase_bff(ibin_o_bff) = irphase_bff(ibin_o_bff) + q

         ibin_n_bff = int((dTheta_n_bff + pi_bff) / (2.d0 * pi_bff) * NIRBINS_bff) + 1
         if (ibin_n_bff .lt. 1) ibin_n_bff = 1
         if (ibin_n_bff .gt. NIRBINS_bff) ibin_n_bff = NIRBINS_bff
         irphase_bff(ibin_n_bff) = irphase_bff(ibin_n_bff) + p
      end if

      ! ---------------------------------------------------------
      ! M-fold Target Angle & Coordinate Transformation
      ! ---------------------------------------------------------
      n_skip_bff = 0
      do i = 1, nrings_bff
         if (conf_bff(i) .eq. M_bff) exit
         n_skip_bff = n_skip_bff + conf_bff(i)
      end do
      if (n_skip_bff + M_bff .gt. nelec) n_skip_bff = 0 

      C_o = 0.d0; S_o = 0.d0
      C_n = 0.d0; S_n = 0.d0
      
      do k = 1, M_bff
         j = idx_o(n_skip_bff + k)
         theta_j = datan2(xold(2,j), xold(1,j))
         C_o = C_o + dcos(M_bff * theta_j)
         S_o = S_o + dsin(M_bff * theta_j)

         j = idx_n(n_skip_bff + k)
         theta_j = datan2(xnew(2,j), xnew(1,j))
         C_n = C_n + dcos(M_bff * theta_j)
         S_n = S_n + dsin(M_bff * theta_j)
      end do

      Theta_body_o = (1.d0 / M_bff) * datan2(S_o, C_o)
      Theta_body_n = (1.d0 / M_bff) * datan2(S_n, C_n)

      ! Transform configurations into the Body-Fixed Frame using Lab-Frame origin
      do i = 1, nelec
         xbody_o(i) = xold(1,i) * dcos(-Theta_body_o) - xold(2,i) * dsin(-Theta_body_o)
         ybody_o(i) = xold(1,i) * dsin(-Theta_body_o) + xold(2,i) * dcos(-Theta_body_o)

         xbody_n(i) = xnew(1,i) * dcos(-Theta_body_n) - xnew(2,i) * dsin(-Theta_body_n)
         ybody_n(i) = xnew(1,i) * dsin(-Theta_body_n) + xnew(2,i) * dcos(-Theta_body_n)
      end do

      ! 5. UNCONDITIONAL Binning in the Body-Fixed Frame
      do ie2 = 1, nelec
         ! Old configuration binning
         ix1roto = nint(delxi(1) * xbody_o(ie2))
         ix2roto = nint(delxi(2) * ybody_o(ie2))

         if (ix1roto .ge. -NAX .and. ix1roto .le. NAX .and. ix2roto .ge. -NAX .and. ix2roto .le. NAX) then
            bffden2d_t(ix1roto,ix2roto) = bffden2d_t(ix1roto,ix2roto) + q
            if (ie2 .le. nup) then
               bffden2d_u(ix1roto,ix2roto) = bffden2d_u(ix1roto,ix2roto) + q
            else
               bffden2d_d(ix1roto,ix2roto) = bffden2d_d(ix1roto,ix2roto) + q
            endif
         end if

         ! New configuration binning
         ix1rotn = nint(delxi(1) * xbody_n(ie2))
         ix2rotn = nint(delxi(2) * ybody_n(ie2))

         if (ix1rotn .ge. -NAX .and. ix1rotn .le. NAX .and. ix2rotn .ge. -NAX .and. ix2rotn .le. NAX) then
            bffden2d_t(ix1rotn,ix2rotn) = bffden2d_t(ix1rotn,ix2rotn) + p
            if (ie2 .le. nup) then
               bffden2d_u(ix1rotn,ix2rotn) = bffden2d_u(ix1rotn,ix2rotn) + p
            else
               bffden2d_d(ix1rotn,ix2rotn) = bffden2d_d(ix1rotn,ix2rotn) + p
            endif
         end if
      end do

      return
      end

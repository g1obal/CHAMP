      subroutine den2dwrt(passes,r1ave)
! Written by A.D.Guclu, modified by Cyrus Umrigar for MPI
! routine to print out 2d-density related quantities
! called by finwrt from vmc,dmc,dmc_elec

      use mpi_mod
      use dets_mod
      use contr3_mod
      use pairden_mod
      use fourier_mod
      use contrl_per_mod
      use zigzag_mod
      use const_mod
      implicit real*8(a-h,o-z)
      character*20 file1,file2,file3,file4,file5,file6,file7,file8

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      
! Local variables for BFF and Inter-Ring Phase Correlation writing
      double precision :: range_half_bff, del_phase_bff, angle_bff, alg_weight
      integer :: in1_bff, ir
      
! Normalization variables of pairden
      double precision :: term_u_fix, term_d_fix
      double precision :: w_up, w_dn, total_hits

      if(icoosys.eq.1) then
        del1=1/delxi(1)
        del2=1/delxi(2)
        dely=del2 
        delxt=del1  
        nax1=NAX
        nax2=NAX
      else
        del1=1/delradi
        del2=1/delti
        dely=del1 
        delxt=del2
        nax1=nmeshr
        nax2=nmesht
      endif
      
      if(index(mode,'mov1').ne.0) then
        alg_weight = 1.d0 / dble(nelec)
      else
        alg_weight = 1.d0
      endif
      
      if(ifixe.le.-2) then
        ! Normalization of pair density
        if (pair_hits_u .gt. 0.d0) then
          term_u_fix = 1.d0 / (pair_hits_u * del1 * del2)
        else
          term_u_fix = 0.d0
        endif
        
        if (pair_hits_d .gt. 0.d0) then
          term_d_fix = 1.d0 / (pair_hits_d * del1 * del2)
        else
          term_d_fix = 0.d0
        endif
        
        ! Calculate exact local spin fractions (W_up and W_dn)
        total_hits = pair_hits_u + pair_hits_d
        if (total_hits .gt. 0.d0) then
          w_up = pair_hits_u / total_hits
          w_dn = pair_hits_d / total_hits
        else
          w_up = 0.d0
          w_dn = 0.d0
        endif
      endif

      term=1.d0/(passes*del1*del2)
      
      ! Reference mesh point
      imfix1 = nint(delxi(1) * xfix(1)) !GO
      imfix2 = nint(delxi(2) * xfix(2)) !GO

      if(ifixe.gt.0) then          ! fixed electron pair-density
        if(ifixe.le.nup) then
          file1='pairden_ut'
          file2='pairden_ud'
          file3='pairden_uu'
         else
          file1='pairden_dt'
          file2='pairden_dd'
          file3='pairden_du'
        endif

        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
        endif

! verify the normalization later...
        do in1=-nax1,nax1
          do in2=-nax2,nax2
            write(41,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,den2d_t(in1,in2)*term
            write(42,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,den2d_d(in1,in2)*term
            write(43,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,den2d_u(in1,in2)*term
          enddo
! following spaces are for gnuplot convention:
          write(41,*)
          write(42,*)
          write(43,*)
        enddo
        close(41)
        close(42)
        close(43)
      endif

      if(ifixe.eq.-1 .or. ifixe.eq.-3) then          ! 2d density
        if(index(mode,'vmc').ne.0) then
          file1='den2d_t_vmc'
          file2='den2d_d_vmc'
          file3='den2d_u_vmc'
          file4='pot_t_vmc'
          file5='pot_d_vmc'
          file6='pot_u_vmc'
         else
          file1='den2d_t_dmc'
          file2='den2d_d_dmc'
          file3='den2d_u_dmc'
          file4='pot_t_dmc'
          file5='pot_d_dmc'
          file6='pot_u_dmc'
        endif

        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
          open(44,file=file4,status='unknown')
          open(45,file=file5,status='unknown')
          open(46,file=file6,status='unknown')
         else
!c        open(41,file='/dev/null',status='unknown')
!c        open(42,file='/dev/null',status='unknown')
!c        open(43,file='/dev/null',status='unknown')
!         open(41,file='junk41',status='unknown')
!         open(42,file='junk42',status='unknown')
!         open(43,file='junk43',status='unknown')
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
          open(44,status='scratch')
          open(45,status='scratch')
          open(46,status='scratch')
        endif

! verify the normalization later...
        do in1=-nax1,nax1
          do in2=-nax2,nax2
            write(41,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,den2d_t(in1,in2)*term
            write(42,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,den2d_d(in1,in2)*term
            write(43,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,den2d_u(in1,in2)*term
            write(44,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,pot_ee2d_t(in1,in2)/passes
            write(45,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,pot_ee2d_d(in1,in2)/passes
            write(46,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,pot_ee2d_u(in1,in2)/passes
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
          write(44,*)
          write(45,*)
          write(46,*)
        enddo
        close(41)
        close(42)
        close(43)
        close(44)
        close(45)
        close(46)
      endif

      if(ifixe.le.-2) then          ! full pair-density

! up electron:
        if(nup.gt.0) then
          if(index(mode,'vmc').ne.0) then
            file1='pairden_ut_vmc'
            file2='pairden_ud_vmc'
            file3='pairden_uu_vmc'
          else
            file1='pairden_ut_dmc'
            file2='pairden_ud_dmc'
            file3='pairden_uu_dmc'
          endif

          if(idtask.eq.0) then
            open(41,file=file1,status='unknown')
            open(42,file=file2,status='unknown')
            open(43,file=file3,status='unknown')
           else
            open(41,status='scratch')
            open(42,status='scratch')
            open(43,status='scratch')
          endif
          !do in0=0,NAX !GO
            !if(icoosys.eq.1) then
              !r0=in0*dely
              !if(in0.lt.nint(xfix(1)/dely) .or. in0.gt.nint(xfix(2)/dely)) cycle
              !write(41,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0
              !write(42,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0
              !write(43,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0
              write(41, '("Grid index:",i4,1x,i4,", xfix=",3G20.8E3,", W_up=",G15.6,", W_dn=",G15.6)') imfix1, imfix2, xfix(1), xfix(2), xfix(3), w_up, w_dn
              write(42, '("Grid index:",i4,1x,i4,", xfix=",3G20.8E3,", W_up=",G15.6,", W_dn=",G15.6)') imfix1, imfix2, xfix(1), xfix(2), xfix(3), w_up, w_dn
              write(43, '("Grid index:",i4,1x,i4,", xfix=",3G20.8E3,", W_up=",G15.6,", W_dn=",G15.6)') imfix1, imfix2, xfix(1), xfix(2), xfix(3), w_up, w_dn           
            !else
              !r0=in0*dely + xfix(1)
              !if(r0.gt.xfix(2)) exit
              !write(41,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0-rmean
              !write(42,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0-rmean
              !write(43,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0-rmean
            !endif
            do in1=-NAX,NAX
              do in2=-NAX,NAX
                write(41,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,xx0probut(0,in1,in2)*term_u_fix
                write(42,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,xx0probud(0,in1,in2)*term_u_fix
                write(43,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,xx0probuu(0,in1,in2)*term_u_fix
              enddo
              write(41,*)
              write(42,*)
              write(43,*)
            enddo
          !enddo
          close(41)
          close(42)
          close(43)
        endif

! down electron:
        if(ndn.gt.0) then
          if(index(mode,'vmc').ne.0) then
            file1='pairden_dt_vmc'
            file2='pairden_dd_vmc'
            file3='pairden_du_vmc'
          else
            file1='pairden_dt_dmc'
            file2='pairden_dd_dmc'
            file3='pairden_du_dmc'
          endif

          if(idtask.eq.0) then
            open(41,file=file1,status='unknown')
            open(42,file=file2,status='unknown')
            open(43,file=file3,status='unknown')
           else
            open(41,status='scratch')
            open(42,status='scratch')
            open(43,status='scratch')
          endif

          !do in0=0,NAX
            !if(icoosys.eq.1) then
              !r0=in0*dely
              !if(in0.lt.nint(xfix(1)/dely) .or. in0.gt.nint(xfix(2)/dely)) cycle
              !write(41,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0
              !write(42,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0
              !write(43,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0
              write(41, '("Grid index:",i4,1x,i4,", xfix=",3G20.8E3,", W_up=",G15.6,", W_dn=",G15.6)') imfix1, imfix2, xfix(1), xfix(2), xfix(3), w_up, w_dn
              write(42, '("Grid index:",i4,1x,i4,", xfix=",3G20.8E3,", W_up=",G15.6,", W_dn=",G15.6)') imfix1, imfix2, xfix(1), xfix(2), xfix(3), w_up, w_dn
              write(43, '("Grid index:",i4,1x,i4,", xfix=",3G20.8E3,", W_up=",G15.6,", W_dn=",G15.6)') imfix1, imfix2, xfix(1), xfix(2), xfix(3), w_up, w_dn
            !else
              !r0=in0*dely + xfix(1)
              !if(r0.gt.xfix(2)) exit
              !write(41,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0-rmean
              !write(42,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0-rmean
              !write(43,'(''# Grid point:'',i4,''  r0 ='',G20.8E3)') in0,r0-rmean
            !endif
            do in1=-NAX,NAX
              do in2=-NAX,NAX
                write(41,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,xx0probdt(0,in1,in2)*term_d_fix
                write(42,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,xx0probdd(0,in1,in2)*term_d_fix
                write(43,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,xx0probdu(0,in1,in2)*term_d_fix
              enddo
              write(41,*)
              write(42,*)
              write(43,*)
            enddo
          !enddo
          close(41)
          close(42)
          close(43)
        endif

      endif

      if(ifourier.eq.1 .or. ifourier.eq.3) then          ! 2d fourier t. of the density
        if(index(mode,'vmc').ne.0) then
          file1='fourierrk_t_vmc'
          file2='fourierrk_d_vmc'
          file3='fourierrk_u_vmc'
        else
          file1='fourierrk_t_dmc'
          file2='fourierrk_d_dmc'
          file3='fourierrk_u_dmc'
        endif
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
        endif

! verify the normalization later...
        if(iperiodic.eq.1) then
          term = 1/(passes*dely)
          naxmin = -NAX
          naxmax = NAX
        elseif(rring.eq.0.d0) then
          term=1/(passes*del1)
          naxmin = 1
          naxmax = NAX
        else ! for rings, use same r_min and r_max as density
          term=delradi/passes
          naxmin = -nmeshr
          naxmax = nmeshr
        endif
        do in1=naxmin,naxmax
          if(iperiodic.eq.1) then
            ri = in1*dely
            ri2 = 1.0  ! no need to normalize by 1/r^2 for wires
          elseif(rring.eq.0.d0) then
            ri=(in1*1.d0-0.5d0)*del1
            ri2 = ri*ri
          else
            ri = in1/delradi + rmean
            ri2 = ri*ri
          endif
          do in2=0,nmeshk1
            fk=delk1*in2
            write(41,'(2G20.8E3,G20.8E3)') ri,fk,fourierrk_t(in1,in2)*term/ri2
            write(42,'(2G20.8E3,G20.8E3)') ri,fk,fourierrk_d(in1,in2)*term/ri2
            write(43,'(2G20.8E3,G20.8E3)') ri,fk,fourierrk_u(in1,in2)*term/ri2
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
        enddo
        close(41)
        close(42)
        close(43)
      endif

      if(ifourier.eq.2 .or. ifourier.eq.3) then
        if(index(mode,'vmc').ne.0) then
          file1='fourierkk_t_vmc'
          file2='fourierkk_d_vmc'
          file3='fourierkk_u_vmc'
        else
          file1='fourierkk_t_dmc'
          file2='fourierkk_d_dmc'
          file3='fourierkk_u_dmc'
        endif
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
        endif

! verify the normalization later...
        do in1=-NAK2,NAK2
          do in2=-NAK2,NAK2
            write(41,'(2G20.8E3,G20.8E3)') in1*delk2,in2*delk2,fourierkk_t(in1,in2)/passes
            write(42,'(2G20.8E3,G20.8E3)') in1*delk2,in2*delk2,fourierkk_d(in1,in2)/passes
            write(43,'(2G20.8E3,G20.8E3)') in1*delk2,in2*delk2,fourierkk_u(in1,in2)/passes
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
        enddo
        close(41)
        close(42)
        close(43)
      endif

      if(izigzag.gt.0) then
        if(iperiodic.eq.0) then  !ring
          zzcorr(:) = zzcorr(:) - r1ave*zzcorr1(:) + r1ave*r1ave*zzcorr2(:)
          zzcorrij(:) = zzcorrij(:) - r1ave*zzcorrij1(:) + r1ave*r1ave*zzcorrij2(:)
          yycorr(:) = yycorr(:) - r1ave*yycorr1(:) + r1ave*r1ave*yycorr2(:)
          yycorrij(:) = yycorrij(:) - r1ave*yycorrij1(:) + r1ave*r1ave*yycorrij2(:)
        else if(iperiodic.eq.1) then
          r1ave = 0.d0
        endif
        if(index(mode,'vmc').ne.0) then
          file1='zzcorr_vmc'
          file2='zzcorrij_vmc'
          file5='znncorr_vmc'
          file6='zn2ncorr_vmc'
          file7='yycorr_vmc'
          file8='yycorrij_vmc'
        else
          file1='zzcorr_dmc'
          file2='zzcorrij_dmc'
          file5='znncorr_dmc'
          file6='zn2ncorr_dmc'
          file7='yycorr_dmc'
          file8='yycorrij_dmc'
        endif
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(45,file=file5,status='unknown')
          open(46,file=file6,status='unknown')
          open(47,file=file7,status='unknown')
          open(48,file=file8,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(45,status='scratch')
          open(46,status='scratch')
          open(47,status='scratch')
          open(48,status='scratch')
        endif
        do in2=0,nax2
          zznorm = sum(zzpairden_t(:,in2))
          !if(in2.eq.0) then
          !  write(41,'(G20.8E3,G20.8E3)') in2*delxt,zzcorr(in2)/(zznorm+passes)/delxt
          !else
          !  write(41,'(G20.8E3,G20.8E3)') in2*delxt,zzcorr(in2)/passes/delxt
          !endif
          write(41,'(G20.8E3,G20.8E3)') in2*delxt,zzcorr(in2)/passes/delxt
          write(47,'(G20.8E3,G20.8E3)') in2*delxt,yycorr(in2)/passes/delxt
          write(45,'(G20.8E3,G20.8E3)') in2*delxt*10./dble(nelec),znncorr(in2)/passes
          write(46,'(G20.8E3,G20.8E3)') in2*delxt*10./dble(nelec),zn2ncorr(in2)/passes
        enddo
        do ine = 0,nelec-1
          write(42,'(i8,G20.8E3)') ine,zzcorrij(ine)/passes
          write(48,'(i8,G20.8E3)') ine,yycorrij(ine)/passes
        enddo

        if(izigzag.eq.2) then
          if(index(mode,'vmc').ne.0) then
            file3='zzpairden_t_vmc'
            file4='zzpairdenij_t_vmc'
          else
            file3='zzpairden_t_dmc'
            file4='zzpairdenij_t_dmc'
          endif
          if(idtask.eq.0) then
            open(43,file=file3,status='unknown')
            open(44,file=file4,status='unknown')
           else
            open(43,status='scratch')
            open(44,status='scratch')
          endif

          do in1=-nax1,nax1
            do in2 = -nax2,nax2
              write(43,'(2G20.8E3,G20.8E3)') in1*zzdelyr,in2*delxt,zzpairden_t(in1,in2)/(passes*zzdelyr*delxt)
            enddo
            do ine = 0,nelec-1
              write(44,'(G20.8E3,i8,G20.8E3)') in1*zzdelyr,ine,zzpairdenij_t(in1,ine)/(passes*zzdelyr)
            enddo
            write(43,*)
            write(44,*)
          enddo
          close(43)
          close(44)
        endif
        close(41)
        close(42)
        close(45)
        close(46)
        close(47)
        close(48)
      endif
      
      ! ---------------------------------------------------------
      ! BFF Density and Inter-Ring Phase Correlation Output
      ! ---------------------------------------------------------
      if(ifixe.le.-2 .and. M_bff.gt.0) then
        if(index(mode,'vmc').ne.0) then
          file1='bffden_t_vmc'
          file2='bffden_u_vmc'
          file3='bffden_d_vmc'
          file4='ircorr_vmc'
        else
          file1='bffden_t_dmc'
          file2='bffden_u_dmc'
          file3='bffden_d_dmc'
          file4='ircorr_dmc'
        endif

        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
        else
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
        endif

!       1. Header: Write the parameters and Order Parameters mapped to all parsed rings
        write(41,'("# M_bff: ", i4, ",  nrings_bff: ", i4, ",  conf_bff: ", 20i4)') M_bff, nrings_bff, (conf_bff(ir), ir=1,nrings_bff)
        write(41,'("# M-fold Order Parameters: ", 20G15.6)') (psi_M_acc_bff(ir) * alg_weight / passes, ir=1,nrings_bff)

        write(42,'("# M_bff: ", i4, ",  nrings_bff: ", i4, ",  conf_bff: ", 20i4)') M_bff, nrings_bff, (conf_bff(ir), ir=1,nrings_bff)
        write(42,'("# M-fold Order Parameters: ", 20G15.6)') (psi_M_acc_bff(ir) * alg_weight / passes, ir=1,nrings_bff)

        write(43,'("# M_bff: ", i4, ",  nrings_bff: ", i4, ",  conf_bff: ", 20i4)') M_bff, nrings_bff, (conf_bff(ir), ir=1,nrings_bff)
        write(43,'("# M-fold Order Parameters: ", 20G15.6)') (psi_M_acc_bff(ir) * alg_weight / passes, ir=1,nrings_bff)
        
        ! 2. Density Body (scale the grid term)
        term = term * alg_weight
        
        do in1=-nax1,nax1
          do in2=-nax2,nax2
            write(41,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,bffden2d_t(in1,in2)*term
            write(42,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,bffden2d_u(in1,in2)*term
            write(43,'(2G20.8E3,G20.8E3)') in1*del1,in2*del2,bffden2d_d(in1,in2)*term
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
        enddo

        close(41)
        close(42)
        close(43)

        ! 3. Write Inter-Ring Phase Correlation
        if(nrings_bff .ge. 2) then
          if(idtask.eq.0) then
            open(49,file=file4,status='unknown')
          else
            open(49,status='scratch')
          endif
          
          ! Extract theoretical limits based on M_bff modulo wrapping
          range_half_bff = pi / dble(M_bff)
          del_phase_bff = 2.d0 * range_half_bff / dble(NIRBINS_bff)
          
          do in1_bff = 1, NIRBINS_bff
            angle_bff = -range_half_bff + (dble(in1_bff) - 0.5d0) * del_phase_bff
            ! Divide by del_phase_bff to create a true continuous PDF
            write(49,'(2G20.8E3)') angle_bff, (irphase_bff(in1_bff) * alg_weight) / (passes * del_phase_bff)
          enddo
          
          close(49)
        endif
      endif

      return
      end

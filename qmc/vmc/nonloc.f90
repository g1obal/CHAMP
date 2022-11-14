      subroutine nonloc(x,rshift,rvec_en,r_en,vpsp)
! This routine is called by vmc and dmc. l_do_tmoves=true only if doing dmc and tmoves=true in input.
! Warning: I may still need to update the new tmove parts for the derivatives put in by Julien
! Warning: At present it uses the same points for the reverse T-move to save time, but this may create a tiny bias.
! Warning: Assumes for the reverse T-move probability that the nonlocal ranges of the atoms are non overlapping.
! npotd(ict) is the number of channels, local and nonlocal, for nucleus type ict. In most of the code npotd=npotu and only npotd is used.
! lpotp1=L+1, where L is the local channel.  It would have been more clear to start some arrays at 0, in which case the +1 would not be needed.
! exp(value) is the ratio of the Jastrows at the quadrature points and the initial point.
! w_psi_quad is the wt times the ratio of the wavefns at the quadrature points and the initial point
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use constants_mod
      use control_mod
      use deriv_orb_mod
      use periodic_jastrow_mod !WAS
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use periodic_mod
      use qua_mod
      use contrldmc_mod, only : tau,taunow,tmoves
      use config_dmc_mod, only : xoldw, voldw, psidow, psijow
      use average_mod, only : current_walker
      use variables_mod, only: l_mode_dmc
      use pjase_mod, only: ido_pjas
! Temporary:
      use div_v_dmc_mod
      use delocc_mod
!     use slater_mod, only: detu, detd, slmui, slmdi

      implicit real*8(a-h,o-z)
      integer, parameter :: max_lpotp1=5

      dimension x(3,*),rshift(3,nelec,*),rvec_en(3,nelec,*),r_en(nelec,*)
!    &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension rr_en(nelec,ncent),rr_en2(nelec,ncent),rr_en_sav(ncent),rr_en2_sav(ncent) &
     &,xsav(3),rshift_sav(3,ncent),rvec_en_sav(3,ncent),r_en_sav(ncent)
      dimension xnew(3)               ! local variables for tmoves
      dimension x_tmove_sav(3,nquad*ncent), gpsp_tmove(nquad*ncent), gpsp_tmove_heatbath(nquad*ncent), w_psi_quad(nquad*ncent), pl_costh(max_lpotp1,nquad,nquad)
      integer iwhich_cent_tmove(nquad*ncent)

      if(maxval(npotd).gt.5) stop 'maxval(npotd).gt.max_lpotp1'
      do iq=1,nquad
        do jq=1,iq
          costh=min(1.d0,max(-1.d0,xq(iq)*xq(jq)+yq(iq)*yq(jq)+zq(iq)*zq(jq)))
          do l=1,maxval(npotd)
            pl_costh(l,iq,jq)=yl0(l,costh)
            pl_costh(l,jq,iq)=pl_costh(l,iq,jq)
          enddo
        enddo
      enddo

      do 10 i=1,nelec
        do 10 ic=1,ncent
          call scale_dist(r_en(i,ic),rr_en(i,ic),1)
   10     call scale_dist(r_en(i,ic),rr_en2(i,ic),3)

      if(l_opt_orb_energy) then
        call object_provide_by_index (param_orb_nb_index)
        call object_alloc ('vpsp_ex', vpsp_ex, param_orb_nb)
        vpsp_ex = 0.d0
      endif

      iaccept_tmove=0
      vpsp=0
      do 150 i=1,nelec

! Save position of ith electron and its distances etc. from all nuclei
        do 11 k=1,ndim
   11     xsav(k)=x(k,i)
        do 12 jc=1,ncent
          r_en_sav(jc)=r_en(i,jc)
          rr_en_sav(jc)=rr_en(i,jc)
          rr_en2_sav(jc)=rr_en2(i,jc)
          do 12 k=1,ndim
            rshift_sav(k,jc)=rshift(k,i,jc)
   12       rvec_en_sav(k,jc)=rvec_en(k,i,jc)

! Note if nonlocal parts of psp. are overlapping, ntmove_pts can be multiple of nquad
        ntmove_pts=0
        do 100 ic=1,ncent
          ict=iwctype(ic)

! vps was calculated by calling getvps_xx from nonloc_pot
          iskip=1
          do 15 l=1,npotd(ict)
   15       if(l.ne.lpotp1(ict) .and. dabs(vps(i,ic,l)).gt.1.d-4) iskip=0

          if(iskip.eq.0) then ! If it is within the radius of any of the nonlocal channels

!           if(l_do_tmoves) call rotqua ! If we are doing t-moves then the grid needs to be rotated more often to prevent electrons from being along rays
!           if(l_mode_dmc) call rotqua ! nonloc is never called when l_do_tmoves is true. So the if condition on the previous line is misleading.
            call rotqua ! It is good to choose different rotation for each electron around each nucleus

            ri=one/r_en(i,ic)

            do 60 iq=1,nquad
              ntmove_pts=ntmove_pts+1
              iwhich_cent_tmove(ntmove_pts)=ic
              costh=rvec_en_sav(1,ic)*xq(iq)+rvec_en_sav(2,ic)*yq(iq)+rvec_en_sav(3,ic)*zq(iq)
              costh=costh*ri

              if (l_mode_dmc) then
                quadx(1,ntmove_pts,i,current_walker) = xq(iq)
                quadx(2,ntmove_pts,i,current_walker) = yq(iq)
                quadx(3,ntmove_pts,i,current_walker) = zq(iq)
              endif

              if(iperiodic.eq.0) then
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)
              else
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)+rshift(1,i,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)+rshift(2,i,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)+rshift(3,i,ic)
              endif

! Since we are rotating on sphere around nucleus ic, that elec-nucl distance does not change but distances to other nuclei do
              do 40 jc=1,ncent
                do 38 k=1,ndim
   38             rvec_en(k,i,jc)=x(k,i)-cent(k,jc)
                if(jc.eq.ic) cycle

                if(jc.ne.ic) then
                  if(iperiodic.eq.0) then
                    r_en(i,jc)=0
                    do 39 k=1,ndim
   39                 r_en(i,jc)=r_en(i,jc)+rvec_en(k,i,jc)**2
                    r_en(i,jc)=dsqrt(r_en(i,jc))
                  else
                    call find_image4(rshift(1,i,jc),rvec_en(1,i,jc),r_en(i,jc),rlatt,rlatt_inv)
                  endif

                  call scale_dist(r_en(i,jc),rr_en(i,jc),1)
                  call scale_dist(r_en(i,jc),rr_en2(i,jc),3)
                endif
   40         continue

              iel=i

              electron = iel !JT
              call object_modified_by_index (electron_index) !JT

              call nonlocd(iel,x(1,i),rvec_en,r_en,deter)
              call nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)

              if(ido_pjas.eq.1) then ! periodic Jastrow implemented by WAS
                call nonloc_pjas (iel, x(:,1:nelec), value)
              endif
              !TODO: Update orbe if using fast version of champ

              if(ipr.ge.4) then
                write(6,'(''rr_en,rr_en2'',2d14.6)') rr_en(1,1),rr_en2(1,1)
                write(6,'(''ic,i,iq,deter,value'',3i3,2d14.6)') ic,i,iq,deter,value
              endif

              if (l_mode_dmc) then
                quadr(ntmove_pts,i,current_walker) = &
                deter*exp(value)/psidow(current_walker,1)
              endif

              vpsp_tmove=0
              if(l_do_tmoves) then
                gpsp_tmove(ntmove_pts)=0
                x_tmove_sav(1:3,ntmove_pts)=x(1:3,i)
                if(ipr.ge.1) write(6,'(''iw,ielec='',9i5)') current_walker, i
              endif
              do 50 l=1,npotd(ict)
                if(l.ne.lpotp1(ict)) then
                  vpsp_tmove=vpsp_tmove + vps(i,ic,l)*yl0(l,costh)
                  if(l_do_tmoves) then
                    if(ipr.ge.1) write(6,'(''vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),tau'',9d12.4)') &
     &                vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),tau
!                   gpsp_tmove(ntmove_pts)=gpsp_tmove(ntmove_pts) + vps(i,ic,l)*yl0(l,costh)                  ! linear approx in Casula
                    gpsp_tmove(ntmove_pts)=gpsp_tmove(ntmove_pts) + (exp(-tau*vps(i,ic,l))-1)*yl0(l,costh) ! exact
                    call systemflush(6)
                  endif
                  if(ipr.ge.1) write(6,'(''nonloc: l,yl0(l,costh),deter,exp(value),yl0(l,costh)*deter*exp(value)'',i3,9f20.15)') &
     &            l,yl0(l,costh),deter,exp(value),yl0(l,costh)*deter*exp(value)
                endif
   50         continue ! npotd(ict)
              factor=wq(iq)*exp(value)
              vpsp=vpsp+factor*deter*vpsp_tmove
              if(l_do_tmoves) then
                w_psi_quad_ex=factor/psidow(current_walker,1)
                w_psi_quad(ntmove_pts)=w_psi_quad_ex*deter
!               gpsp_tmove(ntmove_pts)=w_psi_quad(ntmove_pts)*tau*(-gpsp_tmove(ntmove_pts)) ! linear approx in Casula
                gpsp_tmove(ntmove_pts)=w_psi_quad(ntmove_pts)*gpsp_tmove(ntmove_pts)           ! exact
                if(ipr.ge.-1) then
                  write(6,'(''i, ic, ntmove_pts, current_walker, gpsp_tmove, -vps*yl0, wf_ratio='',4i4,9es12.4)') i, ic, ntmove_pts, current_walker, gpsp_tmove(ntmove_pts), gpsp_tmove(ntmove_pts)/((wq(iq)*tau*deter/psidow(current_walker,1))*exp(value)), wq(iq)*tau*(deter/psidow(current_walker,1))*exp(value)
                  call systemflush(6)
                endif
              endif

! JT          For singly-excited wave functions
              if(l_opt_orb_energy) then
                call object_provide_by_index (psid_ex_in_x_index)
                do iex = 1, param_orb_nb
                  vpsp_ex(iex)=vpsp_ex(iex)+w_psi_quad_ex*psid_ex_in_x(iex)*vpsp_tmove
                enddo
              endif

   60       continue ! nquad

! Restore electron i and its distances to value it had before doing angular integration
            do 68 k=1,ndim
   68         x(k,i)=xsav(k)
            do 70 jc=1,ncent
              r_en(i,jc)=r_en_sav(jc)
              rr_en(i,jc)=rr_en_sav(jc)
              rr_en2(i,jc)=rr_en2_sav(jc)
              do 70 k=1,ndim
                rshift(k,i,jc)=rshift_sav(k,jc)
   70           rvec_en(k,i,jc)=rvec_en_sav(k,jc)

          endif  ! iskip, i.e. electron within r_c for this nucleus
  100   continue ! ncent

! Warning: For the moment correlated sampling for forces does not work with t-moves.
! Decide which if any tmove to perform for electron i
        if(l_do_tmoves .and. ntmove_pts.ne.0) then
          gpsp_tmove_sum=1
          do itmove_pts=1,ntmove_pts
            gpsp_tmove_sum=gpsp_tmove_sum+max(gpsp_tmove(itmove_pts),0.d0)
            if(itmove_pts.eq.1) then
              gpsp_tmove_heatbath(itmove_pts)=1
            else
              gpsp_tmove_heatbath(itmove_pts)=gpsp_tmove_heatbath(itmove_pts-1)+max(gpsp_tmove(itmove_pts-1),0.d0)
            endif
          enddo
          if(ipr.ge.1) write(6,'(''gpsp_tmove_sum, gpsp_tmove_heatbath'',100f10.6)') gpsp_tmove_sum, (gpsp_tmove_heatbath(ii),ii=1,ntmove_pts)

          if(abs(gpsp_tmove_sum-1).ge.1.d-6) random=rannyu(0)
          if(gpsp_tmove_heatbath(1).le.random*gpsp_tmove_sum) then ! Do t-move
            iwhich_tmove=ntmove_pts ! Last quadrature point is chosen if we make it through loop without satisfying the if within the loop.
            do itmove_pts=2,ntmove_pts
              if(gpsp_tmove_heatbath(itmove_pts).gt.random*gpsp_tmove_sum) then
                iwhich_tmove=itmove_pts-1
                exit
              endif
            enddo

! Starting from the T-move point, compute reverse proposal probability
            ic=iwhich_cent_tmove(iwhich_tmove)      ! which nucleus is t-move done around
            ict=iwctype(ic)                         ! which nucleus type
            iq=1+mod(iwhich_tmove-1,nquad)          ! which quadrature point on nucleus ic
            ioffset=((iwhich_tmove-1)/nquad)*nquad  ! offset for quadrature points on this nucleus
            gpsp_tmove_sum_rev=1
!           write(6,'(''ic,ict,nquad,npotd(ict),lpotp1(ict),iwhich_tmove,ioffset'',9i4)') ic,ict,nquad,npotd(ict),lpotp1(ict),iwhich_tmove,ioffset
            do jq=1,nquad
              psi_ratio=(wq(iq)*w_psi_quad(ioffset+jq))/(wq(jq)*w_psi_quad(iwhich_tmove))
              gpsp=0d0
              do l=1,npotd(ict)
                if(l.ne.lpotp1(ict)) then
!                 gpsp_tmove_sum_rev=gpsp_tmove_sum_rev + max(-wq(jq)*psi_ratio*tau*vps(i,ic,l)*pl_costh(l,iq,jq),0.d0)          ! linear approx in Casula
!                  gpsp_tmove_sum_rev=gpsp_tmove_sum_rev + max(wq(jq)*psi_ratio*(exp(-tau*vps(i,ic,l))-1)*pl_costh(l,iq,jq),0.d0) ! exact
                  gpsp=gpsp + wq(jq)*psi_ratio*(exp(-tau*vps(i,ic,l))-1)*pl_costh(l,iq,jq) ! exact
!           write(6,'(''jq,l,wq(jq),psi_ratio,(exp(-tau*vps(i,ic,l))-1),pl_costh(l,iq,jq)'',2i4,9es12.4)') jq,l,wq(jq),psi_ratio,(exp(-tau*vps(i,ic,l))-1),pl_costh(l,iq,jq)
                endif
              enddo
              gpsp_tmove_sum_rev=gpsp_tmove_sum_rev+max(0d0,gpsp)
            enddo

            p=min(1.d0,gpsp_tmove_sum/gpsp_tmove_sum_rev)  ! acceptance probability for our T-move algorithm
!           write(6,'(''ic,ict,iq,gpsp_tmove_sum,gpsp_tmove_sum_rev,p'',3i4,3es12.4)') ic,ict,iq,gpsp_tmove_sum,gpsp_tmove_sum_rev,p
            if(rannyu(0).lt.p .or. tmoves_type=='casula1') then ! In Casula's algorithm moves are always accepted and causes large time-step error
              iaccept_tmove=1
              xnew(1:3)=x_tmove_sav(1:3,iwhich_tmove)
              xoldw(1:3,i,current_walker,1)=x_tmove_sav(1:3,iwhich_tmove)
              call hpsiedmc(i,current_walker,xnew,psidow(current_walker,1),psijow(current_walker,1),voldw(1,1,current_walker,1))
              call systemflush(6)
              call jassav(i)
              if(ibasis.eq.3) then                        ! complex calculations
                call cdetsav(i)
              else
                call detsav(i)
              endif
            endif ! accept/reject for t-move

          endif ! do tmove or not for this electron

        endif ! l_do_tmoves
        call systemflush(6)

  150 continue ! nelec

! Note that above we called hpsiedmc but not hpsi.  So, the velocity is updated after the tmove but not the energy.
! So, in the reweighting we are using the energies from before and after the drift-dfus-acc/rej step and not after tmove.
! This is slightly different from what is presented on pg. 5 of the Casula et al. paper, if I take it literally.
      if(l_do_tmoves) then !FP
        if(ibasis.eq.3) then                     !complex basis set
          call cwalksav_det(current_walker)
        else
          call walksav_det(current_walker)
        endif
        call walksav_jas(current_walker)
!     else Warning to be done?
      endif

      call object_modified_by_index (vpsp_ex_index) ! JT

!      if (l_do_tmoves.and.(current_walker.eq.100)) stop 'stop after tmoves' !TA

      return
      end
!-----------------------------------------------------------------------

      function yl0(l,costh)
! (2L+1)*P_L(costh)
! (2L+1)*P_L(costh)
! This is not quite Y_L0 but sqrt[4pi*(2L+1)]Y_L0 = (2L+1)*P_L(costh)
! Note that the associated P_L^m and the unassociated P_L Legendre polynomials are the same for m=0.
! l is actually L+1.

      implicit real*8(a-h,o-z)

      if(l.eq.1) then
        yl0=1.d0
       elseif(l.eq.2) then
        yl0=3.d0*costh
       elseif(l.eq.3) then
        yl0=2.5d0*(3*costh*costh-1)
       elseif(l.eq.4) then
        yl0=3.5d0*costh*(5*costh*costh-3)
       elseif(l.eq.5) then
        yl0=1.125d0*(35*costh**4-30*costh**2+3)
       else
        stop 'yl0 implemented to l=4 only (Warning: l is l+1)'
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine nonlocd(iel,x,rvec_en,r_en,determ)
! Written by Claudia Filippi, modified by Cyrus Umrigar
      use all_tools_mod
      use control_mod
      use eloc_mod
      use dorb_mod
      use slatn_mod
      use orbe_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contr3_mod
      use phifun_mod
      use const_mod
      use slatn2_mod
      use slater_mod, only: detu, detd, slmui, slmdi
      implicit real*8(a-h,o-z)

      dimension x(3),rvec_en(3,nelec,*),r_en(nelec,*)
!    &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension ratio(ndet)

!     determ=0

! get orbitals for all electron iel
! Note x is 3-dim but rvec_en,r_en have info about all electrons.
! So, iel is passed to select elements of rvec_en,r_en,phin and for IO.
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call orbitals_loc_anae(iel,rvec_en,r_en,orbe)
         else
          call orbitals_loc_nume(x,orbe)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pwe(iel,x,orbe)
         else
          call orbitals_period_nume(x,orbe)
        endif

      endif

      call object_modified_by_index (orbe_index) !JT

      if(iel.le.nup) then

        ikel=nup*(iel-1)

        do idet=1,ndetup
           ratio(idet)=0
           do j=1,nup
              ratio(idet)=ratio(idet)+slmui(j+ikel,idet)*orbe(iworbdup(j,idet))
           enddo
        enddo

        if(.not. l_opt_exp) then
           do idet=1,ndetup
              detn(idet)=detu(idet)*ratio(idet)
           enddo
        else
           do idet=1,ndetup
              detn(idet)=detu(idet)*ratio(idet)
              do i=1,nup
                 if(i.ne.iel) then
                    ik=nup*(i-1)
                    sum=0
                    do j=1,nup
                       sum=sum+slmui(j+ik,idet)*orbe(iworbdup(j,idet))
                    enddo
                    sum=sum/ratio(idet)
                    do j=1,nup
                       slmin(j+ik,idet)=slmui(j+ik,idet)-slmui(j+ikel,idet)*sum
                    enddo
                 endif
              enddo
              do j=1,nup
                 slmin(j+ikel,idet)=slmui(j+ikel,idet)/ratio(idet)
              enddo
           enddo
        endif

      else

        ikel=ndn*(iel-nup-1)

        do idet=1,ndetdn
           ratio(idet)=0
           do j=1,ndn
              ratio(idet)=ratio(idet)+slmdi(j+ikel,idet)*orbe(iworbddn(j,idet))
           enddo
        enddo
        if(.not. l_opt_exp) then
           do idet=1,ndetdn
              detn(idet)=detd(idet)*ratio(idet)
           enddo
        else
           do idet=1,ndetdn
              detn(idet)=detd(idet)*ratio(idet)
              do i=1,ndn
                 if(i+nup.ne.iel) then
                    ik=ndn*(i-1)
                    sum=0
                    do j=1,ndn
                       sum=sum+slmdi(j+ik,idet)*orbe(iworbddn(j,idet))
                    enddo
                    sum=sum/ratio(idet)
                    do j=1,ndn
                       slmin(j+ik,idet)=slmdi(j+ik,idet)-slmdi(j+ikel,idet)*sum
                    enddo
                 endif
              enddo
              do j=1,ndn
                 slmin(j+ikel,idet)=slmdi(j+ikel,idet)/ratio(idet)
              enddo
           enddo
        endif

      endif

      call object_modified_by_index (detn_index) !JT
      call object_modified_by_index (slmin_index) !fp

      determ=0
      do 115 icsf=1,ncsf
        do 115 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
          if(iel.le.nup) then
            term=detn(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            term=detu(iwdetup(idet))*detn(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
  115     determ=determ+term

! Derivatives wrt to csf_coefs for optimizing them
      if(index(mode,'fit').ne.0 .or. igradhess.ge.1 .or. l_opt_csf) then
        do 140 iparm=1,nparmcsf
          icsf=iwcsf(iparm)
          deti_new(iparm)=0
          do 140 idet_in_csf=1,ndet_in_csf(icsf)
            idet=iwdet_in_csf(idet_in_csf,icsf)
            if(iel.le.nup) then
              term=detn(iwdetup(idet))*detd(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)
             else
              term=detu(iwdetup(idet))*detn(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)
            endif
  140       deti_new(iparm)=deti_new(iparm)+term
      endif

      return
      end
!-----------------------------------------------------------------------

!WAS  subroutine nonlocj(iel,x,rshift,rr_en,rr_en2,value)
      subroutine nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)
! fso are pair contributions in the Jastrow exponent before electron iel moves
! fsn are pair contributions in the Jastrow exponent after  electron iel moves
! value=sum(fsn-fso), i.e., the change in the Jastrow exponent when electron iel moves
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use control_mod
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contrl_per_mod
      use jaspar_mod
      use bparm_mod
      use periodic_mod
      use jaso_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,nelec,ncent),rr_en(nelec,ncent),rr_en2(nelec,ncent),fsn(nelec,nelec),dx(3)

      dimension r_en(nelec,ncent)

      fsumn=0

      if(nelec.lt.2) goto 47

      do 45 jj=1,nelec

        if(jj.eq.iel) goto 45
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif

        sspinn=1
        ipar=0
        isb=1
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2b.eq.2) then
            isb=2
           elseif(nocuspb.eq.0) then
            sspinn=half
          endif
          ipar=1
        endif

        do 10 k=1,ndim
   10     dx(k)=x(k,jj)-x(k,iel)

        if(iperiodic.eq.0) then
          rij=0
          do 20 k=1,ndim
   20       rij=rij+dx(k)**2
          rij=dsqrt(rij)
         else
          call find_image3(dx,rij,rlatt_sim,rlatt_sim_inv)
        endif

! e-e terms
        call scale_dist(rij,u,2)

        fsn(i,j)=psibnl(u,isb,ipar)

! e-e-n terms
! The scaling is switched in psinl, so do not do it here.
!     if(isc.ge.12) call scale_dist(rij,u,3)
      call scale_dist(rij,u,4)

        do 40 ic=1,ncent
          it=iwctype(ic)
!WAS 40   fsn(i,j)=fsn(i,j) + psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
   40     fsn(i,j)=fsn(i,j) + psinl(u,rshift(1,i,ic),rshift(1,j,ic),r_en(i,ic),r_en(j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
        fsumn=fsumn+fsn(i,j)-fso(i,j)

   45 continue

! e-n terms
   47 fsn(iel,iel)=0
      do 50 ic=1,ncent
        it=iwctype(ic)
   50   fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en(iel,ic),it)

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)

      value=fsumn

      return
      end

      subroutine pairden2d(p,q,xold,xnew)

! Written by A.D.Guclu jun2005.
! Heavily edited by Gokhan Oztarhan Feb 2022.

      use dets_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      dimension xold(3,nelec),xnew(3,nelec)
      
      logical insideo, insiden
      
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
          if (ier .le. nup) then
            pair_hits_u = pair_hits_u + q
          else
            pair_hits_d = pair_hits_d + q
          endif
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
          if (ier .le. nup) then
            pair_hits_u = pair_hits_u + p
          else
            pair_hits_d = pair_hits_d + p
          endif
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
        
      enddo

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

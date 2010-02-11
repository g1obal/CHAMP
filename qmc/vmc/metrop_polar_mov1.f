      subroutine metrop_polar_mov1(ipass)
c     subroutine metrop5(ipass)
c Written by Cyrus Umrigar
c Uses the accelerated Metropolis method described in:
c 1) Accelerated Metropolis Method, C.J. Umrigar, PRL 71, 408 (1993).
c 2) Variational Monte Carlo Basics and Applications to Atoms and Molecules,
c    C.J. Umrigar, in {\it Quantum Monte Carlo Methods in Physics and Chemistry},
c    edited by M.~P. Nightingale and C.~J. Umrigar. NATO ASI Series, Series C,
c    Mathematical and Physical Sciences, Vol. C-525,
c    (Kluwer Academic Publishers, Boston, 1999)
      use all_tools_mod
      use constants_mod
      use contrl_per_mod, only: iperiodic
      use atom_mod
      use config_mod
      use dets_mod
      use const_mod
      use const2_mod
      use dim_mod
      use forcepar_mod
!      use doefp_mod
      use delocc_mod
      use denupdn_mod
      use stepv_mod
      use kinet_mod
      use estsig_mod
      use estsum_mod
      implicit real*8(a-h,o-z)

      parameter (d3b2=1.5d0)
      parameter (eps=1.d-10)
      parameter (g3b2=.886226925452758d0)
c     parameter (g5b2=1.329340388179137d0)
c g3b2, g5b2 are gamma3/2), gamma(5/2)

c The foll. additions have been made:
c 1) Linear and Exponential forms of Tij.
c 2) Make theta_max a function of r

c The foll. addtions still need to be made:
c 1) Quadratic, gaussian, Morse and Exp(-zeta*r)+co*Exp(-r) forms of Tij
c    Last 2 are prob. best. Slater has been done and si minor improvement
c 2) Generalize to molecules. This requires geometric rejections.

      common /stats_vmc/ rejmax


      dimension xaxis(3),yaxis(3),zaxis(3),div_vn(nelec)

      area(ri,r1,r2,v)=dabs((one/sqrt(ri))*
     &(r2**d3b2*(two*(one-v*ri)/three+.4d0*v*r2)
     &-r1**d3b2*(two*(one-v*ri)/three+.4d0*v*r1)))

      deltri=one/deltar
c The transition probability is an approximation to psi(new)/psi(old)
c and is positive definite.

      do 300 i=1,nelec
        fxop=one
c Although rmino is saved, recalculate it, otherwise there is a build
c up of errors via zaxis
        rmino(i)=dsqrt(rvmino(1,i)**2+rvmino(2,i)**2+rvmino(3,i)**2)
c Choose lower and upper values of r sampling
        rbot=rmino(i)*deltri
        rtop=rmino(i)*deltar
c Calculate magnitude of the velocity in the radial direction
        voldr=zero
        do 10 k=1,ndim
   10     voldr=voldr+vold(k,i)*rvmino(k,i)
        voldr=voldr/rmino(i)
        voldr=max(voldr,-2*znuc(iwctype(1)))
c Place x-axis along direction of angular change and
c Calculate the velocity in the phi direction
        voldp=zero
        do 20 k=1,ndim
          zaxis(k)=rvmino(k,i)/rmino(i)
          xaxis(k)=vold(k,i)-voldr*zaxis(k)
   20     voldp=voldp+xaxis(k)**2
        voldp=dsqrt(voldp)
        if(voldp.lt.eps) then
          xaxis(1)=eps*(one-zaxis(1)**2)
          xaxis(2)=eps*(-zaxis(1)*zaxis(2))
          xaxis(3)=eps*(-zaxis(1)*zaxis(3))
          voldp=eps*sqrt(one-zaxis(1)**2)
        endif
        do 30 k=1,ndim
   30     xaxis(k)=xaxis(k)/voldp

c y-axis is cross-product of z and x axes
        yaxis(1)=zaxis(2)*xaxis(3)-zaxis(3)*xaxis(2)
        yaxis(2)=zaxis(3)*xaxis(1)-zaxis(1)*xaxis(3)
        yaxis(3)=zaxis(1)*xaxis(2)-zaxis(2)*xaxis(1)

c Temporary test of fbias
        voldr=voldr*fbias
        voldp=voldp*fbias

        if(voldr.ge.zero) then
c         write(6,*) 'lo'
c Use linear approx for radial fn
c Determine the maximum value of radial function for rejection sampling
          rmax=third*(rmino(i)-one/voldr)
c         rmax=-one/voldr
          if(rmax.lt.rbot) rmax=rbot
          if(rmax.gt.rtop) rmax=rtop
          fmax=sqrt(rmax)*abs(one+voldr*(rmax-rmino(i)))
c         fmax=sqrt(rmax/rmino(i))*dexp(voldr*(rmax-rmino(i)))
c         fmax=sqrt(rmax)*dexp(voldr*rmax)
          fmax2=sqrt(rtop)*abs(one+voldr*(rtop-rmino(i)))
          fmax=max(fmax,fmax2)
c   Sample sqrt(r_f/r_i)*(1+v*(r_f-r_i)) by rejection
   40     rtry=((deltar-deltri)*rannyu(0)+deltri)*rmino(i)
            ratio=sqrt(rtry)*abs(one+voldr*(rtry-rmino(i)))/fmax
            rejmax=max(rejmax,ratio)
            if(ratio.gt.rannyu(0)) goto 50
c           if(sqrt(rtry/rmino(i))*dexp(voldr*(rtry-rmino(i)))/fmax.gt.
c    &      rannyu(0)) goto 50
c           if(sqrt(rtry)*dexp(voldr*rtry)/fmax.gt.rannyu(0)) goto 50
          goto 40
   50     fxop=fxop*sqrt(rtry/rmino(i))*abs(one+voldr*(rtry-rmino(i)))

c   Calculate the integral of T
          rzero=rmino(i)-one/voldr
          if(rzero.lt.rbot .or. rzero.gt.rtop) then
            areao=area(rmino(i),rbot,rtop,voldr)
           else
            areao=area(rmino(i),rbot,rzero,voldr) +
     &      area(rmino(i),rzero,rtop,voldr)
          endif

         else
c Use exponential approx for radial fn
c Determine the maximum value of radial function for rejection sampling
          rmax=-half/voldr
          if(rmax.lt.rbot) rmax=rbot
          if(rmax.gt.rtop) rmax=rtop
          fmax=sqrt(rmax)*exp(voldr*rmax)
c         fmax=sqrt(rmax/rmino(i))*dexp(voldr*(rmax-rmino(i)))
c         fmax=sqrt(rmax)*dexp(voldr*rmax)
          fmax2=sqrt(rtop)*exp(voldr*rtop)
          fmax=max(fmax,fmax2)
c   Sample r*exp(v*r) by rejection
cc   Warning: need to put in deltar
   42     rtry=((deltar-deltri)*rannyu(0)+deltri)*rmino(i)
            ratio=sqrt(rtry)*exp(voldr*rtry)/fmax
            rejmax=max(rejmax,ratio)
            if(ratio.gt.rannyu(0)) goto 52
          goto 42
   52     fxop=fxop*sqrt(rtry/rmino(i))*exp(voldr*(rtry-rmino(i)))

c   Calculate the integral of T
          zrbot=-voldr*rbot
          zrtop=-voldr*rtop
          zebot=dsqrt(zrbot**3)*dexp(-zrbot)
          zetop=dsqrt(zrtop**3)*dexp(-zrtop)
          g32bot=gammai(d3b2,zrbot,zebot,iflagb)
          g32top=gammai(d3b2,zrtop,zetop,iflagt)
          if(iflagb*iflagt.eq.1) then
            g32dif=g32top-g32bot
           else
            g32dif=g32top-g32bot+g3b2
          endif
          areao=g32dif
     &    /(dexp(voldr*rmino(i))*(dsqrt(rmino(i)*(-voldr)**3)))

c Warning: The replacemnt of the foll 2 lines by the above lines has not
c yet been tested
c         areao=(gammai(d3b2,zrtop,zetop)-gammai(d3b2,zrbot,zebot))
c    &    /(dexp(voldr*rmino(i))*(dsqrt(rmino(i)*(-voldr)**3)))

        endif

c Sample cos(theta)
        raver=half*(rmino(i)+rtry)
        deltt=deltat+(two-deltat)/(one+25*raver*raver)
        costht=one-deltt*rannyu(0)
c       costht=one-deltt*half*(one+rannyu(0))
        zprime=rtry*costht
c       theta=dacos(costht)
c       sintht=dsin(theta)
        sintht=dsqrt(one-costht*costht)
c Truncate phi variation if it goes through zero
c Sample phi by rejection. Note it is OK to have a) theta or sin(theta)
c and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
c are consistent
c       term=dmin1(voldp*rtry*theta,one)
c       term=dmin1(voldp*rtry*sintht,one)
c       term=dmin1(voldp*raver*sintht,one)
        term=voldp*raver*sintht
        fmax=one+term
   60   phitry=pi*rannyu(0)
c         if((one+term*dcos(phitry))/fmax.gt.rannyu(0)) goto 70
          if(abs(one+term*dcos(phitry))/fmax.gt.rannyu(0)) goto 70
        goto 60
c  70   fxop=fxop*(one+term*dcos(phitry))
   70   fxop=fxop*abs(one+term*dcos(phitry))
        if(term.gt.one) areao=areao*(one+(term-one)**2/(two*term))
c Calculate x and y coordinates in local coordinate system
        xprime=rtry*sintht*dcos(phitry)
        yprime=dsqrt(rtry*rtry-xprime**2-zprime**2)
        if(rannyu(0).lt.half) yprime=-yprime

c Convert back to original coordinate system
        do 80 k=1,ndim
          xnew(k,i)=xaxis(k)*xprime+yaxis(k)*yprime+zaxis(k)*zprime
   80     rvminn(k,i)=xnew(k,i)
        rminn(i)=rtry

        if(ipr.ge.1) then
        rtest=dsqrt(rvminn(1,i)**2+rvminn(2,i)**2+rvminn(3,i)**2)
        rtest2=dsqrt(xprime**2+yprime**2+zprime**2)
        write(6,'(''rtest,rtest2,rtry'',9d14.6)')rtest,rtest2,rtry,
     &  rtest-rtry,rtest2-rtry
        write(6,'(''vold='',9d12.4)') (vold(k,i),k=1,ndim)
        write(6,'(''voldr,voldp='',9d12.4)') voldr,voldp
        write(6,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
        write(6,'(''rmino(i),rmax,rzero'',9f9.4)') rmino(i),rmax,rzero
        write(6,'(''rtry,costht,sintht,phitry'',9f9.4)') rtry,costht,
     &  sintht,phitry
        write(6,'(''fxop'',9f9.4)') fxop
        endif

c Write warning msg. if electron is going far away
        if(rminn(i).gt.100.d0 .and. ndim.eq.3 .and. iperiodic.eq.0) then
          write(6,'(''Warning: rminn(i) too large, i, rminn(i) ='',i4,d12.4)') i,rminn(i)
          write(6,'(''Warning: xold,xnew='',9es12.4)') (xold(k,i),k=1,ndim),(xnew(k,i),k=1,ndim)
          if(rminn(i).gt.1000.d0) stop 'rminn(i) too large'
        endif

c calculate psi etc. at new configuration

      call hpsi(xnew,psidn,psijn,vnew,div_vn,d2,pen,pein,enew(1),denergy,1)
      psi2n(1)=2*(dlog(dabs(psidn))+psijn)

c form the Jackson Feenberg kinetic energy expression

      tjfn=-d2*half*hb

c If error is large then save config. to use in optimizing routine

c     if(dabs((enew(1)-etrial)/etrial).gt.1.d0.and.iwrit8.le.1000) then
c        iwrit8=iwrit8+1
c        write(8,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
c    &   (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
cc       write(8,'(10f8.4)') ((xnew(k,jj),k=1,ndim),jj=1,nelec)
c     endif
c     if(dabs(psin*(enew(1)-etrial)/etrial).gt.1.d-7) then
c        write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
c    &   (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
c     endif

c calculate probability for reverse transition

      fxnp=one
c Choose lower and upper values of r sampling
        rbot=rminn(i)*deltri
        rtop=rminn(i)*deltar
c Calculate magnitude of the velocity in the radial direction
        vnewr=zero
        do 110 k=1,ndim
  110     vnewr=vnewr+vnew(k,i)*rvminn(k,i)
        vnewr=vnewr/rminn(i)
        vnewr=max(vnewr,-2*znuc(iwctype(1)))

c Place x-axis along direction of angular change and
c Calculate the velocity in the phi direction
        vnewp=zero
        do 120 k=1,ndim
          xaxis(k)=vnew(k,i)-vnewr*rvminn(k,i)/rminn(i)
  120     vnewp=vnewp+xaxis(k)**2
        vnewp=dsqrt(vnewp)
        if(vnewp.lt.eps) then
          xaxis(1)=eps*(one-rvminn(1,i)**2)
          xaxis(2)=eps*(-rvminn(1,i)*rvminn(2,i))
          xaxis(3)=eps*(-rvminn(1,i)*rvminn(3,i))
          vnewp=eps*sqrt(one-rvminn(1,i)**2*(two-rminn(i)**2))
        endif
        do 130 k=1,ndim
  130     xaxis(k)=xaxis(k)/vnewp

c Temporary test of fbias
        vnewr=vnewr*fbias
        vnewp=vnewp*fbias

        if(vnewr.ge.zero) then
c         write(6,*) 'ln'
c Use linear approx for radial fn
          fxnp=fxnp*sqrt(rmino(i)/rtry)*abs(one+vnewr*(rmino(i)-rtry))

c Calculate the integral of T
          rzero=rminn(i)-one/vnewr
          if(rzero.lt.rbot .or. rzero.gt.rtop) then
            arean=area(rminn(i),rbot,rtop,vnewr)
           else
            arean=area(rminn(i),rbot,rzero,vnewr) +
     &      area(rminn(i),rzero,rtop,vnewr)
          endif

         else

c Use exponential approx for radial fn
            fxnp=fxnp*sqrt(rmino(i)/rtry)*exp(vnewr*(rmino(i)-rtry))
c   Calculate the integral of T
        zrbot=-vnewr*rbot
        zrtop=-vnewr*rtop
        zebot=dsqrt(zrbot**3)*dexp(-zrbot)
        zetop=dsqrt(zrtop**3)*dexp(-zrtop)
        g32bot=gammai(d3b2,zrbot,zebot,iflagb)
        g32top=gammai(d3b2,zrtop,zetop,iflagt)
        if(iflagb*iflagt.eq.1) then
          g32dif=g32top-g32bot
         else
          g32dif=g32top-g32bot+g3b2
        endif
        arean=g32dif
     &  /(dexp(vnewr*rminn(i))*(dsqrt(rminn(i)*(-vnewr)**3)))

c       arean=(gammai(d3b2,zrtop,zetop)-gammai(d3b2,zrbot,zebot))
c    &  /(dexp(vnewr*rminn(i))*(dsqrt(rminn(i)*(-vnewr)**3)))

        endif
c       write(6,'(''fs'',2f7.1,2f7.2)') voldr,vnewr,voldp,vnewp

c Truncate phi variation if it goes through zero
c Note it is OK to have a) theta or sin(theta)
c and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
c are consistent
c       term=dmin1(vnewp*rmino(i)*theta,one)
c       term=dmin1(vnewp*rmino(i)*sintht,one)
c       term=dmin1(vnewp*raver*sintht,one)
        term=vnewp*raver*sintht
c Determine cos(phi)
        cosphi=zero
        rnorm=zero
        do 160 k=1,ndim
c 160   cosphi=cosphi+(zaxis(k)-rvminn(k,i)/rminn(i))*xaxis(k)
c       cosphi=cosphi/(two*dsin(half*theta))
        term2=zaxis(k)-costht*rvminn(k,i)/rminn(i)
        rnorm=rnorm+term2*term2
  160   cosphi=cosphi+term2*xaxis(k)
        cosphi=cosphi/dsqrt(rnorm)
c       fxnp=fxnp*(one+term*cosphi)
        fxnp=fxnp*abs(one+term*cosphi)
        if(term.gt.one) arean=arean*(one+(term-one)**2/(two*term))


c p is the probability of accepting new move

      p=(rtry/rmino(i))**2*exp(psi2n(1)-psi2o(1))*
     &dabs((fxnp*areao)/(fxop*arean))
c     write(6,'(19d12.4)') rmino(i),rtry,psi2o, psi2n,fxop,fxnp,areao
c    &,arean,p,enew(1)-eold
      pp=(rtry/rmino(i))**2*exp(psi2n(1)-psi2o(1))*dabs(fxnp/fxop)
cr    write(6,'(''psi2n,psi2o,fxnp,fxop,psi2n/psi2o,fxnp/fxop,p'',9f12.6
cr   1)')psi2n,psi2o,fxnp,fxop,psi2n/psi2o,fxnp/fxop,p
c     write(6,'(''fxop,fxnp,p'',9d12.5)') fxop,fxnp,p

        if(ipr.ge.1) then
        write(6,'(''rminn,rvminn,vnew,vnewr'',9f10.4)')
     &  rminn(i),(rvminn(k,i),k=1,ndim),(vnew(k,i),k=1,ndim),vnewr
        write(6,'(''vnew='',9d12.4)') (vnew(k,i),k=1,ndim)
        write(6,'(''vnewr,vnewp='',9d12.4)') vnewr,vnewp
        write(6,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
        write(6,'(''rminn(i),rmax,rzero'',9f9.4)') rminn(i),rmax,rzero
        write(6,'(''rtry,costht,sintht,phitry,cosphi'',9f9.4)') rtry,
     &  costht,sintht,phitry,cosphi
        write(6,'(''fxop,fxnp,areao,arean,psi2n,psi2o,pp,p'',9f9.4)')
     &              fxop,fxnp,areao,arean,psi2n(1),psi2o(1),pp,p
        endif

      p=dmin1(one,p)
      q=one-p

c form expected values of e, pe, etc.

      esum1=          p*enew(1)+q*eold(1)
      esum(1)=esum(1)+p*enew(1)+q*eold(1)
      pesum=pesum+p*pen+q*peo
      peisum=peisum+p*pein+q*peio
      tpbsum=tpbsum+p*(enew(1)-pen)+q*(eold(1)-peo)
      tjfsum=tjfsum+p*tjfn+q*tjfo
      do 210 k=1,ndim
c       do 210 i=1,nelec
  210     r2sum=r2sum+p*xnew(k,i)**2+q*xold(k,i)**2
!      if(nefp.gt.0) then
!        call sample_efp(0,xold,eold(1),q)
!        call sample_efp(1,xnew,enew(1),p)
!      endif
      call acues1

c Calculate as a function of the distance to the nucleus
c 1) acceptance,  2) force-bias truncation probability,
c 3) kinetic energy and it's fluctuation
c The K.E. is not quite correct, since we should use p times new
c and q times old, and keep track of which bin the old was in
c The reason why I changed
c itryo=min(int(delri*rold)+1,NRAD)  to
c itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
c is that 2147483647 is the largest 32-bit integer and 1 more than that gives -2147483648.
      rold=dsqrt(xold(1,i)**2+xold(2,i)**2+xold(3,i)**2)
      rnew=dsqrt(xnew(1,i)**2+xnew(2,i)**2+xnew(3,i)**2)
c     itryo=min(nint(20*rold+1),80)
c     itryn=min(nint(20*rnew+1),80)
      itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
      itryn=int(min(delri*rnew+1,dfloat(NRAD))+eps)
      try(itryo)=try(itryo)+1
      suc(itryo)=suc(itryo)+p
      if(try(itryo).lt.0.) write(6,'(''itryo,try'',i5,d13.5)')itryo,
     &try(itryo)
      if(suc(itryo).lt.0.) write(6,'(''itryo,suc'',i5,d13.5)')itryo,
     &suc(itryo)
      if(voldp*raver*sintht.gt.one) trunfb(itryo)=trunfb(itryo)+1
      if(i.le.nup) then
        rprobup(itryo)=rprobup(itryo)+q
        rprobup(itryn)=rprobup(itryn)+p
       else
        rprobdn(itryo)=rprobdn(itryo)+q
        rprobdn(itryn)=rprobdn(itryn)+p
      endif
      rprob(itryo)=rprob(itryo)+q
      rprob(itryn)=rprob(itryn)+p
      ekin(itryo)=ekin(itryo)+q*ekineo(i)
      ekin2(itryo)=ekin2(itryo)+q*ekineo(i)**2
      ekin(itryn)=ekin(itryn)+p*ekinen(i)
      ekin2(itryn)=ekin2(itryn)+p*ekinen(i)**2

c     eksum=zero
c     do 220 ii=1,nelec
c 220   eksum=eksum+ekine(ii)
c     write(6,'(''ke='',9d13.5)') eksum,enew(1)-pen


c accept new move with probability p

      if(rannyu(0).gt.p) then
        do 230 k=1,ndim
  230     xnew(k,i)=xold(k,i)
c Warning: have to set igeometrical
        if(igeometrical.eq.0) call distancese_restore(i)
       else

c move is accepted so update positions etc.
        rmino(i)=rminn(i)
        do 240 k=1,ndim
          xold(k,i)=xnew(k,i)
          rvmino(k,i)=rvminn(k,i)
          do 240 ii=1,nelec
  240       vold(k,ii)=vnew(k,ii)
        accsum=accsum+one
        eold(1)=enew(1)
        peo=pen
        peio=pein
        tjfo=tjfn
        psi2o(1)=psi2n(1)
        psido=psidn
        psijo=psijn
        ekineo(i)=ekinen(i)

!        if(nefp.gt.0) call efpsav
      endif

  300 continue

      call object_modified_by_index (xold_index)  !JT

c Warning: does not have correlated sampling for forces yet.
      do 380 ifr=1,nforce
        if(ifr.eq.1) then
          esum1s(ifr)=eold(ifr)
         else
          wstro=exp(psi2o(ifr)-psi2o(1))
          wsum1s(ifr)=wstro
          esum1s(ifr)=eold(ifr)*wstro
        endif
  380 continue

      call acusig

      return
      end

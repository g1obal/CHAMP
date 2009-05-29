      subroutine readps
c Written by Claudia Filippi
c Read pseudopotential in Fahy format
      use all_tools_mod
      use atom_mod
      use pseudo_mod
      use qua_mod
      use pseudo_fahy_mod
      implicit real*8(a-h,o-z)

      character*20 filename,atomtyp

c nquad = number of quadrature points
c nlang = number of non-local potentials
c rcmax = cutoff radius for non-local potential
c npotl = number of mesh point for local potential
c dradl = spacing of uniform mesh for local potential

      call alloc ('potl', potl, MPS_GRID, nctype)
      call alloc ('ptnlc', ptnlc, MPS_GRID, nctype, MPS_L)
      call alloc ('ptnlc', ptnlc, MPS_GRID, nctype, MPS_L)
      call alloc ('dradl', dradl, nctype)
      call alloc ('drad', drad, nctype)
      call alloc ('rcmax', rcmax, nctype)
      call alloc ('npotl', npotl, nctype)
      call alloc ('nlrad', nlrad, nctype)

      do 20 ic=1,nctype

      if(ic.lt.10) then
        write(atomtyp,'(i1)') ic
       elseif(ic.lt.100) then
        write(atomtyp,'(i2)') ic
       else
        stop 'readps_tm, problem atomtyp'
      endif

      filename='ps.data.'//atomtyp(1:index(atomtyp,' ')-1)
      open(3,file=filename,status='old',form='formatted')

      read(3,*) nquad
      write(6,'(''quadrature points'',i4)') nquad

      read(3,*) nlang,rcmax(ic)
c If the local pseudopot component is not set in input, set it here
      if(lpotp1(ic).lt.0) then
        lpotp1(ic)=nlang+1
        write(6,'(''local pseudopot component is'',i3)') lpotp1(ic)
      endif

c local potential
      read(3,*)
      read(3,*) npotl(ic),nzion,dradl(ic)
      if(npotl(ic).gt.MPS_GRID) stop 'npotl gt MPS_GRID'
      if(nzion.ne.int(znuc(iwctype(ic)))) stop 'nzion ne znuc'

      read(3,*) (potl(i,ic),i=1,npotl(ic))

      do 5 i=1,npotl(ic)
  5     write(33,*) (i-1)*dradl(ic),potl(i,ic)

c non-local potential
      read(3,*)
      read(3,*) nlrad(ic),drad(ic)
      if(nlrad(ic).gt.MPS_GRID) stop 'nlrad gt MPS_GRID'

      if(drad(ic)*(nlrad(ic)-1).le.rcmax(ic)) then
        write(6,'(''non-local table max radius = '',
     &  f10.5,'' too small for cut-off = '',f10.5)')
     &  drad(ic)*(nlrad(ic)-1),rcmax(ic)
        stop
      endif

      do 10 l=1,nlang
        read(3,*)
        read(3,*) (ptnlc(i,ic,l),i=1,nlrad(ic))
      do 10 i=1,nlrad(ic)
 10     write(34,*) (i-1)*drad(ic),ptnlc(i,ic,l)

      close(3)
 20   continue


c     call gesqua (nquad,xq0,yq0,zq0,wq)
c     call gesqua (nquad,xq0,yq0,zq0,wq)
c     call gesqua (nquad,xq,yq,zq,wq)

c     write(6,'(''quadrature points'')')
c     do 30 i=1,nquad
c30     write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      return
      end
c-----------------------------------------------------------------------
      subroutine getvps_fahy(rad,iel)
c Written by Claudia Filippi
c compute Fahy-pseudopotential for electron iel
      use atom_mod
      use const_mod
      use pseudo_mod
      use pseudo_fahy_mod
      implicit real*8(a-h,o-z)

      dimension rad(nelec,ncent)

      do 10 ic=1,ncent
        ict=iwctype(ic)
        r=rad(iel,ic)
c local potential
        if(r.lt.(npotl(ict)-1)*dradl(ict)) then
          ri=r/dradl(ict)
          ir=int(ri)
          ri=ri-dfloat(ir)
          ir=ir+1
          vps(iel,ic,lpotp1(ict))=potl(ir+1,ict)*ri+(1.d0-ri)*potl(ir,ict)
         else
          vps(iel,ic,lpotp1(ict))=-znuc(ict)/r
        endif
c non-local pseudopotential
        do 10 l=1,npotd(ict)
          if(l.ne.lpotp1(ict)) then
            if(r.lt.rcmax(ict)) then
              ri=r/drad(ict)
              ir=int(ri)
              ri=ri-dfloat(ir)
              ir=ir+1
              vps(iel,ic,l)=ptnlc(ir+1,ict,l)*ri+(1.d0-ri)*ptnlc(ir,ict,l)
             else
              vps(iel,ic,l)=0.0d0
            endif
          endif
   10 continue

      return
      end

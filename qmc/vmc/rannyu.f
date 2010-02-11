c NYU linear congruential random number generator.
c R_{n+1} = M*R_n mod 2^48
c Since there is no additive constant the period is less than 2^48 = 2.81474976710656 x 10^14
c Uses 48 bits, rather than the usual 32 bits since it was developed in 1960s
c for CDC computers which had 48-bit integers.
c Seed and multiplier are 48-bit integers stored as 4 twelve-bit integers.
c Multiplier is 11^13 = 34522712143931
c M is octal so
c M(1-4)= 502,1521,4071,2107 since 502*8^12+1521*8^8+4071*8^4+2107 = 34522712143931

      subroutine setrn(iseed)
c Set seed and make sure it is an odd integer.
      implicit real*8(a-h,o-z)
      common /rnyucm/ m(4),l(4)
      integer iseed(4)
         do 10 i=1,4
         l(i)=iseed(i)
   10    continue
      l(4)=2*(l(4)/2)+1
      return
      end
c-----------------------------------------------------------------------

      block data
      implicit real*8(a-h,o-z)
      common /rnyucm/ m(4),l(4)
      data m / 502,1521,4071,2107/
      data l /   0,   0,   0,   1/
      end
c-----------------------------------------------------------------------

c     function rannyu(idum)
c     implicit real*8(a-h,o-z)
c     parameter(two=2.d0)
c     common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4
c     i1=l1*m4+l2*m3+l3*m2+l4*m1
c     i2=l2*m4+l3*m3+l4*m2
c     i3=l3*m4+l4*m3
c     i4=l4*m4
c     l4=mod(i4,2**12)
c     i3=i3+i4/2**12
c     l3=mod(i3,2**12)
c     i2=i2+i3/2**12
c     l2=mod(i2,2**12)
c     l1=mod(i1+i2/2**12,2**12)
c     rannyu=two**(-12)*(dfloat(l1)+
c    &       two**(-12)*(dfloat(l2)+
c    &       two**(-12)*(dfloat(l3)+
c    &       two**(-12)*(dfloat(l4)))))
c     return
c     end

      function rannyu(idum)
      implicit real*8(a-h,o-z)
c     On bat it is more efficient to not precompute 2**12 and 2**-12
c     Therefore use original verion of code instead of this one.
      parameter(itwo12=4096,two12i=2.44140625d-4)
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4
      i1=l1*m4+l2*m3+l3*m2+l4*m1
      i2=l2*m4+l3*m3+l4*m2
      i3=l3*m4+l4*m3
      i4=l4*m4
      l4=mod(i4,itwo12)
      i3=i3+i4/itwo12
      l3=mod(i3,itwo12)
      i2=i2+i3/itwo12
      l2=mod(i2,itwo12)
      l1=mod(i1+i2/itwo12,itwo12)
      rannyu=two12i*(dfloat(l1)+
     &       two12i*(dfloat(l2)+
     &       two12i*(dfloat(l3)+
     &       two12i*(dfloat(l4)))))
      return
      end
c-----------------------------------------------------------------------

      subroutine savern(iseed)
c Return current seed.
      implicit real*8(a-h,o-z)
      common /rnyucm/ m(4),l(4)
      integer iseed(4)
         do 10 i=1,4
         iseed(i)=l(i)
   10    continue
      return
      end

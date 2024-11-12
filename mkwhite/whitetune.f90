  program whitetune
  parameter(nw=5, nsqr=nw*nw)
  integer nx,ny,nz,mx,local_nz
  parameter(nid=32)
  parameter(nx=2048,ny=nx,nz=nx,local_nz = nx/nid)
  parameter(mx = 2*(nx/2+1))

  external gasdev
  real gasdev,xxx
  integer xksq
  real rl,im,amp,xisq

  integer iseed,myid,nid,iiseed
  integer nsize
  real invsqrt2

  real pk(nsqr)
  integer npk(nsqr)
  real, dimension(:), allocatable:: dist
  integer, dimension(:), allocatable:: indx
! integer indx(mx*ny*local_nz)
! real dist(mx*ny*local_nz)


  invsqrt2 = 1./sqrt(2.)
  nsize = mx*ny*local_nz



! do kk = 1, 100000000
!    do kk =  1,800000000
   do kk =  800000001,2000000000
     do i = 1, nsqr
        pk(i) = 0
        npk(i) = 0
     enddo
     iseed = -kk
     xxx = ran3(iseed)
     do i = 0, nw
     do j = 0, nw
     do k = 0, nw
        rl = gasdev(iseed)
        im = gasdev(iseed)
        xksq =  (i*i+j*j+k*k)
        if(xksq.gt.0 .and. xksq.le.nsqr) then
           amp = (rl*rl + im*im)/2.
           npk(xksq) = npk(xksq) + 1
           pk(xksq) = pk(xksq) + amp
        endif
     enddo
     enddo
     enddo
     xisq = 0
     icount = 0
     do i = 1, nsqr
        if(npk(i).gt.0) then
        pk(i) = pk(i)/npk(i)
        xisq = xisq + (pk(i)-1.)*(pk(i)-1.)
        icount = icount + 1
        endif
     enddo
     xisq = xisq/icount
     xisq = sqrt(xisq)
     if(xisq .lt. 0.26) then
        write(*,'(i10,11(1x,g12.5))') kk,xisq,pk(1),pk(2),pk(3),pk(4),pk(5),pk(6),pk(7),pk(8),pk(9),pk(10)
     endif

  enddo


  stop
  end
!c--

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
!CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
!C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.


      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.


  program whitetune
  use grafic_types


  integer, parameter :: RECORD_FLAG = 4

  parameter(nw=5, nsqr=nw*nw)
  integer nx,ny,nz,mx,local_nz
! parameter(nid=250)
! parameter(nx=2048,ny=nx,nz=nx)
  parameter(nid=8)
  parameter(nx=512,ny=nx,nz=nx)
  parameter(mx = 2*(nx/2+1))

  external gasdev
  real gasdev,xxx,work1
  integer xksq,ierror
  real rl,im,amp,xisq

  integer iseed,myid,nid,iiseed
  integer nsize,sendsize
  integer(RECORD_FLAG) :: taille_tampon
  real(sp), dimension(:), allocatable:: den
  character(len=128) :: filename
  real invsqrt2
  real*8 rinvcube,rngc

  integer(i8b) :: plan, iplan
  integer(i8b) :: offset,length
  
  integer total_local_size
  integer local_z_start, local_ny,local_y_start
  real denpad(2,ny,nz), tampon(nx*ny)
  integer status(MPI_STATUS_SIZE)

  real pk(nsqr), npk(nsqr)


  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)

  call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nx,ny,nz, &
       fftw_real_to_complex, fftw_estimate+fftw_in_place)
  call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nx,ny,nz, &
       fftw_complex_to_real, fftw_estimate+fftw_in_place)
  call rfftwnd_f77_mpi_local_sizes(plan,local_nz,local_z_start,local_ny, &
       local_y_start,total_local_size)

  do i = 1, nsqr
     pk(i) = 0
     npk(i) = 0
  enddo

  invsqrt2 = 1./sqrt(2.)
  nsize = mx*ny*local_nz

  if(nsize .gt.0) then
  allocate(den(nsize))
  else
  allocate(den(nx*ny))
  endif

! iseed = -152090039
!  iseed = -887843645
! iseed = -807601913
! iseed = -1807601913
  iseed = -56
  iiseed = iseed

  iseed = iseed + myid
  xxx = ran3(iseed)


  rinvcube = 1.d0/dble(nx)/dble(ny)/dble(nz)
! conversion factor from k(=2pi*u/Lbox) to u
  rngc = dsqrt(dble(nx)*dble(ny)*dble(nz))
  !  renormalization is postponed inside mpi4c.f90
! rngc = 1.d0

  if(myid .ne.0) then
      do i = 0, nx/2
      do j = 0, ny-1
      do k = 0, local_nz-1
         rl = gasdev(iseed)*invsqrt2*rngc
         im = gasdev(iseed)*invsqrt2*rngc
         den(2*i+1+mx*(j+ny*k)) = rl
         den(2*i+2+mx*(j+ny*k)) = im
      enddo
      enddo
      enddo
  else if(myid.eq.0) then
!    do i = 0, nw
!    do j = 0, nw
!    do k = 0, nw
!       rl = gasdev(iseed)*invsqrt2*rngc
!       im = gasdev(iseed)*invsqrt2*rngc
!       den(2*i+1+mx*(j+ny*k)) = rl
!       den(2*i+2+mx*(j+ny*k)) = im
!       ii = i*i + j*j + k*k
!       if(ii.le.nsqr) then
!          pk(ii) = pk(ii) + (rl*rl + im*im)/2.
!          npk(ii) = npk(ii) + 1
!       endif
!    enddo
!    enddo
!    enddo
     do i = 0, nx/2
     do j = 0, ny-1
     do k = 0, local_nz-1
!       if(i.le.nw .and. j.le.nw.and.k .le.nw) goto 345
        rl = gasdev(iseed)*invsqrt2*rngc
        im = gasdev(iseed)*invsqrt2*rngc
        den(2*i+1+mx*(j+ny*k)) = rl
        den(2*i+2+mx*(j+ny*k)) = im
!345     continue
      enddo
      enddo
      enddo
     den(1) = 0
     den(2) = 0
  else
     print *, 'Error. ', myid
  endif

  if(myid.eq.0) then
     print *, iiseed
     do i = 1, nsqr
        if(npk(i).gt.0) then
           pk(i) = pk(i)/npk(i)
           print *, i, sqrt(real(i)),pk(i),npk(i)
        endif
     enddo
  endif

! call MPI_Finalize(error)
! stop



  do k = 1, local_nz
  do j = 1, ny
  do i = 1, 2
          denpad(i,j,local_z_start+k) = den(i+mx*(j-1+ny*(k-1)))
  enddo
  enddo
  enddo
  nghz = nz/2
  nghy =ny/2

  if(myid .eq. 0) then
         do i = 1, nid-1
           call MPI_RECV(mlocalstartz,1,MPI_INTEGER,i,0,MPI_COMM_WORLD, &
     &              status,ierror)
           call MPI_RECV(mlocalnz,1,MPI_INTEGER,i,0,MPI_COMM_WORLD, &
     &              status,ierror)
           call MPI_RECV(denpad(1,1,mlocalstartz+1),2*ny*mlocalnz, &
     &              MPI_REAL,i,0,MPI_COMM_WORLD,status,ierror)
         enddo

  else
     call MPI_SEND(local_z_start,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,   &
     &             ierror)
         call MPI_SEND(local_nz,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,   &
     &             ierror)
         call MPI_SEND(denpad(1,1,local_z_start+1),2*ny*local_nz,   &
     &     MPI_REAL,0,0,MPI_COMM_WORLD,ierror)
   endif
   if(myid .eq. 0) then
         do 12 k = 2, nghz
            denpad(1,1,nz+2-k) = denpad(1,1,k)
            denpad(2,1,nz+2-k) =-denpad(2,1,k)
            do 12 j=2,nghy
            denpad(1,ny+2-j,nz+2-k) = denpad(1,j,k)
12          denpad(2,ny+2-j,nz+2-k) =-denpad(2,j,k)

         do 13 j=2,nghy
            denpad(1,ny+2-j,1) = denpad(1,j,1)
13          denpad(2,ny+2-j,1) =-denpad(2,j,1)

         do 14 k=2,nghz
         do 14 j=nghy+2,ny
            denpad(1,ny+2-j,nz+2-k) = denpad(1,j,k)
14          denpad(2,ny+2-j,nz+2-k) =-denpad(2,j,k)
   endif
   sendsize = 2*ny*nz
   call MPI_BCAST(denpad,sendsize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
   if(myid.eq.0) print *,'just after broadcasting the conjugate'
   do k = 1, local_nz
   do j = 1, ny
   do i = 1, 2
          den(i+mx*(j-1+ny*(k-1))) = denpad(i,j,local_z_start+k)
   enddo
   enddo
   enddo
   call rfftwnd_f77_mpi(iplan,1,den,work1,use_work, FFTW_NORMAL_ORDER)

   filename = 'white.out'
   if(myid.eq.0) then
      open(1,file=trim(filename),form='unformatted')
      write(1) nx,ny,nz,iiseed
!     offset = 4*4+2*RECORD_FLAG + nz*int(nx*ny*sp + 2*RECORD_FLAG,8)-4
!     call fseek(0, offset,SEEK_SET,ierror)
!     write(1) iiseed
      close(1)
   endif

   offset = 4*4 + 2*RECORD_FLAG
   taille_tampon = nx*ny*sp
   offset = offset + local_z_start*int(nx*ny*sp + 2*RECORD_FLAG,8)

   do k = 1, total_local_size
      den(k) = den(k)*rinvcube
   enddo



   if(0) then
      do ii = 0, nid-1
         if(myid.eq.ii) then
         open(1,file=trim(filename), form='unformatted', access='append')
         do k = 1, local_nz
            do j = 1, ny
               do i = 1, nx
                   tampon(i+nx*(j-1)) = den(i+mx*(j-1+ny*(k-1)))
               enddo
            enddo
            write(1) (tampon(i),i=1, nx*ny)
         enddo
         close(1)
         endif
         call mpi_barrier(MPI_COMM_WORLD,ierror)
      enddo
   else 
      do k = 0, nid-1
         length = RECORD_FLAG
         offset = 4*4 + 2*RECORD_FLAG
         offset = offset + local_z_start*int(nx*ny*sp + 2*RECORD_FLAG,8)
         if(myid.eq.ii) then
            do j = 1, local_nz
               call f77_parallel_write(trim(filename),len_trim(filename),length, &
     &             offset, taille_tampon)
               offset = offset + RECORD_FLAG + nx*ny*sp
               call f77_parallel_write(trim(filename),len_trim(filename),length, &
     &             offset, taille_tampon)
               offset = offset + RECORD_FLAG
            enddo
         endif
         call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo

      do k = 1, local_nz
         do j = 1, ny
            do i = 1, nx
                tampon(i+nx*(j-1)) = den(i+mx*(j-1+ny*(k-1)))
            enddo
         enddo
         length = RECORD_FLAG
         call f77_parallel_write(trim(filename),len_trim(filename),length, &
      &          offset, taille_tampon)
         offset = offset + RECORD_FLAG
         length = nx*ny*sp
         call f77_parallel_write(trim(filename),len_trim(filename),length, &
      &          offset, tampon)
         offset = offset + length
         length = RECORD_FLAG
         call f77_parallel_write(trim(filename),len_trim(filename), length &
               & ,offset,taille_tampon)
         offset = offset + RECORD_FLAG
   enddo

   endif
!  deallocate(tampon)

   call mpi_barrier(MPI_COMM_WORLD, ierror)


   call mpi_finalize(ierror)
   if(myid.eq.0) print *, 'well ended!!!!'

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


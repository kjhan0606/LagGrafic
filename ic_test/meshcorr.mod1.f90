!---  Get the particle position data and calculates the power spectrum
!---  by Juhan Kim @ KIAS in 20/06/2024

! ******************************************************************************
! *                                                                            *
! *                              meshcorr                                      *
! *                                                                            *
! ******************************************************************************
!--- meshcorr: main program
      program meshcorr
      real rminmass
      integer redflag
      integer readflag
      integer iorder
      real threshold(1)
      character*80 args
      data threshold /0./

      INCLUDE "meshcorrps.inc"

      call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierror)
      call getarg(1,infile)
      call getarg(2,outfile)
      call getarg(3,args)
      read(args,*) imax
      jmax = imax
      kmax = imax
      mmax = imax


      outfile2 = 'mode.'//outfile



      call rfftw3d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,imax,jmax, &
     &       kmax, &
     &       FFTW_REAL_TO_COMPLEX,FFTW_MEASURE+FFTW_IN_PLACE)
      call rfftw3d_f77_mpi_create_plan(iplan,MPI_COMM_WORLD,imax,jmax, &
     &       kmax, &
     &       FFTW_COMPLEX_TO_REAL,FFTW_MEASURE+FFTW_IN_PLACE)
      call rfftwnd_f77_mpi_local_sizes(plan,local_nz,local_z_start, &
     &       local_ny_after_transpose,local_y_start_after_transpose, &
     &       total_local_size)

!--  flag is for reading or not
      flag = 0
      print *,'setting up fftw stuffs'

      do 555 ii = 1, 1
      initflag = 0
      ni = imax
      nj = jmax
      nk = kmax
      smoothflag = 0
      call readinandpkmeasure(rminmass,redflag)
      smoothflag = 1
       close(9)
555   continue 
      call  MPI_FINALIZE(ierror)
      stop
      end
!--- draw: find threshold density corresponding to desired volume (or mass) 
!--- fraction, and define density contour



      subroutine readinandpkmeasure(rminmass,redflag)

      INCLUDE "meshcorrps.inc"

      real rminmass
      integer redflag,iseed

!--- rho is really a local array, but declaring it common stops the Ultrix
!--- fort compiler from making a huge executable.  Feel free to declare it
!--- local if that fits your compiler/operating system better
!     common /classi/ rho(2*(imax/2+1),jmax,kmax)
!     real*4 den(2*(imax/2+1),jmax,kmax/nrank)
      real, allocatable :: den(:,:,:)
      integer i,j,k
      integer*8 high
      integer*8 thigh
      real*8 summm,tsummm
      real tmass
      double precision ls1,ls0,lmt,lvt
      double precision s1,s0,mt,vt
      real gasdev
      external gasdev
      integer icount

      allocate(den(2*(imax/2+1),jmax,kmax/nrank))

      if (flag.eq.0) then
      nsize = 2*(imax/2+1)*jmax*(kmax/nrank)
      itag = 10
      if(myrank.eq. nrank-1) then
         open(1,file=infile,form='unformatted')
         read(1) ni,nj,nk,iseed
         print *, ni,nj,nk,imax,jmax,kmax
         nstep = nk/nrank
         icount = 0
         do jj = 0, nrank-1
           do k =1, kmax/nrank
           icount = icount + 1
           read(1) ((den(i,j,k),i=1, imax),j=1,jmax)
           print *, 'icount = ',icount, den(1,1,k)
           enddo
           if(jj .ne. nrank-1) then
              call MPI_Send(den,nsize,MPI_REAL,jj,itag,&
     &             MPI_COMM_WORLD,ierror)
              print *,'sent mesh data to ', jj
           endif
         enddo
         close(1)
      else
         call MPI_Recv(den,nsize,MPI_REAL,nrank-1,itag,MPI_COMM_WORLD, &
     &            istatus,ierror)
      endif
      summm = 0
      do k = 1, kmax/nrank
      do j = 1, jmax
      do i = 1, imax
      summm = summm + den(i,j,k)
      enddo
      enddo
      enddo
      print *,myrank,'local sum=', summm
      call MPI_REDUCE(summm,tsummm,1,MPI_REAL8,MPI_SUM, &
     &        0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(tsummm,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      tsummm = tsummm/(dble(kmax)*dble(jmax)*dble(imax))
      if(myrank .eq. 0)print *,myrank,'av(dTb) = ', tsummm

      redshift = 15.d0
      Tnow = 2.725d0
      Tcmb = Tnow*(1.d0+redshift)*1000.d0
!     if(myrank .eq. 0)print *,'Tcmb(mK) = ', Tcmb

      do k = 1, kmax/nrank
      do j = 1, jmax
      do i = 1, imax
      den(i,j,k) = (den(i,j,k)-tsummm)
      enddo
      enddo
      enddo
      endif
      call pkmeasure(den)

      deallocate(den)
      return
      end

!--- sum: calculate genus of boundary surface by summing angle deficits
!
!---    At each vertex, sum evaluates the contribution to the angle deficit
!--- of the boundary surface using the look-up tables d(1,0:255) and 
!--- d(2,0:255), depending on whether the geometry of the vertex is type I
!--- or type II (see GMD).  The angle deficit at a vertex is a function of the
!--- density values (high or low) of the eight cells surrounding the vertex.
!--- When written out in binary notation, the subscript config of d(type,config)
!--- describes the configuration of the cells surrounding the vertex, 
!--- 1 indicating high density and 0 low density.  Angle deficits are 
!--- evaluated as if the cube rows were offset from each other in a regular 
!--- pattern, in order to eliminate ambiguities in the location of the 
!--- boundary surface (see GMD fig. 10).  However, the results are not 
!--- sensitive to, or systematically dependent on, this prescription.
!--- Only vertices that lie entirely within the survey region are counted.

!--- This version, vsum, has a simplified loop structure that makes it
!--- faster on a vector machine.


!--- denread: read density from file "rho.smooth" into array rho(i,j,k)
!---  Calculates the smooth density field from particle realization
!---  by using CIC scheme


      subroutine pkmeasure(den)
      parameter(nw=10, nsqr= nw*nw)
      integer pflag
      INCLUDE "meshcorrps.inc"
      real*4 den(2*(imax/2+1),jmax,kmax/nrank)
      real lambda 
      character*80 iname,tmpchar
      real*8 rnorm,tnorm
      real zheight,zstart
      real work
      real rincube
!     integer*2 tabp0(0:mmax+1)
!     integer*2 tabpz(0:mmax+1)
!     real pk(mmax+1),npk(mmax+1),cor(mmax+1),pksig(mmax+1)
!     real tpk(mmax+1),tnpk(mmax+1),tcor(mmax+1),tpksig(mmax+1)
!     real pkconv(mmax+1),tpkconv(mmax+1)
      integer*2, allocatable :: tabp0(:), tabpz(:)
      real, allocatable, dimension(:) :: pk, npk, cor, pksig
      real, allocatable, dimension(:) :: tpk, tnpk, tcor, tpksig
      real, allocatable, dimension(:) :: pkconv, tpkconv

      real bpk(nsqr),nbpk(nsqr)
      real tbpk(nsqr),tnbpk(nsqr)
      integer*8 iii


      allocate(tabp0(0:mmax+1))
      allocate(tabpz(0:mmax+1))
      allocate(pk(mmax+1), npk(mmax+1),cor(mmax+1),pksig(mmax+1))
      allocate(tpk(mmax+1),tnpk(mmax+1),tcor(mmax+1),tpksig(mmax+1))
      allocate(pkconv(mmax+1),tpkconv(mmax+1))


      do i = 1, nsqr
         bpk(i) = 0
         nbpk(i) = 0
      enddo

      pi = 4.*atan(1.0)
      npix=imax
      xside = float(npix)
      n2 = npix/2



      rnorm = 0
      do k = 1, kmax/nrank
      do j = 1, jmax
      do i = 1, imax
         rnorm =  rnorm + den(i,j,k)
      enddo
      enddo
      enddo
      call MPI_REDUCE(rnorm,tnorm,1,MPI_REAL8,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      call MPI_Bcast(tnorm,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      rnorm = tnorm/ (dble(imax)*dble(jmax)*dble(kmax))
      print *,'rnorm= ', rnorm



!--- forward transform
!     call src3ft(den,npix,npix,npix,npix+2,npix,1,ier)
!     if(ier.ne.0)print *,'*** FFT error *** : ier=',ier
!cc   Calculate a real-to-complex FFT in transposed order
      call rfftwnd_f77_mpi(plan,1,den(1,1,1),work,use_work, &
     &               FFTW_TRANSPOSED_ORDER)

! ---- This is inserted for the correlation measuring
      nx = imax
      ny = jmax
      nz = kmax
      rngc = real(nx)*real(ny)*real(nz)
      rngm6 = 1./rngc/rngc
      xyz = real(nx)/real(nz)
      pi = 3.141592654
      ng = nx
      rng = ng
      ngx = nx
      ngy = ny
      ngz = nz
      rngx = ngx
      rngy = ngy
      rngz = ngz
      xyz = rng/rngz
      xyy = rng/rngy
      ngxh = ngx/2
      ngh = ngxh
      ngyh = ngy/2
      ngzh = ngz/2
      do 90 i=1,ngx+1
      pk(i) = 0.
      pkconv(i) = 0.
      pksig(i) = 0.
      npk(i) = 0
      tpk(i) = 0.
      tpkconv(i) = 0.
      tnpk(i) = 0
   90 cor(i) = 0.


      do j = 1, local_ny_after_transpose, +1
      wj = xyy*(j - 1 + local_y_start_after_transpose)
      if(wj.gt.ngh) wj = wj - ngx
!     if(wj.eq.0) then
!     tscdeconvy = 1
!     else 
!     tscdeconvy = (sin(pi*wj/ng)/(pi*wj/ng))**3
!     endif
      do k = 1, ngz, +1
      wk = xyz*(k - 1)
      if(wk.gt.ngh) wk = wk - ngx
!     if(wk.eq.0) then
!     tscdeconvz = 1
!     else
!     tscdeconvz = (sin(pi*wk/ng)/(pi*wk/ng))**3
!     endif
      do i = 0, ngh, +1
      wave = sqrt(i*i+wj*wj+wk*wk)
!---  density array in k-space is transposed order: (i,j,k)-->(i,k,j)
      delksq = (den(2*i+1,k,j)**2+den(2*i+2,k,j)**2)*rngm6
      kbin1 = nint(wave)+1
      wt1 = wave-kbin1+1
      kbin2 = kbin1+sign(1.,wt1)
      wt1 = abs(wt1)
!     if(i.eq.0) then
!     tscdeconvx = 1
!     else
!     tscdeconvx = (sin(pi*real(i)/ng)/(pi*real(i)/ng))**3
!     endif
!     tscdeconv = tscdeconvx * tscdeconvy * tscdeconvz

      pk(kbin1) = pk(kbin1)+(1.-wt1)*delksq
!     pkconv(kbin1) = pkconv(kbin1)+(1.-wt1)*delksq/tscdeconv/tscdeconv
      pkconv(kbin1) = pkconv(kbin1)+(1.-wt1)*delksq
      pksig(kbin1) = pksig(kbin1)+(1.-wt1)*delksq**2
      npk(kbin1) = npk(kbin1)+1.-wt1
      pk(kbin2) = pk(kbin2)+wt1*delksq
!     pkconv(kbin2) = pkconv(kbin2)+wt1*delksq/tscdeconv/tscdeconv
      pkconv(kbin2) = pkconv(kbin2)+wt1*delksq
      pksig(kbin2) = pksig(kbin2)+wt1*delksq**2
      npk(kbin2) = npk(kbin2)+wt1
      
      iii = i*i+wj*wj+wk*wk

      if(iii .gt.0 .and. iii .le. nsqr) then
         bpk(iii) = bpk(iii) + delksq
         nbpk(iii) = nbpk(iii) + 1
      endif

      enddo
      enddo
      enddo

      print *, xyy,xyz,kbin1,kbin2,wt1,delksq

!--   Gathering pk in all processors and summing
      call MPI_REDUCE(pk,tpk,nx+1,MPI_REAL,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(bpk,tbpk,nsqr,MPI_REAL,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(pkconv,tpkconv,nx+1,MPI_REAL,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(pksig,tpksig,nx+1,MPI_REAL,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(npk,tnpk,nx+1,MPI_REAL,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(nbpk,tnbpk,nsqr,MPI_REAL,MPI_SUM,0, &
     &           MPI_COMM_WORLD,ierror)
      if(myrank.eq.0) then
         do i = 1, nx+1
!           print *, i, tpk(i)/tnpk(i)
            pk(i) = tpk(i)
            pkconv(i) = tpkconv(i)
            pksig(i) = tpksig(i)
            npk(i) = tnpk(i)
         enddo
      endif
      call MPI_BCAST(pk,nx+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(pkconv,nx+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(pksig,nx+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(npk,nx+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)

      rth = 8.* rngx/size
      a = 16.

      if(myrank.eq.0) then
         std = 0.
         fact = 2.*pi*rth/rngx
         do 110 i=2,ngx+1
         if(npk(i).gt.0)then
            pk(i)=pk(i)/npk(i)
            pkconv(i)=pkconv(i)/npk(i)
            pksig(i)=sqrt(pksig(i)/npk(i)-pk(i)**2)
         else
            pk(i) = 0.
            pkconv(i) = 0.
            pksig(i) = 0.
         end if

         thk = fact*(i-1.)
  110    std = std+pk(i)*(i-1.)**2*(sin(thk)-thk*cos(thk))**2/thk**6
         std = sqrt(4.*pi*std*9.)
         do i=2,ngx+1
            seper = 2.*pi*(i-1.)/rngx
            do k=2,ngx+1
               cor(i) = cor(i)+sin((k-1.)*seper)*(k-1.)*pk(k)
            enddo
            cor(i) = 2.*rngx*cor(i)/(i-1.)
         enddo
         if(cor(4).gt.0..and.cor(3).gt.0.)then
            corn2 = (alog(cor(4))-alog(cor(3)))/(alog(4.)-alog(3.))
         endif
         corn1 = corn2
         open(7,file=outfile)
         write(7,*)a,std,corn1
         do i = 1, ngh
            write(7,444) i,i-1,pkconv(i),cor(i),pksig(i),npk(i)
         enddo
         close(7)
         open(8,file=outfile2)
         do i = 1, nsqr
            if(tnbpk(i) .gt.0) then
              tbpk(i) = tbpk(i)/tnbpk(i)
              write(8,*) sqrt(real(i)), tbpk(i), tnbpk(i)
            endif
         enddo
         close(8)
      endif

444   format(2(i5,1x),4(g14.7,1x))
      

!---  stop
      call MPI_Finalize(ierror)
      stop
   
      end

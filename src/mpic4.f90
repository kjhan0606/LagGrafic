!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getPhi1Phi2(irand,iseed,itide,m1s,m2s,m3s,m1off, &
     &    m2off,m3off,hanning,filename,astart,pk,dx,xoff,phi1,phi2,fm, &
     &    plan,iplan,local_nz,local_z_start,total_local_size, &
     &    headt,headc,small_kfile_name,ipad, anorml)
  !  Generates phi1 and phi2 for the 2LPT.

  use grafic_types
  use grafic_io
  use transform
  use normalize
  use random
  use paste

  implicit none
  include 'grafic1.inc'
  integer n12,n22,n32,n23,n2p1,npow
  parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
  parameter (n2p1=2*(np1/2+1))
! parameter (npow=30720)
  parameter (npow=4096)
  integer irand,iseed,itide,m1s,m2s,m3s,m1off,m2off,m3off
  logical hanning
  character(len=128) :: filename, small_kfile_name
  integer :: ipad
  real :: anorml
  type(taille) :: headt
  type(cosmo) :: headc
  real astart,pk,dx,xoff,fm
  external pk
  integer(i8b) :: plan, iplan
  integer :: local_nz, local_z_start
  integer :: total_local_size
#ifdef DOUB
  real(dp), dimension(total_local_size) :: phi1,phi2
#else
  real(sp), dimension(total_local_size) :: phi1,phi2
#endif
  integer :: myid, nproc, ierr
  integer :: seeds(MaxRandNumStreams,IRandNumSize)
  integer(i8b) :: index
  complex(dpc) :: z, ctemp, ctemp1,ctemp2

  real(dp) :: debug

  real lfm
  real(dp) :: twopi,avg,lavg,sigma,lsigma,chisq,lchisq
  parameter (twopi=6.283185307179586d0)
  integer j,i1,i2,i3,k1,k2,k3,k23,jp,modes
  integer(i8b) :: ndof
  real dk1,dk2,dk3,d3k,akmax,akmaxf,ak,ak23,akk,fact
  real(dp) :: xr
  real ak1,ak2,ak3,ak33,anu,dq,tf,tsav(0:npow),theta
  logical inflag, white_in
  integer mkmax
  parameter (mkmax=8196)
  integer ncamblines,ired,jndex
  real  kcamb(mkmax),pktcamb(mkmax,2),pkbcamb(mkmax,2),pkccamb(mkmax,2)
  real  ft(mkmax), fb(mkmax),fc(mkmax),hubble
  double precision dlnkh, minkh,akh
  common /camb/ kcamb,pktcamb,pkbcamb,pkccamb,ft,fb,fc,ncamblines,ired,dlnkh,minkh

  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)

  hubble = headc%h0/100.


  if (irand.lt.0.or.irand.gt.2) then
     if (myid==0) print*,'Error in ic4! Need irand=0,1,2'
     call mpi_finalize(ierr)
     stop
  end if
  if (itide.lt.-1.or.itide.gt.1) then
     if (myid==0) print*,'Error in ic4! Need itide=-1,0,1'
     call mpi_finalize(ierr)
     stop
  end if
  dk1=twopi/(np1*dx)
  dk2=twopi/(np2*dx)
  dk3=twopi/(np3*dx)
  d3k=dk1*dk2*dk3
  akmax=twopi/dx
  !  Precompute transfer function table for interpolation.
  !  N.B. must go to at least sqrt(3)*akmax/2 unless use Hanning filter,
  !  in which case go to akmax/2.
  if (hanning) then
     akmaxf=akmax/2.0
  else
     akmaxf=akmax
  end if
  do j=0,npow
     ak=j*akmaxf/npow
     tsav(j)=sqrt(pk(ak,astart)*d3k)
     if (hanning) tsav(j)=tsav(j)*cos(0.5*ak*dx)
     tsav(j)=tsav(j)/(ak*ak) ! velocity potential: phi
  end do
  tsav(0) = 0 ! zero padding
  !  Get white noise sample.
  if (irand.lt.2) then
     if (myid==0) print*,'Warning: Generating new random numbers in getPhi1Phi2!'
     call rans(nproc,iseed,seeds) ! Generate seeds for parallel random number generator
     call mpi_barrier(mpi_comm_world,ierr)
     do i3=1,local_nz
        do i2=1,np2
           do i1=1,np1
              index = int((i3-1)*np2+i2-1,8)*n2p1+i1
              call gaussdev(seeds(myid+1,:),xr)  ! actually f is not density but minus density
#ifdef DOUB
              phi1(index)=xr
#else
              phi1(index)=real(xr,kind=sp)
#endif
           end do
        end do
     end do
     ! f=1.0

     !  Output white noise field.
     if (irand.eq.1 .and. ipad.eq.0) then
        if (myid==0) then 
           call grafic_write_header_white(filename,np1,np2,np3,iseed)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call grafic_write(phi1,local_nz,local_z_start,np2,np1,filename,white_in=.true.)
     end if
  else
     !  irand=2 to read white noise from an input file
     if (myid==0) then
        print*,'Reading random numbers used in ic4 from ',trim(filename)
        call grafic_read_header_white(filename,i1,i2,i3,iseed)
     endif
     call mpi_barrier(mpi_comm_world,ierr)
     call mpi_bcast(i1,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(i2,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(i3,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,ierr)
     if (i1.ne.np1.or.i2.ne.np2.or.i3.ne.np3) then
        if (myid==0) then
           print*,'Error in ic4! file has np1,np2,np3=',i1,i2,i3
           print*,' Expected ',np1,np2,np3
        endif
        call mpi_finalize(ierr)
        stop
     end if
     if (myid==0) print*,'  Random numbers generated with iseed=',iseed
     call mpi_barrier(mpi_comm_world,ierr)
     call grafic_read(phi1,local_nz,local_z_start,np2,np1,filename,&
     &          padding_in=.false.,&
     &          white_in=.true.)
     ! Make sure input (big) white noise box has stdev=1
!    if(anorml .ne. 0.0) then
     if(myid==0) print *, 'renormalizing the input white noise'
     call mpnorm(local_nz,np3,np2,np1,total_local_size,phi1)
!    endif

  end if
  
  !  Compute mean.
  lavg=0.0
  do i3=1,local_nz
     do i2=1,np2
    	do i1=1,np1
           index = int((i3-1)*np2+i2-1,8)*n2p1+i1
#ifdef DOUB
           lavg=lavg+phi1(index)
#else
           lavg=lavg+real(phi1(index),kind=dp)
#endif
    	end do
     end do
  end do
  ! Mean is needed on every cpu
  !print*,'lavg = ',lavg
  call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  avg=avg/(real(np1*np2,dp)*np3)
  fact=sqrt(1.0*np1*np2*np3)


  !  Enforce zero mean for simulations with periodic boundary conditions.
  !  Compute chisq for this sample, too.
  lchisq=0.0
  do i3=1,local_nz
     do i2=1,np2
    	do i1=1,np1
           !  Subtract mean.
           index=int((i3-1)*np2+i2-1,8)*n2p1+i1
#ifdef DOUB
           phi1(index)=phi1(index)-avg
#else
           phi1(index)=phi1(index)-real(avg,kind=sp)
#endif
           lchisq=lchisq+phi1(index)**2
           !  Standard deviation is fact, but divide by np1*np2*np3=fact**2 to
           !  normalize f for FFT
#ifdef DOUB
           phi1(index)=phi1(index)/fact
#else
           phi1(index)=phi1(index)/real(fact,kind=sp)
#endif
    	end do
     end do
  end do

  call mpi_reduce(lchisq,chisq,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

  ndof=int(np1*np2,8)*np3-1
  anu=(chisq-ndof)/sqrt(float(ndof))
  if (myid==0) print*,'ic4 white noise: chisq, dof, nu=',real(chisq),ndof,anu
  !  Transform noise to k-space.
  call fft_mpi(plan,phi1,total_local_size)


! if(.true.) call psmeasure(f,total_local_size,local_nz,local_z_start,iseed)


  ! Paste if needed large scales, before applying transfer function
  ! This means that the large scales have been prewhitened ...
  if (ipad .eq. 1) then
     call grid_paste(local_z_start,local_nz,headt,headc,phi1,small_kfile_name)
     if (irand .eq. 1) then ! Output padded white noise file
        if (myid==0) then 
           call grafic_write_header_white(filename,np1,np2,np3,iseed)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call fft_mpi(iplan,phi1,total_local_size)
        phi1 = phi1/fact ! Renormalize to compensate backward FFT
        call grafic_write(phi1,local_nz,local_z_start,np2,np1,filename,white_in=.true.)
        call fft_mpi(plan,phi1,total_local_size)
        phi1 = phi1/fact ! Renormalize to compensate forward FFT
     endif
  endif


  !  Generate unconstrained sample in Fourier transform space.
  do k3=1,local_nz
     ak3=(k3+local_z_start-1)*dk3
     if (k3+local_z_start-1.gt.n32) ak3=ak3-akmax
     ak33=ak3*ak3
     do k2=1,np2
        ak2=(k2-1)*dk2
        if (k2-1.gt.n22) ak2=ak2-akmax
        ak23=ak2*ak2+ak33
        k23=k2-1+(k3-1)*np2
        do k1=1,n12+1 ! Complex 
           !  Do k1=n12+1 separately below.
           ak1=(k1-1)*dk1
           akk=ak1*ak1+ak23
           ak=sqrt(akk)
           !  Evaluate transfer function.
           dq=npow*ak/akmaxf
           if (dq.ge.npow) then
              tf=0.0
           else
              jp=int(dq)
              dq=dq-jp
              ! white noise is \delta
              ! here tsav = phi(1) = P(k)/k^2 
              tf=(1.0-dq)*tsav(jp)+dq*tsav(jp+1)
           end if
           !  Shift using offsets
           theta=ak1*xoff
           theta=theta+ak2*xoff
           theta=theta+ak3*xoff
           z=cmplx(cos(theta),sin(theta),kind=dpc)
           
           index = int((k3-1)*np2+k2-1,8)*n2p1+2*k1-1 !Real index
#ifdef DOUB
           ctemp = cmplx(phi1(index),phi1(index+1),kind=dpc)
#else
           ctemp = cmplx(real(phi1(index),kind=dp),real(phi1(index+1),kind=dp))
#endif
           ctemp = ctemp*z*tf

           !  Double the contribution to account for modes with k1 > n12+1 (k1 < 0).
           modes=2
           if (k1.eq.1) modes=1
           ! fill the output array, phi1
#ifdef DOUB
           phi1(index) = real(ctemp)
           phi1(index+1) = aimag(ctemp)
#else
           phi1(index) = real(ctemp,kind=sp)
           phi1(index+1) = real(aimag(ctemp),kind=sp)
#endif
        enddo
     enddo
  enddo
  call mpi_barrier(mpi_comm_world,ierr)
  !  Enforce zero mean.
  if (local_z_start==0) then
     phi1(1)=0
     phi1(2)=0
  endif

  if (myid==0) print*,'passed input power '

! print *, phi1(100), phi1(200), phi1(34), phi1(55)

! stop

  ! input phi(1), output phi(2) both of which are in Fourier space
  call twolpt(plan,iplan,phi1,total_local_size,local_nz,local_z_start,dx,phi2)
  if (myid==0) print*,'passed generating 2lpt data '


  return
end subroutine getPhi1Phi2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mpic4(icomp,idim,irand,iseed,itide,m1s,m2s,m3s,m1off, &
     &    m2off,m3off,hanning,filename,astart,mtype,dx,xoff,f,g,fm, &
     &    plan,iplan,local_nz,local_z_start,total_local_size, &
     &    headt,headc,small_kfile_name,ipad, anorml, phi1,phi2,omei)
  !  Generate an unconstrained sample of (rho,psi1,psi2,psi3,phi) for
  !  idim=0,1,2,3,4.
  !  Input: idim, irand, iseed, itide, m?s, m?off, hanning, filename,
  !         astart, pk, dx, xoff
  !  irand=0: use randg to generate white noise, don't save.
  !  irand=1: use randg to generate white noise, then save in real space
  !    in filename.
  !  irand=2: read filename to get random numbers.
  !  iseed: 9-digit integer random number seed.  Beware that rand8 does not
  !    give the same random numbers on 32-bit and 64-bit machines!
  !  itide=0 to use full subvolume for computing f.
  !  itide=1 to set xi=0 inside subvolume so as to get outer field.
  !  itide=-1 to set xi=0 outside subvolume so as to get inner field.
  !  m?s = size of next-level subvolume to split if itide.ne.0.
  !  m?off = offset of next-level subvolume to split if itide.ne.0
  !  hanning=T to apply hanning filter to f.
  !  hanning=F to not apply hanning filter to f.
  !  filename = file containing random numbers in real space.
  !  astart = expansion factor
  !  pk(ak,astart) = power spectrum function for wavenumber ak
  !  dx = grid spacing.
  !  xoff = offset to evaluate fields (e.g. use to shift baryon or cdm fields).
  !  Output: f=fc (sampled field in real space), fm (maximum absolute value of f).
  !  N.B. f and fc must point to the same place in memory - they are listed
  !  separately in the subroutine call because f77 will not allow equivalencing
  !  pointers.  The calling routine must pass the same pointer for each.
  use grafic_types
  use grafic_io
  use transform
  use normalize
  use random
  use paste

  implicit none
  include 'grafic1.inc'
  integer n12,n22,n32,n23,n2p1,npow
  parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
  parameter (n2p1=2*(np1/2+1))
  parameter (npow=30720)
  integer icomp, idim,irand,iseed,itide,m1s,m2s,m3s,m1off,m2off,m3off
  logical hanning
  character(len=128) :: filename, small_kfile_name
  character(len=*) :: mtype
  integer :: ipad
  real :: anorml
  type(taille) :: headt
  type(cosmo) :: headc
  real astart,pk,dx,xoff,fm,omei
  real zp1,Hastart,dladt
  integer(i8b) :: plan, iplan
  integer :: local_nz, local_z_start
  integer :: total_local_size
#ifdef DOUB
  real(dp), dimension(total_local_size) :: f,g,phi1,phi2
#else
  real(sp), dimension(total_local_size) :: f,g,phi1,phi2
#endif
  integer :: myid, nproc, ierr
  integer :: seeds(MaxRandNumStreams,IRandNumSize)
  integer(i8b) :: index
  complex(dpc) :: z, ctemp, ctemp1,ctemp2

  real(dp) :: debug

  real lfm
  real(dp) :: twopi,avg,lavg,sigma,lsigma,chisq,lchisq
  parameter (twopi=6.283185307179586d0)
  integer j,i1,i2,i3,k1,k2,k3,k23,jp,modes
  integer(i8b) :: ndof
  real dk1,dk2,dk3,d3k,akmax,akmaxf,ak,ak23,akk,fact
  real(dp) :: xr
  real ak1,ak2,ak3,ak33,anu,dq,f1,f2,tsav(0:npow),theta
  real D2f2Ha, D1f1Ha
  real D1,D2
  logical inflag, white_in
  external pk
  integer mkmax
  parameter (mkmax=8196)
  integer ncamblines,ired,jndex
  real  kcamb(mkmax),pktcamb(mkmax,2),pkbcamb(mkmax,2),pkccamb(mkmax,2)
  real, target::  ft(mkmax),fb(mkmax),fc(mkmax)
  real  hubble
  double precision dlnkh, minkh,akh
  common /camb/ kcamb,pktcamb,pkbcamb,pkccamb,ft,fb,fc,ncamblines,ired,dlnkh,minkh
  real omegam,omegav,h0,omegavz
  real w0,wa,cs2
  common /cosmoparms/ omegam,omegav,h0,w0,wa,cs2

  real, pointer :: growthRate(:)
  real*8 rinvcube


  rinvcube = 1.d0/(dble(np1)*dble(np2)*dble(np3))


  D1 = 1 ! This is because the input P(k) is at z_init not at z=0.
  D2 = -3./7. ! This scaling is correct irrepective of the wave number


  if(mtype == 'baryon')  then
     growthRate => fb
  else if (mtype == 'total') then
     growthRate => ft
  else if (mtype == 'cdm') then
     growthRate => fc
  else
     print *,'unknown matter type in mpi4c', mtype
     call mpi_finalize(ierr)
     stop
  endif

  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)

  hubble = headc%h0/100.


  if (idim.lt.0.or.idim.gt.4) then
     if (myid==0) print*,'Error in ic4! Need idim=0,1,2,3,4'
     call mpi_finalize(ierr)
     stop
  end if
  if (irand.lt.0.or.irand.gt.2) then
     if (myid==0) print*,'Error in ic4! Need irand=0,1,2'
     call mpi_finalize(ierr)
     stop
  end if
  if (itide.lt.-1.or.itide.gt.1) then
     if (myid==0) print*,'Error in ic4! Need itide=-1,0,1'
     call mpi_finalize(ierr)
     stop
  end if
  dk1=twopi/(np1*dx)
  dk2=twopi/(np2*dx)
  dk3=twopi/(np3*dx)
  d3k=dk1*dk2*dk3
  akmax=twopi/dx
  !  Precompute transfer function table for interpolation.
  !  N.B. must go to at least sqrt(3)*akmax/2 unless use Hanning filter,
  !  in which case go to akmax/2.
  if (hanning) then
     akmaxf=akmax/2.0
  else
     akmaxf=akmax
  end if

! peculiar velocity (in proper km/second) 
  zp1 = 1./astart
  omegavz = omegav*zp1**(3*(1+w0+wa))*exp(-3*wa*(zp1-1)/zp1)
  Hastart = headc%h0*sqrt(omegam*zp1**3+omegavz+(1-omegam-omegav)*zp1**2)*astart

  f2=2.*omei**(4./7.)
  D2f2Ha = D2*f2*Hastart

! if(myid ==0) then
!    do k1 = 1, ncamblines
!       print *, kcamb(k1), growthRate(k1), fb(k1), fc(k1)
!    enddo
! endif
! stop

  lchisq=0.0
  lsigma=0.0
  !  Generate unconstrained sample in Fourier transform space.
  do k3=1,local_nz
     ak3=(k3+local_z_start-1)*dk3
     if (k3+local_z_start-1.gt.n32) ak3=ak3-akmax
     ak33=ak3*ak3
     do k2=1,np2
        ak2=(k2-1)*dk2
        if (k2-1.gt.n22) ak2=ak2-akmax
        ak23=ak2*ak2+ak33
        k23=k2-1+(k3-1)*np2
        do k1=1,n12+1 ! Complex 
           index = int((k3-1)*np2+k2-1,8)*n2p1+2*k1-1 !Real index

           !  Do k1=n12+1 separately below.
           ak1=(k1-1)*dk1
           akk=ak1*ak1+ak23
           if(akk.le.0.) then
              f(index) = 0
              f(index+1) = 0
              g(index) = 0
              g(index+1) = 0
           else
              ak=sqrt(akk)
              akh = ak/hubble
              dq = log(akh/minkh)/dlnkh + 1
              jp=int(dq)
              dq=dq-jp
              f1=(1.0-dq)*growthRate(jp)+dq*growthRate(jp+1)
              D1f1Ha = D1*f1*Hastart

!             print *, jp, dq, growthRate(jp), f1,f2
!             print *, phi1(index), phi2(index)
!             stop

           
!--------------------------------------------
!---       this is for the velocity or density: f
#ifdef DOUB
              ctemp1 = cmplx(phi1(index),phi1(index+1),kind=dpc)
              ctemp2 = cmplx(phi2(index),phi2(index+1),kind=dpc)
#else
              ctemp1 = cmplx(real(phi1(index),kind=dp),real(phi1(index+1),kind=dp))
              ctemp2 = cmplx(real(phi2(index),kind=dp),real(phi2(index+1),kind=dp))
#endif

              if(idim.eq.0 .or. idim .eq. 4) then
                 ctemp = ctemp1*akk ! density. note the plus sign because of the velocity potential
              else if (idim.eq.1) then
                 ctemp=cmplx(0.0,1.0)*ak1*(D1f1Ha*ctemp1 +D2f2Ha*ctemp2) ! velocity
              else if (idim.eq.2) then
                 ctemp=cmplx(0.0,1.0)*ak2*(D1f1Ha*ctemp1 +D2f2Ha*ctemp2)
              else if (idim.eq.3) then
                 ctemp=cmplx(0.0,1.0)*ak3*(D1f1Ha*ctemp1 +D2f2Ha*ctemp2)
              end if
              ! fill the output array, f
#ifdef DOUB
              f(index) = real(ctemp)
              f(index+1) = aimag(ctemp)
#else
              f(index) = real(ctemp,kind=sp)
              f(index+1) = real(aimag(ctemp),kind=sp)
#endif
!--------------------------------------------
!---       this is for the displacement or density: g
              if(idim.eq.0 .or. idim .eq. 4) then
                 ctemp = ctemp1*akk ! density. note the plus sign because of the velocity potential
              else if (idim.eq.1) then
                 ctemp=cmplx(0.0,1.0)*ak1*(D1*ctemp1 +D2*ctemp2) ! displacement
              else if (idim.eq.2) then
                 ctemp=cmplx(0.0,1.0)*ak2*(D1*ctemp1 +D2*ctemp2)
              else if (idim.eq.3) then
                 ctemp=cmplx(0.0,1.0)*ak3*(D1*ctemp1 +D2*ctemp2)
              end if
              ! fill the output array, f
#ifdef DOUB
              g(index) = real(ctemp)
              g(index+1) = aimag(ctemp)
#else
              g(index) = real(ctemp,kind=sp)
              g(index+1) = real(aimag(ctemp),kind=sp)
#endif
           endif

        enddo
     enddo
  enddo
  call mpi_barrier(mpi_comm_world,ierr)
  !  Enforce zero mean.
  if (local_z_start==0) then
     f(1)=0
     f(2)=0
     g(1)=0
     g(2)=0
  endif

  if (myid==0) print*,'passed input power '


  !      Transform back to velocity in real space.
  call fft_mpi(iplan,f,total_local_size)
  !      Transform back to position in real space.
  call fft_mpi(iplan,g,total_local_size)

! f = f *rinvcube
! g = g *rinvcube

  
  
  lfm=0.0
  do i3=1,local_nz
     do i2=1,np2
        do i1=1,np1
           index=int((i3-1)*np2+i2-1,8)*n2p1+i1
#ifdef DOUB
           lfm=max(lfm,abs(f(index)))
#else
           lfm=max(lfm,real(abs(f(index)),kind=dp))
#endif
        end do
     end do
  end do
  ! Reduce max on root proc
  call mpi_reduce(lfm,fm,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  if (myid==0) then
     print*,'Statistics of ic4 for idim, itide=',idim,itide
     print*,'   Mean sigma, sampled sigma, maximum=',real(chisq), &
          &    real(sigma),fm
  endif

  return
end subroutine mpic4

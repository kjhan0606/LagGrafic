!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  input: f && output: g (or phi(2))
subroutine twolpt(plan,iplan,phi1, total_local_size,local_nz,local_z_start,dx,phi2)

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
  integer i,j,k,i1,i2,i3,k1,k2,k3,k23,jp,modes
  real dkx,dky,dkz

  real astart,pk,dx,xoff,fm
  integer(i8b) :: plan, iplan
  integer :: local_nz, local_z_start
  integer :: total_local_size
#ifdef DOUB
! real(dp), dimension(total_local_size) :: f,g
  real(dp), dimension(:), allocatable :: w1,w2
  real(dp), dimension(total_local_size):: phi1,phi2

#else
! real(sp), dimension(total_local_size) :: f,g
  real(sp), dimension(:), allocatable :: w1,w2
  real(sp), dimension(total_local_size):: phi1,phi2
#endif
  integer(i8b) :: index
  complex(dpc) :: z, ctemp
  integer :: myid, nproc, ierr
  integer :: seeds(MaxRandNumStreams,IRandNumSize)
  real(dp) ky,kz,sqrrincube,wave,wave2,invwave2,fourpi2,Lx,Ly,Lz,twopi
  parameter(twopi=6.283185307179586d0)

  fourpi2 = twopi*twopi

  allocate(w1(total_local_size))
  allocate(w2(total_local_size))

  Lx = dx*np1
  Ly = dx*np2
  Lz = dx*np3

  dkx = twopi/Lx
  dky = twopi/Ly
  dkz = twopi/Lz



  sqrrincube = dsqrt(1.d0/np1/np2/np3)**2




!  xx-yy component

  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n12
           index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
           w1(index)   = -fourpi2*i*i*phi1(index)/(Lx*Lx)
           w1(index+1) = -fourpi2*i*i*phi1(index+1)/(Lx*Lx)
           w2(index)   = -fourpi2*ky*ky*phi1(index)/(Ly*Ly)
           w2(index+1) = -fourpi2*ky*ky*phi1(index+1)/(Ly*Ly)
        enddo
     enddo
  enddo


  call fft_mpi(iplan,w1,total_local_size)
  call fft_mpi(iplan,w2,total_local_size)
  do index = 1, total_local_size
     phi2(index) = w1(index)*w2(index)*sqrrincube
  enddo


!  xx-zz component
  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n12
           index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
           w2(index)   = -fourpi2*kz*kz*phi1(index)/(Lz*Lz)
           w2(index+1) = -fourpi2*kz*kz*phi1(index+1)/(Lz*Lz)
        enddo
     enddo
  enddo
  call fft_mpi(iplan,w2,total_local_size)
  do index = 1, total_local_size
     phi2(index) = phi2(index) + w1(index)*w2(index)*sqrrincube
  enddo


!  y-y component
  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n12
           index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
           w1(index)   = -fourpi2*ky*ky*phi1(index)/(Ly*Ly)
           w1(index+1) = -fourpi2*ky*ky*phi1(index+1)/(Ly*Ly)
        enddo
     enddo
  enddo
  call fft_mpi(iplan,w1,total_local_size)
  do index = 1, total_local_size
     phi2(index) = phi2(index) + w1(index)*w2(index)*sqrrincube
  enddo


!  x-y component
  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n12
           index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
           w1(index)   = -fourpi2*ky*i*phi1(index)/(Ly*Lx)
           w1(index+1) = -fourpi2*ky*i*phi1(index+1)/(Ly*Lx)
        enddo
     enddo
  enddo
  call fft_mpi(iplan,w1,total_local_size)
  do index = 1, total_local_size
     phi2(index) = phi2(index) - w1(index)*w1(index)*sqrrincube
  enddo


!  x-z component
  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n12
           index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
           w1(index)   = -fourpi2*kz*i*phi1(index)/(Lz*Lx)
           w1(index+1) = -fourpi2*kz*i*phi1(index+1)/(Lz*Lx)
        enddo
     enddo
  enddo
  call fft_mpi(iplan,w1,total_local_size)
  do index = 1, total_local_size
     phi2(index) = phi2(index) - w1(index)*w1(index)*sqrrincube
  enddo


!  y-z component
  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n12
           index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
           w1(index)   = -fourpi2*kz*ky*phi1(index)/(Lz*Ly)
           w1(index+1) = -fourpi2*kz*ky*phi1(index+1)/(Lz*Ly)
        enddo
     enddo
  enddo
  call fft_mpi(iplan,w1,total_local_size)
  do index = 1, total_local_size
     phi2(index) = phi2(index) - w1(index)*w1(index)*sqrrincube
  enddo

  call fft_mpi(plan,phi2,total_local_size)
  ! final factor of -1/k^2
! do k=1,local_nz
!    kz = (k-1) + local_z_start
!    if(kz .gt. n32) kz = kz-np3
!    do j=1,np2
!       ky = j - 1
!       if(ky .gt. n22) ky = ky - np2
!       do i=0,n12
!          index = 2*i+1 + n2p1*int(j-1+np2*(k-1),8)
!          wave2 = -(i*i/(Lx*Lx)+ky*ky/(Ly*Ly)+kz*kz/(Lz*Lz))*fourpi2
!          invwave2 = 1.d0/wave2
!          if(wave2 .lt.0) then
!            phi2(index)     = phi2(index)*invwave2
!            phi2(index+1)   = phi2(index+1)*invwave2
!          else
!            phi2(index)     = 0
!            phi2(index+1)   = 0
!          endif
!	   ! return to the gradient
!          if(idim.eq.1) then
!            ctemp = cmplx(phi2(index),phi2(index+1),kind=dpc)
!            ctemp = ctemp*cmplx(0.,1.)*twopi*i/Lx
!            phi2(index) = real(ctemp)
!            phi2(index+1) = aimag(ctemp)
!          else if(idim.eq.2)then
!            ctemp = cmplx(phi2(index),phi2(index+1),kind=dpc)
!            ctemp = ctemp*cmplx(0.,1.)*twopi*ky/Ly
!            phi2(index) = real(ctemp)
!            phi2(index+1) = aimag(ctemp)
!          else if(idim.eq.3) then
!            ctemp = cmplx(phi2(index),phi2(index+1),kind=dpc)
!            ctemp = ctemp*cmplx(0.,1.)*twopi*kz/Lz
!            phi2(index) = real(ctemp)
!            phi2(index+1) = aimag(ctemp)
!          endif
!       enddo
!    enddo
! enddo

  deallocate(w1,w2)
  return
  end

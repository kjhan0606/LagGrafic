!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine psmeasure(f,total_local_size,local_nz,local_z_start,iseed)

  use grafic_types
  use grafic_io
  use transform
  use normalize
  use random
  use paste

  implicit none
  include 'grafic1.inc'
  integer n12,n22,n32,n23,n2p1,npow,idim
  parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
  parameter (n2p1=2*(np1/2+1))
  integer i,j,k,i1,i2,i3,k1,k2,k3,k23,jp,modes,iseed
  real dkx,dky,dkz

  real astart,dx,xoff,fm
  integer(i8b) :: plan, iplan
  integer :: local_nz, local_z_start
  integer :: total_local_size
#ifdef DOUB
  real(dp), dimension(total_local_size) :: f
#else
  real(sp), dimension(total_local_size) :: f
#endif
  integer(i8b) :: index
  complex(dpc) :: z, ctemp
  integer :: myid, nproc, ierr
  real(dp) ky,kz,sqrrincube,wave,wave2,invwave2,fourpi2,Lx,Ly,Lz,twopi
  parameter(twopi=6.283185307179586d0)
  double precision pk(10),delksq,wt1,npk(10)
  integer kbin1, kbin2

  fourpi2 = twopi*twopi

  Lx = dx*np1
  Ly = dx*np2
  Lz = dx*np3

  dkx = twopi/Lx
  dky = twopi/Ly
  dkz = twopi/Lz


  call mpi_comm_rank(mpi_comm_world,myid,ierr)

!c  sqrrincube = dsqrt(1.d0/np1/np2/np3)**2
  sqrrincube = 1

  do k = 1, 10
     pk(k) = 0
     npk(k) = 0
  enddo

  do k=1,local_nz
     kz = (k-1) + local_z_start
     if(kz .gt. n32) kz = kz-np3
     do j=1,np2
        ky = j - 1
        if(ky .gt. n22) ky = ky - np2
        do i=0,n2p1
           index = 2*i+1 + n2p1*(j-1+np2*(k-1))
           wave = sqrt(i*i+ky*ky+kz*kz)
           delksq = f(index)*f(index) + f(index+1)*f(index+1)
           delksq = delksq * sqrrincube
           kbin1 = nint(wave)
           wt1 = wave-kbin1
           if(wt1 .ge. 0) then
           kbin2 = kbin1 + 1
           else
           kbin2 = kbin1 - 1
           endif
           if(kbin1 < 10) then
              pk(kbin1+1) = pk(kbin1+1) + (1-wt1)*delksq
              npk(kbin1+1) = npk(kbin1+1) + 1-wt1
           endif
           if(kbin2 < 10) then
              pk(kbin2+1) = pk(kbin2+1) + wt1*delksq
              npk(kbin2+1) = npk(kbin2+1) + wt1
           endif

        enddo
     enddo
  enddo

  do i = 1, 10
     pk(i) = pk(i)/npk(i)
  enddo
  if(myid.eq.0) then
     open(12,file='PS.out', access='append', form='formatted')
     write(12,'(i8,10(1x,g12.5))')iseed,pk(1),pk(2),pk(3),pk(4),pk(5),pk(6),pk(7),pk(8),pk(9),pk(10)
     close(12)
  endif

  call mpi_finalize(ierr)
  stop

  end

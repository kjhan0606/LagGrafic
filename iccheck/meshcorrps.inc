c234567
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
      integer nrank,myrank,np
      integer imax,jmax,kmax
      integer ng,nspace
      integer initflag
      parameter(nrank = 32)
      parameter(boxsize= 1048.576)
      parameter (imax=2048,jmax=imax,kmax=imax)
c     parameter(mmax=imax,mc=mmax*mmax*mmax,mp=mmax/nrank*mmax*mmax*8)
      parameter(mmax=imax,mp=mmax/nrank*mmax*mmax*8)
      parameter (mspace=1)
      integer*2 hilo(1,1,1)
      integer*8 deficit,Ldeficit,nvert,Lnvert
      integer flag,ni,nj,nk,smoothflag
      real delta,guess,label,mass,mtot,rhoc,shift,vmfrac,vol,vtot,
     &     guess0,shift0
      character *1 vm,ip
      integer sendsize,recvsize
      integer stag,rtag
      integer istatus(MPI_STATUS_SIZE)
      integer dest,src
      real esmooth
      integer iarg
      common /block1/ deficit,Ldeficit,nvert,ni,nj,nk,flag,iarg,
     &       delta,guess,label,mass,esmooth,smoothflag,
     &       mtot,rhoc,shift,vmfrac,vol,vtot,guess0,shift0,initflag,
     &       hilo,vm,ip
      integer local_nz,local_z_start,local_ny_after_transpose
      integer local_y_start_after_transpose,total_local_size
      integer ierror
      integer*8 plan,iplan
      real size
      common /myfftwblock/ plan,iplan, local_nz,local_z_start,
     &          local_ny_after_transpose,
     &          local_y_start_after_transpose,
     &          total_local_size,myrank,ng,nspace,size
      real hubble,omega,omegaL,bias,smoothing
	  character*80 infile,outfile,outfile2
	  common /filenames/ infile,outfile,outfile2

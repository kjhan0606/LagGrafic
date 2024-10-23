!234567
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
      integer np
      integer ng,nspace
      integer initflag
      integer myrank,nrank, imax,jmax,kmax,mmax
	  common /mpi_size/ nrank, myrank,imax, jmax,kmax, mmax
      parameter (mspace=1)
      integer*2 hilo(1,1,1)
      integer*8 deficit,Ldeficit,nvert,Lnvert
      integer flag,ni,nj,nk,smoothflag
      real delta,guess,label,mass,mtot,rhoc,shift,vmfrac,vol,vtot, &
     &     guess0,shift0
      character *1 vm,ip
      integer sendsize,recvsize
      integer stag,rtag
      integer istatus(MPI_STATUS_SIZE)
      integer dest,src
      real esmooth
      integer iarg
      common /block1/ deficit,Ldeficit,nvert,ni,nj,nk,flag,iarg, &
     &       delta,guess,label,mass,esmooth,smoothflag, &
     &       mtot,rhoc,shift,vmfrac,vol,vtot,guess0,shift0,initflag, &
     &       hilo,vm,ip
      integer local_nz,local_z_start,local_ny_after_transpose
      integer local_y_start_after_transpose,total_local_size
      integer ierror
      integer*8 plan,iplan
      real size
      common /myfftwblock/ plan,iplan, local_nz,local_z_start, &
     &          local_ny_after_transpose, &
     &          local_y_start_after_transpose, &
     &          total_local_size,ng,nspace,size
      real hubble,omega,omegaL,bias,smoothing
	  character*80 infile,outfile,outfile2
	  common /filenames/ infile,outfile,outfile2

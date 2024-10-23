!c  grafic1.inc
!c  User must set parameters then make grafic1
!c
!c  Parameters for offsetting baryon and CDM fields by 0.5 grid spacing.
!c  The origin of the baryon (offvelb) and cdm (offvelc) with respect to (0,0,0).  
!c  In this version, we change the meaning of this parameter from the offset of velocity
!c  to the offset of different species with the origin of grid (0,0,0) 
!c  The origin of the matter will be offvelb/c*(0.5,0.5,0.5)*dx
!c  offvelb/c should be +/-1 or 0.
	integer offvelb,offvelc
!c	parameter (offvelb=0,offvelc=+1)
  	parameter (offvelb=+1,offvelc=0)
!c
!c  Initial density fluctuation amplitude on the ultimate subgrid scale.
!c  For large grids (> 10^7) one may wish to decrease this; for small
!c  refinement grids (< 10^6) and multiple levels one may wish to increase it.
	real sigstart
!c	parameter (sigstart=0.2)
!c
!c  Top grid mesh size.
!c
	integer np1,np2,np3
!c	parameter (np1=128,np2=128,np3=128)
!c  parameter (np1=256,np2=256,np3=256)
  	parameter (np1=512,np2=512,np3=512)
!c	parameter (np1=2048,np2=2048,np3=2048)
!c	parameter (np1=1024,np2=1024,np3=1024)

!c
!c np[1-3]:      particle lattice size for computing initial conditions.
!c                 All must be divisible by 2.  They need not be equal.
!c  The memory required by grafic1 is np1*np2*np3 words plus a few.

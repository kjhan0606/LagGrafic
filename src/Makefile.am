bin_PROGRAMS = mpgrafic

# set the include path found by configure
INCLUDES= $(all_includes)

# the library search path.
mpgrafic_LDFLAGS = -L/usr/local/lib
mpgrafic_SOURCES = grafic_types.f90 parallel_io.c grafic_io.f90 \
	transform.f90  paste.f90 random.f90 util.f90 time.f90 \
	mpnorm.f90 mpic4.f90 mppower.f90 mpgrafic.f90 twolpt.f90
mpgrafic_LDADD = -l@FFTWPREC@rfftw_mpi -l@FFTWPREC@fftw_mpi \
	 -l@FFTWPREC@rfftw -l@FFTWPREC@fftw
AM_FC = @MPIFC@
AM_FCFLAGS = -cpp -D@PRECISION@

clean:
	rm -rf $(bin_PROGRAMS) *.o *.mod

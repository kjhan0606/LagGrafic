#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#


#AR = ar -X64 rcv
AR = ar rcv
RANLIB = ranlib
#OPT = -O3 -q64 -qtune=pwr4 -qarch=pwr4
#OPT = -g 
#OPT = -g  -mcmodel=medium
OPT = -g  -mcmodel=medium
#OPT = -O1


#############################################
# List of compilation directives
#############################################
# 
# -DNTREEMIN=nnn  - minimum number of particles in tree - default NTREEMIN=4
# -DGROUPNUM=nnn  - tree grouping number - default GROUPNUM=10
# -DEPSSQ=fff     - softening length squared for PM force correction 
#                 - default=4.0e-8
# -DNMEG=nnn      - number of megabytes of storage per processor - default=40
# -DINCLUDE_TREE_FORCE - define if you want tree code corrections
# -DTREE               - define if you want tree code corrections
# -DCOMPACT               - run in single mass mode for maximum memory usage
# -DTREEFIX            - correct treebuilding code for closely spaced particles
# -DMOVIE              - include image rendering code
# -DACCOUT             - test mode - read in file 'rvtmp' and dump accelerations
# -DDEBUG              - output debugging information
# -DVERBOSE            - output even more debugging information
# -DBIT64              - use on machines using 64 bit addressing 
#                       - i.e. Compaq alpha
# -DAIX                - use on IBM SP3 machines
# -DDYNAMIC_LOAD_BALANCE  - use on systems with heterogeneous processors
# -DOBSERVATION_MODE   - save halo data at every ObservationStride'th step
# -DSAVESLICE          - save slice density file at every time step
######################################################
# CURRENTLY NOT CHECKED - April 23, 2001 - jjd
#############################################################
# -DENOUGH_SPACE   : if there is enough space for memory
#    -DFDA_FORCE_ARRAY_SWP  : force array + swap space
#    -DTSC_SWAP : tsc swap
#############################################################
#
#####################################
# System specific definitions
#####################################
# LINUX
#####################################
#FC = mpif77
#CC = mpicc
#INCLUDES = -I/user/kjhan/local/include
#
##FDFLAGS = -DINCLUDE_TREE_FORCE 
#DFLAGS = -DNMEG=350 -DCORRELATION_FUNCTION \
#	-DTREEFIX -Dkjhr -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE \
#	-DOBSERVATION_MODE  -DPMSEEDFORCE 
#	 #  -DSAVESLICE  -DSURVEY \
#	# -DINCLUDE_TREE_FORCE -DTREE
#FFLAGS = $(FDFLAGS) $(OPT)
#CFLAGS = $(OPT) $(DFLAGS)  
#LDFLAGS = $(OPT)
#LIBS = -L/user/kjhan/local/lib/ -L./PREWORK -L./WORKNEW \
#	-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw \
#	-lm -lpre -lhfind 
##################################
# PG Compiler
##################################
#FC = mpif77 
#CC = mpicc
#INCLUDES = -I/home/g004cbp/local64/include
#
## 64 bits flag #############
#FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DVarPM -DINITIALREADING 
##-WF,-DACCOUT -DINDVEPS
#DFLAGS = -Bstatic -DNMEG=800L -DCORRELATION_FUNCTION -Dkjhr  -DBIT64 \
#	-DLONGINDEX -DINCLUDE_TREE_FORCE -DTREE \
#	-DTREEFIX  -D_LARGE_FILES -DSAVESLICE \
#	-DNEW_SAVING_IN_PREWORK -DOBSERVATION_MODE   -DVarPM -DSIMPLESURVEY\
#	-DINITIALREADING
#
#FFLAGS = $(FDFLAGS) $(OPT) 
#CFLAGS = $(OPT) $(DFLAGS)  
#LDFLAGS = $(OPT) 
#LIBS = -L/home/kjhan/fftw/lib/ -L./PREWORK -L./WORKNEW \
#	-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw \
#	-lm -lhfind -lpre
##################################
# IBM SP3
##################################
FC = mpiifort
CC = mpiicc
INCLUDES = -I/home/kjhan/local/include/

# 64 bits flag #############
#FDFLAGS = -WF,-DAIX -WF,-DBIT64 -WF,-DINCLUDE_TREE_FORCE -WF,-DVarPM


FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DVarPM

#-WF,-DACCOUT -DINDVEPS
DFLAGS = -DNMEG=1100L -DCORRELATION_FUNCTION -Dkjhr  -DBIT64 \
        -DLONGINDEX -DINCLUDE_TREE_FORCE -DTREE \
        -DTREEFIX  -D_LARGE_FILES -DSAVESLICE \
        -DNEW_SAVING_IN_PREWORK -DOBSERVATION_MODE   -DVarPM \
		-DDEBUG -DPMSEEDFORCE  -DSIMPLESURVEY

#		-DSIMPLESURVEY\
#        -DINITIALREADING  
# -DPMSEEDFORCE -DDEBUG

FFLAGS = $(FDFLAGS) $(OPT)
CFLAGS = $(OPT) $(DFLAGS)
LDFLAGS = $(OPT)
LIBS = -L/home/kjhan/local/lib \
       -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw \
       -lm 


#FC = mpxlf_r
#CC = mpcc_r
#INCLUDES = -I/home/g004cbp/local64/include

# 64 bits flag #############
#FDFLAGS = -WF,-DAIX -WF,-DBIT64 -WF,-DCOMPACT -WF,-DINCLUDE_TREE_FORCE #-WF,-DACCOUT
#DFLAGS = -DNMEG=5500L -DCORRELATION_FUNCTION -Dkjhr  -DBIT64 -DCOMPACT \
#	-DMULTIMASS -DLONGINDEX \
#	-DINCLUDE_TREE_FORCE -DTREE \
#	-DTREEFIX  -D_LARGE_FILES -DSAVESLICE \
#	-DOBSERVATION_MODE -DPMSEEDFORCE  -DNEW_SAVING_IN_PREWORK  -DSURVEY \
#	 #-DACCOUT
#FFLAGS = $(FDFLAGS) $(OPT) -qextname 
#CFLAGS = $(OPT) $(DFLAGS)  
#LDFLAGS = $(OPT) 
#LIBS = -L/home/g004cbp/local64/lib/ -L./PREWORK -L./WORKNEW \
#	-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw \
#	-lm -lmpi -lhfind -lpre

# 32 bits flag #############
#FFLAGS = $(FDFLAGS) $(OPT) -qextname -q32 -qtune=pwr3 -qarch=pwr3 -bmaxdata:2000000000 -bmaxstack:2000000000
#CFLAGS = $(OPT) $(DFLAGS)  -q32 -qtune=pwr3 -qarch=pwr3 -bmaxdata:2000000000 -bmaxstack:2000000000
#LDFLAGS = $(OPT) -qextname -q32 -qtune=pwr3 -qarch=pwr3 -bmaxdata:2000000000 -bmaxstack:2000000000

#################################
# Compaq Alpha
#################################
#FC = f77
#CC = cc 
#INCLUDES = -I/home/kjhan/dolphin/fftw/include
#FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT
#DFLAGS = -DBIT64 -DNMEG=256 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT -DTREEFIX
#FFLAGS = $(FDFLAGS)  -fast -nofor_main
#CFLAGS = $(DFLAGS)  -fast
#LDFLAGS =  -fast -nofor_main
#LIBS = -L/home/kjhan/dolphin/fftw/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm -lz -lmpi -lm
##################################

#--- C Compiler information
#  Leave the rest untouched


#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c .f .F

#rules to build binary from source


.c.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.f.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.for.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.c.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.f.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.for.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.F.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

#--- Targets

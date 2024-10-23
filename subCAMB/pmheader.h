#define FILEHEADER "p1024.z"
#define indextype long  long
#define MPI_INDEX MPI_LONG_LONG
typedef struct pmparticletype {
#ifdef MULTIMASS
	float mass;
#endif
#ifdef INDTIME
	int tflag;
#endif

	float x,y,z,vx,vy,vz;

#ifdef INDEX
	indextype indx;
#endif
} pmparticletype;

typedef struct pmvrparticletype {
	float vx,vy,vz;
} pmvrparticletype;
typedef struct treeparticletype {
#ifdef MULTIMASS
	float mass;
#endif
#ifdef INDTIME
	int tflag;
#endif
	float x,y,z,vx,vy,vz;
#ifdef INDEX
	indextype indx;
#endif
	struct treeparticletype *next;
} treeparticletype;
float gettime();
#define TIMER_START(A) cputime0[A] = gettime();
#define TIMER_STOP(A)  cputime1[A] = gettime();
#define ELAPSED_TIME(A) (cputime1[A] - cputime0[A])

#define LCDM 0
#define VAR_DE 1

typedef struct simparameters{
	int myid,nid,IC;
	float omep,omeplam,omepb,omei,wlam0,wlam1;
	float CsDE;
#if CAMB == DE_enabled
	float  omepnu,omepnusterile,numnu;
	int numMsvnu, neutrino_hierarchy;
#endif
	float pcorr;
	float hubble,npow;
	float boxsize,bias8, biascmb,As,sigma8zi;
	int Pk_rescaling;
	float amax,anow;
	int nx,ny,nz,nspace;
	double lnx,lny,lnz;
	indextype mx,my,mz,mxmy;
	float zinit,astep,rsmooth,rth;
	float theta,sphere_radius,particle_radius;
	float ken,ktot,const0,poten0;
	float fact1,fact2,pfact;
	float mass;
	float fNL,gNL;
	int nstep,stepcount,stepnum,iseed,nskip;
	int powreadflag;
	char rvfilename[190],rvprefix[190];
	char powfilename[190],inpapkfilename[190];
	char DarkEnergyModel[30];
	float zmin,zmax;
	int local_nz,local_z_start;
	long long np;
	int xyzshiftflag, ptypesize, nviewer;
	char Viewerfile[190];
} SimParameters;
#ifdef DEFINE_SIM_PARA
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN SimParameters simpar;
EXTERN float cputime0[100],cputime1[100];
#undef EXTERN



#ifdef XYZDBL
#	ifdef INDEX
/*
#           define XofP(p) fmod((double)((p)->x)+(double)(((p)->indx%simpar.mx)*simpar.nspace)+simpar.nx,simpar.nx)
#           define YofP(p) fmod((double)((p)->y)+(double)((((p)->indx%simpar.mxmy)/simpar.mx)*simpar.nspace)+simpar.ny,simpar.ny)
#           define ZofP(p) fmod((double)((p)->z)+(double)(((p)->indx/simpar.mxmy)*simpar.nspace)+simpar.nz,simpar.nz)
 *
 */
#		define XofP(p) fmod((double)((p)->x)+(double)(((p)->indx%simpar.mx))+simpar.nx,simpar.nx)
#		define YofP(p) fmod((double)((p)->y)+(double)((((p)->indx%simpar.mxmy)/simpar.mx))+simpar.ny,simpar.ny)
#		define ZofP(p) fmod((double)((p)->z)+(double)(((p)->indx/simpar.mxmy))+simpar.nz,simpar.nz)
#	else
#		error the pair definition of XYZDBL and !INDEX  is not allowed in the current version.
#	endif
#else
#	define XofP(p) ((p)->x)
#	define YofP(p) ((p)->y)
#	define ZofP(p) ((p)->z)
#endif
#define MofP(p,lvsimpar) (lvsimpar.mass)



void appapsq_(float *, float *, float *, float *, float *, float *, float *, float *);

/* This is the Poisson Correction factor needed for the CPL and Quintessence model */
float interppcork(float );
float interppcorr(float );


#if DEMODEL == 0
#define getpcorr(a) (simpar.pcorr)
#elif DEMODEL == 1
#define getpcorr(a) (interppcorr(a))
#endif

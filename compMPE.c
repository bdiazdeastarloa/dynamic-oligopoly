/* ------------------------------------------------------------------------
 * Computation of a Markov-Perfect Equilibrium:
 * Continuous-time dynamic stochastic game with R&D and learning by doing.
 *
 * Based on original code for Windows by Uli Doraszelski, Feb 2006. 
 * 
 * Extended and edited for Mac/Linux by Bernardo Diaz de Astarloa, 
 * Jan 2015 @ Penn State. 
 *
 * Extension includes GSL and NLopt libraries to compute pricing policies 
 * using FOCs.
 * --------------------------------------------------------------------- */


#include "mex.h"
#include "matrix.h"
#include "mat.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdarg.h>
/* Load GSL libraries for root-finding routine. */
/* Un-comment this ifdef FOCNODERIV
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
*/
/* Load NLopt library for optimization routine. */
/* Un-comment this ifdef FOCDERIV */
#include <nlopt.h>

/* Compiler directives. To be modified by user. */
#define NMAX 10                             /* Max. # of firms allowed in the market. */
#define DMAX 20                             /* Max. # of commons states allowed. */
/* #define GJ 1 */                                /* Iteration method: Gauss-Jacobi. */
#define FOCDERIV 1                          /* Use NLopt, derivative-based. */
/* #define FOCNODERIV 1 */                        /*Use GSL: bracke ting-based. */
/* #define BUBBLE */
/* #define BACKWARD */
/* #define ALTERNATING */
/* #define PRINTITER */
/* #define DEBUG */

/* More compiler directives. Not to be modified by user. */
#define STRLEN 40                           /* Max. # of string characters for input MAT-file. */

/* Global variables. */
unsigned char D;                            /* # (demand,import price) common states. */
unsigned char dsize;                        /* # demand common states. */
unsigned char p0size;                       /* # import price common states. */
int N;                                      /* # firms. */
unsigned char M;                            /* # productivity states (i_j). */
unsigned long S;                            /* State space: # of unique states of the world. */
mxArray *binom_ptr; unsigned long *binom;   /* State space: Binomial coefficients for decoding/encoding. */
mxArray *state_ptr; unsigned char *state;   /* State space: List of unique states of the world. */
mxArray *pfor_ptr; double *pfor;            /* State space: Import price grid. */
mxArray *mkt_ptr; double *mkt;              /* State space: Market size. */
mxArray *mgc_ptr; double *mgc;              /* State space: Marginal cost vector (productivity). */

double rho;                                 /* Discount factor. */
double alpha;                               /* Demand: Price coefficient. */
double Phi_lo;                              /* Scrap value: lower bound. */
double Phi_hi;                              /* Scrap value: upper bound. */
double Phie_lo;                             /* Entry cost: lower bound. */
double Phie_hi;                             /* Entry cost: upper bound. */

mxArray *trans_ptr; double *trans;          /* Transition rates: (demand,import price) states. */
double delta;                               /* Transition rates: Incumbent productivity. */
double alpha1;                              /* Transition rates: Incumbent productivity. */
double alpha2;                              /* Transition rates: Incumbent productivity. */
double eta1;                                /* Transition rates: Incumbent lbd. */
double eta2;                                /* Transition rates: Incumbent lbd. */
double etax1;                               /* Transition rates: Incumbent exit. */
double etax2;                               /* Transition rates: Incumbent exit (exponential spec). */
double etae1;                               /* Transition rates: Potential entrant entry. */
double etae2;                               /* Transition rates: Potential entrant entry (exponential spec). */
unsigned char ie0;                          /* Transition rates and initial state: Potential entrant. */

char method[STRLEN];                        /* Program control: Iteration method. */
int maxiter;                                /* Program control: Maximum number of iterations. */
double tol;                                 /* Program control: Tolerance. */
int steps;                                  /* Program control: Number of steps in modified policy iteration. */
double lambdaV;                             /* Program control: Weight of updated value function. */
double lambdax;                             /* Program control: Weight of updated investment policy. */
double lambdap;                             /* Program control: Weight of updated pricing policy. */
double lambday;                             /* Program control: Weight of updated entry/exit policy. */

mxArray *V0_ptr; double *V0;				/* Starting values: Value function. */
mxArray *x0_ptr; double *x0;				/* Starting values: Investment policy. */
mxArray *p0_ptr; double *p0;                /* Starting values: Pricing policy. */
mxArray *z0_ptr; double *z0;                /* Starting values: Exit/entry policy. */

mxArray *info_ptr; double *info_;           /* Output: Performance information. */
mxArray *V1_ptr; double *V1;                /* Output: Value function. */
mxArray *x1_ptr; double *x1;                /* Output: Investment policy. */
mxArray *p1_ptr; double *p1;                /* Output: Pricing policy. */
mxArray *z1_ptr; double *z1;                /* Output: Entry/exit policy. */

typedef struct focpar {                     /* Pricing policy: FOC parameters. */
    double low;                             /* FOC: Bracketing lower bound. */
    double high;                            /* FOC: Bracketing upper bound. */
    unsigned char *inc;                     /* FOC: List of firms' states. */
    double cj;                              /* FOC: Marginal cost. */
    double dem;                             /* FOC: Demand. */
    double pfor;                            /* FOC: Current foreign price. */  
    int j;                                  /* FOC: Firm being solved. */ 
    double *price;                          /* FOC: Current price vector. */ 
    double *dV;                             /* FOC: Vector of differences in value functions. */
} focpar;                       

typedef struct profitpar {                  /* Profit function parameters. */
    unsigned char *inc;                     /* Profit: List of firms' states. */
    double cj;                              /* Marginal cost. */
    double dem;                             /* Profit: Demand. */
    double pfor;                            /* Profit: Current foreign price. */  
    int j;                                  /* Profit: Firm being solved. */ 
    double *price;                          /* Profit: Current price vector. */ 
} profitpar; 

/* Function prototypes. */
void init(const mxArray *prhs[]);
void iterGS(void);
void transdema(double *, unsigned char *, const unsigned char );
void transincu(double *, unsigned char *, const unsigned char, const double, double, const double);
void transentr(double *, unsigned char *, const double, const double, const unsigned char);
void summation_V1(double *, double *, const unsigned char *, double *[], unsigned char *[], const double *);
double makeprofits(int , profitpar *);

/* Hazard rates functional form definitions. */
double hazlbd(const double);
double hazrd(const double);
double hazex(const double);
double hazen(const double);

/* Optimization routine functions. */
double pstarfind(focpar);                   
double foc(double , void *);               
double price_obj(unsigned n, const double *x, double *grad, void *param);

/* Encoding functions. */
unsigned long encode(const unsigned char *);
void decode(unsigned char *, const unsigned long);
void sort(unsigned char *, int *);
void change(unsigned char *, const int, const unsigned char, int *);
void checkstate(const unsigned char *, const char *);
double checkparm(double, double, double, const char *, const char *);

/* Other functions */
int match(const char *[], const char *);
mxArray *read(MATFile *, const char *);
void write(MATFile *, const char *, const mxArray *);
double *create(void);
void destroy(double *);
void copy(double *, const double *);


/* Macros to retrieve state-dependent objects. */
#define	lookup_state(s) (state+(N+1)*(s))
#define	lookup_V0(s) (V0+N*(s)-1)
#define	lookup_x0(s) (x0+N*(s)-1)
#define	lookup_y0(s) (z0+N*(s)-1)
#define	lookup_p0(s) (p0+N*(s)-1)
#define	lookup_V1(s) (V1+N*(s)-1)
#define	lookup_x1(s) (x1+N*(s)-1)
#define	lookup_y1(s) (z1+N*(s)-1)
#define	lookup_p1(s) (p1+N*(s)-1)

/* Need to define 'max' and 'min' which are not defined in C. 
 * Comment this out if in Windows (maybe change to lower case too?). */
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)


/********************************************************************************/
/* Main routine. */
/********************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const char *methods[3] = {"GaussJacobiIteration", "GaussSeidelIteration", ""};
	time_t start_time, stop_time;
	double elapsed;

	/* Check for proper number (and type) of arguments. */	
    if (nrhs != 1)
	    mexErrMsgIdAndTxt( "MATLAB:compMPE:invalidNumInputs","One input required.");
    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:compMPE:inputNotStruct","Input must be a structure.");  
    
    /* Initialize. */
    init(prhs);

    /* Starting time. */
    time(&start_time);
    /* mexPrintf("\nStarting time is %s", ctime(&start_time)); */

    /* Call value function iteration routine. */
    switch (match(methods, method)) {
    case 0:
        mexErrMsgIdAndTxt("MATLAB:compMPE:invalidMethod","Gauss-Jacobi not implemented separately."); break;
    case 1:
        iterGS(); break;
    default:
        mexErrMsgIdAndTxt("MATLAB:compMPE:invalidMethod","Iteration method not recognized."); break;
	}

    /* Stopping time. */
    time(&stop_time);
    /* mexPrintf("\nStopping time is %s", ctime(&stop_time)); */
    elapsed = difftime(stop_time, start_time);
    /* mexPrintf("Elapsed time is %f min.\n", elapsed/60.0); */
    
    /* Allocate output.
     * This is better than writing to a file since it won't create problems
     * when parallelizing. */  
    plhs[0] = mxDuplicateArray(info_ptr); 
    plhs[1] = mxDuplicateArray(V1_ptr);
    plhs[2] = mxDuplicateArray(x1_ptr); 
    plhs[3] = mxDuplicateArray(p1_ptr);
    plhs[4] = mxDuplicateArray(z1_ptr);
    
	return;
}


/**************************************************************************
 * Main computational routines.
 *************************************************************************/

/* Initialize. 
 * Reads and assigns parameter values and structures from MATLAB. */
void init(const mxArray *prhs[])
{
	unsigned long s;

	/* mexPrintf("\nInitializing...\n\n"); */
    
    /* Get parameters and objects from input structure. 
     * In the future include warnings for potential errors when reading 
     * the structure (empty, types, etc.). */
    
    /* Scalars. */
    D = (unsigned char) checkparm(*mxGetPr(mxGetField(prhs[0],0,"D")), 1, DMAX, "D", "init");
    N = (int) checkparm(*mxGetPr(mxGetField(prhs[0],0,"N")), 2, NMAX, "N", "init");
    M = (unsigned char) checkparm(*mxGetPr(mxGetField(prhs[0],0,"M")), 1, UCHAR_MAX-1, "M", "init");
	S = (unsigned long) checkparm(*mxGetPr(mxGetField(prhs[0],0,"S")), 1, ULONG_MAX, "S", "init");        
    rho = *mxGetPr(mxGetField(prhs[0],0,"rho")); 
    alpha = *mxGetPr(mxGetField(prhs[0],0,"alpha"));
    Phi_hi  = *mxGetPr(mxGetField(prhs[0],0,"Phi_hi"));
    Phi_lo  = *mxGetPr(mxGetField(prhs[0],0,"Phi_lo"));
    Phie_hi = *mxGetPr(mxGetField(prhs[0],0,"Phie_hi"));
    Phie_lo = *mxGetPr(mxGetField(prhs[0],0,"Phie_lo"));
    alpha1 = *mxGetPr(mxGetField(prhs[0],0,"alpha1"));
    alpha2 = *mxGetPr(mxGetField(prhs[0],0,"alpha2"));
    delta = *mxGetPr(mxGetField(prhs[0],0,"delta"));
    eta1 = *mxGetPr(mxGetField(prhs[0],0,"eta1"));
    eta2 = *mxGetPr(mxGetField(prhs[0],0,"eta2"));
    etax1 = *mxGetPr(mxGetField(prhs[0],0,"etax1"));
    etax2 = *mxGetPr(mxGetField(prhs[0],0,"etax2"));
    etae1 = *mxGetPr(mxGetField(prhs[0],0,"etae1"));
    etae2 = *mxGetPr(mxGetField(prhs[0],0,"etae2")); 
    ie0 = (unsigned char) checkparm(*mxGetPr(mxGetField(prhs[0],0,"ie0")), 1, M, "ie0", "init");
    dsize = (unsigned char) *mxGetPr(mxGetField(prhs[0],0,"dsize"));
    p0size = (unsigned char) *mxGetPr(mxGetField(prhs[0],0,"p0size"));
    
    /* Program control parameters. */
    maxiter = (int) *mxGetPr(mxGetField(prhs[0],0,"maxiter"));
    tol = *mxGetPr(mxGetField(prhs[0],0,"tol"));
    steps = (int) *mxGetPr(mxGetField(prhs[0],0,"steps"));
	lambdaV = *mxGetPr(mxGetField(prhs[0],0,"lambdaV"));
    lambdax = *mxGetPr(mxGetField(prhs[0],0,"lambdax"));
    lambdap = *mxGetPr(mxGetField(prhs[0],0,"lambdap"));
    lambday = *mxGetPr(mxGetField(prhs[0],0,"lambday"));
    mxGetString(mxGetField(prhs[0],0,"method"), method, STRLEN);

    /* Arrays (matrices). */
    binom_ptr = mxGetField(prhs[0],0,"binom"); binom = (unsigned long*) mxGetData(binom_ptr);
    state_ptr = mxGetField(prhs[0],0,"state"); state = (unsigned char*) mxGetData(state_ptr);
    trans_ptr = mxGetField(prhs[0],0,"trans"); trans = mxGetPr(trans_ptr);
    
    mgc_ptr = mxGetField(prhs[0],0,"mgc"); mgc = mxGetPr(mgc_ptr);
    pfor_ptr = mxGetField(prhs[0],0,"pfor"); pfor = mxGetPr(pfor_ptr);
    mkt_ptr = mxGetField(prhs[0],0,"mkt"); mkt = mxGetPr(mkt_ptr);

    V0_ptr = mxGetField(prhs[0],0,"V0"); V0 = mxGetPr(V0_ptr);
    x0_ptr = mxGetField(prhs[0],0,"x0"); x0 = mxGetPr(x0_ptr);
    p0_ptr = mxGetField(prhs[0],0,"p0"); p0 = mxGetPr(p0_ptr);
    z0_ptr = mxGetField(prhs[0],0,"y0"); z0 = mxGetPr(z0_ptr);
    
    #ifdef GJ
    V1_ptr = mxDuplicateArray(V0_ptr); V1 = mxGetPr(V1_ptr);
    x1_ptr = mxDuplicateArray(x0_ptr); x1 = mxGetPr(x1_ptr);
    p1_ptr = mxDuplicateArray(p0_ptr); p1 = mxGetPr(p1_ptr);
    z1_ptr = mxDuplicateArray(z0_ptr); z1 = mxGetPr(z1_ptr);
    #else
    V1_ptr = V0_ptr; V1 = mxGetPr(V1_ptr);
    x1_ptr = x0_ptr; x1 = mxGetPr(x1_ptr);
    p1_ptr = p0_ptr; p1 = mxGetPr(p1_ptr);
    z1_ptr = z0_ptr; z1 = mxGetPr(z1_ptr);
    #endif

    /* Create solution info matrix. */ 
    info_ptr = mxCreateDoubleMatrix(1,6,mxREAL); info_ = mxGetPr(info_ptr);
    
	/* More error checking: States of the world must be sorted and within range. */
	/* for (s=0; s<S; s++) 
		checkstate(lookup_state(s), "init");*/

	return;
}

/* Gauss Seidel Iteration. 
 * Performs the value function iteration routine. */
void iterGS(void)
{
	int n, k, st, iter = 1, done = 0, pos[NMAX+1];
    int drctn = 1; /* Direction of counter: 1=iterate forward */
	unsigned long s, cntr;
	double SV[NMAX+1], SH[NMAX+1], diff, alldiff[NMAX+1], *Vup, *Vin, Vspec[NMAX+1], profits[NMAX+1];
    double *V, *x, *p, *y, temp1, temp2;
    #ifdef GJ
	double *V_, *x_, *p_, *y_;
    #else
	double V_[NMAX+1], x_[NMAX+1], p_[NMAX+1], y_[NMAX+1];
    #endif
    double tolV = 0.0, tolx = 0.0, tolp = 0.0, toly = 0.0;
    double *rates[NMAX+1], pfoc[NMAX+1];
    unsigned char *di, *agg[D+1], dip[NMAX+1], *nz[NMAX+1];
    unsigned char i, j, ind, we0;
    signed char die;
    double p0state, dstate, cj, q, sh, d, pstar;
    
    /*
    #ifdef GJ
    mexPrintf("\nIterating (Gauss-Jacobi)...\n\n");
    #else
    mexPrintf("\nIterating (Gauss-Seidel)...\n\n");
    #endif
    */
    
    /* If iterating Gauss-Jacobi: copy initial values. */
    #ifdef GJ
    memcpy(V1, V0, N*S*sizeof(double));
    memcpy(x1, x0, N*S*sizeof(double));
    memcpy(p1, p0, N*S*sizeof(double));
    memcpy(z1, z0, N*S*sizeof(double));
    #endif
    
    /* Allocate transition rates and index to nonzero rates.
     * 'D+2' makes sense to compute sums later; this was there's always a zero at
     * the end of the array. */
    rates[0] = mxCalloc(D+2, sizeof(double));
    nz[0] = mxCalloc(D+2, sizeof(unsigned char));
    for (n=1; n<=N; n++) {
        rates[n] = mxCalloc(M+1+2, sizeof(double));
        nz[n] = mxCalloc(M+1+2, sizeof(unsigned char));
    }
    
    /* Create matrix of (demand,import price) indices */
    for (n=1; n<=D; n++) {
        agg[n] = mxCalloc(2, sizeof(unsigned char));
    }
    for (i=1; i<=dsize; i++) {
        for (j=1; j<=p0size; j++) {
            ind = j + (i-1)*(p0size);
            agg[ind][0] = i;
            agg[ind][1] = j;
        }
    }

	#ifdef BACKWARD
	drctn = -1;
	#endif
    
    /* Start iterative procedure. */
    while (!done) {
    	
        /* Starts with a 'for' loop to cycle through states, which includes 
         * another 'for' loop to cycle through # of firms. */
		#ifdef ALTERNATING
		drctn *= -1;
		#endif
		for (cntr=0; s=(drctn>0) ? cntr : S-cntr-1, cntr<S; cntr++) {
            /* mexPrintf("State is %d-%d\n",s,iter); */
            /* Preallocate arrays with zeros. */
            memset(alldiff, 0, (NMAX+1)*sizeof(double));
            memset(pfoc, 0, (NMAX+1)*sizeof(double));
            
            /* Obtain unique state. */
            /* 'di' is N+1 vector (d,i1,...,iN). */
            /* dstate is demand level associated to d. */
            /* p0state is import prices level associated to d. */
			di = lookup_state(s);
            
            dstate = mkt[agg[di[0]][0]-1];
            p0state = pfor[agg[di[0]][1]-1];
            
			/* Store value and policy functions (previous iteration). */
            #ifdef GJ
            V_ = lookup_V0(s);
            x_ = lookup_x0(s);
            p_ = lookup_p0(s);
            y_ = lookup_y0(s);
            #else
            memcpy(V_, lookup_V0(s), (NMAX+1)*sizeof(double));
            memcpy(x_, lookup_x0(s), (NMAX+1)*sizeof(double));
            memcpy(p_, lookup_p0(s), (NMAX+1)*sizeof(double));
            memcpy(y_, lookup_y0(s), (NMAX+1)*sizeof(double));
            #endif
            
            /* Lookup value and policy functions. */
            /* These are pointers to the mxArrays to be written back to MATLAB. */
			V = lookup_V1(s);
			x = lookup_x1(s);
            p = lookup_p1(s);
			y = lookup_y1(s);

			/* Compute transition rates (demand state). */
			transdema(rates[0], nz[0], di[0]);
            
			/* Compute optimal policy. */
			for (n=1; n<=N; n++) {
				/* Done that before? Enforce symmetry and anonymity. */
				if ((n>1) && (di[n]==di[n-1])) {
					V[n] = V[n-1];
					x[n] = x[n-1];
                    p[n] = p[n-1];
					y[n] = y[n-1];
					Vspec[n] = Vspec[n-1];
					memcpy(rates[n], rates[n-1], (M+1+2)*sizeof(const double));
					memcpy(nz[n], nz[n-1], (M+1+2)*sizeof(unsigned char));

				/* Incumbent firm. */
				} else if (di[n] < M+1) {
                  
					/* Compute optimal investment policy. */
					x[n] = 0.0;
					if (di[n] < M) {
						memcpy(dip, di, (NMAX+1)*sizeof(unsigned char));
						change(dip, n, di[n]+1, pos);           
						Vup = lookup_V1(encode(dip));
						if ((diff = Vup[pos[n]]-V[n]) > 0.0)
							x[n] = pow(alpha1*alpha2*diff,1.0/(1.0-alpha2));
					}

                    /* Compute pricing policy. */
                    /* Compute deltaVs = Vup - V. */
                    /* (Given incumbent n, all other incumbents k increasing w. */
                    for (k=1; k<=N && di[k]<M; k++) {
                        memcpy(dip, di, (NMAX+1)*sizeof(unsigned char));
                        change(dip, k, di[k]+1, pos);           
                        Vup = lookup_V1(encode(dip));
                        alldiff[k] = Vup[pos[n]] - V[n];
                    }
                    cj = mgc[di[n]-1];
                    focpar par_arg = {
                        .low   = -20,                           /* FOC: Bracketing lower bound. */    
                        .high  = 200,                           /* FOC: Bracketing upper bound. */
                        .inc   = di,                            /* FOC: List of firm states. */
                        .cj    = cj,                            /* FOC: Marginal cost. */
                        .dem   = dstate,                        /* FOC: Demand. */
                        .pfor  = p0state,                       /* FOC: Current foreign price. */  
                        .j     = n,                             /* FOC: Firm being solved. */ 
                        .price = p,                             /* FOC: Current price vector. */ 
                        .dV    = alldiff,                       /* FOC: Value function differences. */
                    };
                    pstar = pstarfind(par_arg);
                    p[n] = pstar;

					/* Compute optimal exit policy: y = 1 - G(V). */
					y[n] = 0.0;
                    temp1  = (Phi_hi-V[n])/(Phi_hi-Phi_lo);
					y[n] = MAX(0, MIN(temp1,1));

					/* Dampen optimal investment, pricing, and exit policies. */
					x[n] = lambdax*x[n]+(1.0-lambdax)*x_[n]; 
                    p[n] = lambdap*p[n]+(1.0-lambdap)*p_[n]; 
					y[n] = lambday*y[n]+(1.0-lambday)*y_[n]; 
                    
                    /* Compute quantity for firm n LBD rate. */
                    d = 0;
                    for (i=1; i<=N && di[i]<M+1; i++)
                        d += exp((-alpha)*(p[i] - p0state));
                    sh = exp((-alpha)*(p[n] - p0state))/(1+d);
                    q = (dstate)*(sh);

					/* Compute transition rates (firm state). */
					transincu(rates[n], nz[n], di[n], x[n], q, y[n]);

					/* Store value of exit = E(Phi | Phi > V). */
                    temp1 = (Phi_hi + V[n])/2;
                    temp2 = (Phi_hi + Phi_lo)/2;
					Vspec[n] = MIN(Phi_hi, MAX(temp1, temp2));

				/* Potential entrant. */
				} else {
					/* Compute optimal investment policy. */
					x[n] = 0.0;
                    
                    /* Compute optimal pricing policy. */
                    p[n] = 0.0;

					/* Compute optimal entry policy = G(Ve). */
					y[n] = 0.0;
                    
                    /* Alternative entrant productivity: -2 from leader. */
                    k=1;
                    while (di[k]<M || k==N)
                        k++;
                    die = (signed char) di[k];
                    we0 = MAX(1,die-2);
					memcpy(dip, di, (NMAX+1)*sizeof(unsigned char));
					change(dip, n, we0, pos);
					Vin = lookup_V1(encode(dip));
					temp1 = (Vin[pos[n]]-Phie_lo)/(Phie_hi - Phi_lo);
                    y[n] = MAX(0, MIN(1, temp1));
						
					/* Dampen optimal investment, pricing, and entry policies. */
					x[n] = lambdax*x[n]+(1-lambdax)*x_[n];
                    p[n] = lambdap*p[n]+(1-lambdap)*p_[n];
					y[n] = lambday*y[n]+(1-lambday)*y_[n];

					/* Compute transition rates (firm state). */
					transentr(rates[n], nz[n], x[n], y[n], we0);

					/* Store value of entry: - E(Phie | Phie>Vin) + Vin. */
                    temp1 = (Vin[pos[n]] + Phie_lo)/2;
                    temp2 = (Phie_hi + Phie_lo)/2;
					Vspec[n] = Vin[pos[n]] - MAX(Phie_lo, MIN(temp1, temp2));

				}	/* if (di[n] < M+1) */
               
			}	/* for (n=1; n<=N; n++) */
           
			/* Compute sums of values and rates.
             * All optimal policies have been computed already and will not
             * change from this point of the loop on. */
			summation_V1(SV, SH, di, rates, nz, Vspec);

			/* Update value function (preliminary). */
			for (n=1; n<=N; n++) {

				/* Incumbent firm. */
				if (di[n] < M+1) {
                    
                    /* Compute profits.
                     * Note that 'p' passed in '.price' contains the dampened
                     * policies computed above. */
                    cj = mgc[di[n]-1];
                    profitpar prof_par = {
                        .inc = di,
                        .cj  = cj,                         
                        .dem = dstate,                            
                        .pfor= p0state, 
                        .price = p,                           
                    };
                    profits[n] = makeprofits(n, &prof_par);
                    
					/* Update value function (preliminary). */
					V[n] = (profits[n]-x[n]+SV[n])/(rho+SH[n]);

				/* Potential entrant. */
				} else {

					/* Update value function (preliminary). */
					V[n] = SV[n];

				}	/* if (di[n] < M+1) */

				/* Dampen updated value function. */
				V[n] = lambdaV*V[n]+(1.0-lambdaV)*V_[n];

			}	/* for (n=1; n<=N; n++) */

            /* Tolerances. */
			for (n=1; n<=N; n++) {
				tolV = MAX(tolV, fabs((V[n]-V_[n])/(1.0+fabs(V[n])))); 
				tolx = MAX(tolx, fabs((x[n]-x_[n])/(1.0+fabs(x[n]))));
                tolp = MAX(tolp, fabs((p[n]-p_[n])/(1.0+fabs(p[n]))));
				toly = MAX(toly, fabs((y[n]-y_[n])/(1.0+fabs(y[n]))));
			}

        }	/* for (cntr=0; s=(drctn>0) ? cntr : S-cntr-1, cntr<S; cntr++) */

		/* Perform iteration step: Update value function (final). */
        for (st=0; st<steps; st++) {
			#ifdef ALTERNATING
			drctn *= -1;
			#endif
			for (cntr=0; s=(drctn>0) ? cntr : S-cntr-1, cntr<S; cntr++) {

	            /* Obtain unique state. */
				di = lookup_state(s);
                dstate = mkt[agg[di[0]][0]-1];
                p0state = pfor[agg[di[0]][1]-1];

				/* Lookup value and policy functions and current profit. */
				V = lookup_V1(s);
				x = lookup_x1(s);
                p = lookup_p1(s);
				y = lookup_y1(s);

				/* Compute transition rates (demand state). */
				transdema(rates[0], nz[0], di[0]);

				for (n=1; n<=N; n++) {

					/* Done that before? Enforce symmetry and anonymity. */
					if ((n>1) && (di[n]==di[n-1])) {
						Vspec[n] = Vspec[n-1];
						memcpy(rates[n], rates[n-1], (M+1+2)*sizeof(double));
						memcpy(nz[n], nz[n-1], (M+1+2)*sizeof(unsigned char));
					
					/* Incumbent firm. */
					} else if (di[n] < M+1) {

						/* Compute scrap value. */
                        temp1 = (Phi_hi + V[n])/2;
                        temp2 = (Phi_hi + Phi_lo)/2;
                        
                        /* Compute quantities for LBD rates. */
                        d = 0;
                        for (i=1; i<=N && di[i]<M+1; i++)
                                d += exp((-alpha)*(p[i] - p0state));
                        sh = exp((-alpha)*(p[n] - p0state))/(1+d);
						q = (dstate)*(sh);
                        
                        /* Compute transition rates (firm state). */
						transincu(rates[n], nz[n], di[n], x[n], q, y[n]);

						/* Store value of exit. */
						Vspec[n] = MIN(Phi_hi, MAX(temp1, temp2));;

					/* Potential entrant. */
					} else {
                        
                        k=1;
                        while (di[k]<M || k==N)
                            k++;
                        die = (signed char) di[k];
                        we0 = MAX(1,die-2);

						/* Compute transition rates (firm state). */
						transentr(rates[n], nz[n], x[n], y[n], we0);
                        
                        /* Store value of entry. */
                        temp1 = (Vin[pos[n]] + Phie_lo)/2;
                        temp2 = (Phie_hi + Phie_lo)/2;
					    Vspec[n] = -MAX(Phie_lo, MIN(temp1, temp2)) + Vin[pos[n]];

					}	/* if (di[n] < M+1) */

				}	/* for (n=1; n<=N; n++) */

				/* Compute sums of values and rates. */
				summation_V1(SV, SH, di, rates, nz, Vspec);

				/* Update value function (final). */
				for (n=1; n<=N; n++) {

					/* Done that before? Enforce symmetry and anonymity. */
					if ((n>1) && (di[n]==di[n-1])) {
						V[n] = V[n-1];
					
					/* Incumbent firm. */
					} else if (di[n] < M+1) {
                        
                        /* Compute profits. 
                         * We could omit this step if we stored the profits that were computed in
                         * the preliminary update; policies have not changed, so profits won't. 
                         * Probably does not save much time, but it will take up memory. */
                        cj = mgc[di[n]-1];
                        profitpar prof_par = {
                            .inc = di,
                            .cj  = cj,                         
                            .dem = dstate,                            
                            .pfor= p0state, 
                            .price = p,                          
                        };
                        profits[n] = makeprofits(n, &prof_par);

						/* Update value function (final). */
						V[n] = (profits[n]-x[n]+SV[n])/(rho+SH[n]);

					/* Potential entrant. */
					} else {

						/* Update value function (final). */
						V[n] = SV[n];

					}	/* if (di[n] < M+1) */

				}	/* for (n=1; n<=N; n++) */

            }   /* for (cntr=0; s=(drctn>0) ? cntr : S-cntr-1, cntr<S; cntr++) */
        }   /* for (st=0; st<steps; st++) */

        /* Tolerances. */
        /* 
        if (iter-(maxiter/5)*(iter/(maxiter/5))==0)
            mexPrintf("Iteration %d: tol(V, x, p, y)=(%g, %g, %g, %g).\n", iter, tolV, tolx, tolp, toly);
        */
        
        /* Done? */
        if ((tolV < tol) && (tolx < tol) && (tolp < tol) && (toly < tol)) {
            mexPrintf("Converged!\n");
            mexPrintf("Iteration %d: tol(V, x, p, y)=(%g, %g, %g, %g).\n", iter, tolV, tolx, tolp, toly);
            done = 1;
        } else if (iter >= maxiter) {
            mexPrintf("Maximum number of iterations reached!\n");
            mexPrintf("Iteration %d: tol(V, x, p, y)=(%g, %g, %g, %g).\n", iter, tolV, tolx, tolp, toly);
            done = 2;
        }

        /* Updating. */
        if (!done) {
            ++iter;
            tolV = tolx = tolp = toly = 0.0;
            #ifdef GJ
            copy(V0, V1);
            copy(x0, x1);
            copy(p0, p1);
            copy(z0, z1);
            #endif
        }

    }   /* while (!done) */

    /* Performance information. */
    info_[0] = (double) done;
    info_[1] = (double) iter;
    info_[2] = tolV;
    info_[3] = tolp;
    info_[4] = tolx;
    info_[5] = toly;

    /* Free memory of structures. */
    for (n=0; n<=N; n++) {
        mxFree(rates[n]);
        mxFree(nz[n]);
    }
    for (i=1; i<=D; i++) {
        mxFree(agg[i]);
    }
    return;
}


/**************************************************************************
 *                  Auxiliary function definitions.                        
 *************************************************************************/


/**************************************************************************
 * Define hazard rates of firm-specific jump events:
 * LBD, R&D, entry, exit.                        
 *************************************************************************/

/* hazlbd.
 * Compute hazard of LBD-productivity jump. */
double hazlbd(const double q)
{
    double h = eta1*q/(1+eta1*q);
    /* double h = eta1*pow(q,eta2); */
    return h;
}
/* hazlbdprime.
 * Compute derivative of hazard of LBD-productivity jump. */
double hazlbdprime(const double q)
{
    double hq = eta1/pow(1+eta1*q,2); 
    /* double hq = eta2*eta1*pow(q,eta2-1); */
    return hq;
}
/* hazrd.
 * Compute hazard of R&D-productivity jump. */
double hazrd(const double x)
{
    double h = alpha1*pow(x,alpha2);
    return h;
}
/* hazex.
 * Compute hazard of exiting. */
double hazex(const double y)
{
    double h = etax1*y;
    return h;
}
/* hazen.
 * Compute hazard of exiting. */
double hazen(const double y)
{
    double h = etae1*y;
    return h;
}


/**************************************************************************
 * Optimization routines.
 * Depending on the specification of the learning by doing hazard, 
 * a gradient-based optimization routine for the Bellman equation may be 
 * more convenient than a root-finding routine for the first-order 
 * condition. (The function can have local min and max.)
 *
 * Here I include two optional routines:
 * - A GSL root-finding solver.
 * - A NLopt optimization solver.                  
 *************************************************************************/

/* pstarfind. 
 * GSL routine to solve FOC for prices. */
#ifdef FOCNODERIV
double pstarfind(focpar par_arg)
{
  int status;
  double error=99.0;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;                       /* Pointer to solver type */
  gsl_root_fsolver *solv;                               /* Pointer to solver */
  double r = 0;                                         /* Initialize root */
  double p_lo = par_arg.low, p_hi = par_arg.high;       /* Define bounds */
  
  gsl_function F;                           
  F.function = &foc;
  F.params = &par_arg;
  
  /* Assign solver type, allocate and initialize */
  T = gsl_root_fsolver_brent;              
  solv = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(solv, &F, p_lo, p_hi); 
  
  /* Start iterations */  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate(solv);  /* Perform an iteration */ 
      r = gsl_root_fsolver_root(solv);          /* Current estimate of root */
      p_lo = gsl_root_fsolver_x_lower(solv);    /* Current estimate of lower bound */ 
      p_hi = gsl_root_fsolver_x_upper(solv);    /* Current estimate of upper bound */ 

      /* Test for convergence */
      status = gsl_root_test_interval(p_lo, p_hi, 0, 1e-6);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  return r;
  
  /* Free memory */
  gsl_root_fsolver_free(solv);                     
}
#endif


/* pstarfind. 
 * NLopt routine to solve for prices using Bellman. */
#ifdef FOCDERIV
double pstarfind(focpar par_arg)
{
    double lb[1]; 
    lb[0] = -200;                                      /* Lower bound. */
    double tol_x = 1e-5;                               /* Convergence tolerance. */
    double out;
    
    /* Algorithm and dimensionality. */
    nlopt_opt opt= nlopt_create(NLOPT_LD_MMA, 1);
    /* Set bounds (upper bound not needed). */
    nlopt_set_lower_bounds(opt, lb);
    /* Assign objective function to maximize. */
    nlopt_set_max_objective(opt, price_obj, &par_arg);    
    /* Assign tolerances. */
    nlopt_set_xtol_rel(opt, tol_x);
    
    /* Start optimization. */
    double x[1];                                       /* Initial guess. */
    x[0] = 2;
    double maxf;                                       /* Objective value. */

    if (nlopt_optimize(opt, x, &maxf) < 0) {
        printf("NLopt failed!\n");
        out = 10.0;
    }
    else {
        /* printf("Max. found at f(%g) = %0.10g\n", x[0], maxf); */
        out = x[0];
    }
    
    /* Free memory. */
    nlopt_destroy(opt);
    
    return out;
}
#endif


/* foc. 
 * Compute first order condition for pricing policy. 
 * Input to GSL root-finding function. */
double foc(double x, void *param)
{
    int i;
    double d, q, hsV, s[NMAX+1], pp[NMAX+1], hq[NMAX+1];
    
    focpar *par   = (focpar *) param;
    unsigned char *inc = par->inc;
    double p_0    = par->pfor;
    double cj     = par->cj;
    double dem    = par->dem;
    int j         = par->j;
    double *price = par->price;
    double *diff_ = par->dV;
    
    /* Need to create local copy not to modify last iteration's price. */
    for (i=1; i<=N; i++) {
        pp[i] = *(price+i);
        hq[i] = 0;
        s[i]  = 0;
    }
    pp[j] = x;
    d = 0.0;
    hsV = 0.0;
    q = 0;
    /* Compute denominator of mkt shares. */
    for (i=1; i<=N && inc[i]<M+1; i++) 
        d += exp((-alpha)*(pp[i] - p_0));
   /* Compute market shares and sum of DVs. */
    for (i=1; i<=N && inc[i]<M+1; i++) { 
        s[i] = exp((-alpha)*(pp[i] - p_0))/(1+d);
        q    = dem*s[i];
        hq[i] = hazlbd(q);
        hsV += hq[i]*(s[i])*(diff_[i]);
    }
    
    /* Return FOC */
    return (1/alpha) - (1-s[j])*(x-cj) - hq[j]*(diff_[j]) + hsV;
}

/* price_obj.
 * Compute objective function for pricing policy and gradient. 
 * Input to NLopt maximization routine. */
double price_obj(unsigned n, const double *x, double *grad, void *param)
{
    int i;
    double d, hsV1, hsV2, q[NMAX+1], s[NMAX+1], pp[NMAX+1], h[NMAX+1], hq[NMAX+1];
    
    focpar *par   = (focpar *) param;
    unsigned char *inc = par->inc;
    double p_0    = par->pfor;
    double cj     = par->cj;
    double dem    = par->dem;
    int j         = par->j;
    double *price = par->price;
    double *diff_ = par->dV;
    
    /* Need to create local copy not to modify last iteration's price. 
     * Also, assign zeros to make sure nothing weird happens to arrays. */
    for (i=1; i<=N; i++) {
        q[i]  = 0.0;
        s[i]  = 0.0;
        h[i]  = 0.0;
        hq[i] = 0.0;
        pp[i] = *(price+i);
    }
    pp[j] = x[0];
    d = 0.0;
    hsV1 = 0.0;
    hsV2 = 0.0;
    /* Compute denominator of mkt shares. */
    for (i=1; i<=N && inc[i]<M+1; i++) 
        d += exp((-alpha)*(pp[i] - p_0));
   /* Compute market shares, hazards and sum of deltaVs. */
    for (i=1; i<=N && inc[i]<M+1; i++) { 
        s[i] = exp((-alpha)*(pp[i] - p_0))/(1+d);
        q[i] = s[i]*dem;
        h[i] = hazlbd(q[i]);
        hq[i] = hazlbdprime(q[i]);
        hsV1 += hq[i]*(s[i])*(diff_[i]);
        hsV2 += (h[i])*(diff_[i]);
    }
    
    /* Assign gradient. */
    if (grad) {
        grad[0] = (1/alpha) - (1-s[j])*(pp[j]-cj) - hq[j]*(diff_[j]) + hsV1;
    }
    /* Return objective function. */
    return q[j]*(pp[j]-cj) + hsV2;
}

/* ----------------------------------------------------------------------*/
/* ----------------------------------------------------------------------*/

/* makeprofits.
 * Compute profit function. */
double makeprofits(int j, profitpar *params)
{
    int i;
    double d, s;
    unsigned char *inc = params->inc;
    double dem = params->dem;
    double p_0 = params->pfor;
    double cj  = params->cj;
    double *pp = params->price;
    
    d = 0.0;
    /* Compute denominator of mkt shares */
    for (i=1; i<=N && inc[i]<M+1; i++)
        d += exp((-alpha)*(pp[i] - p_0));
    
    s = exp((-alpha)*(pp[j] - p_0))/(1+d);
    
    return (dem)*(s)*(pp[j] - cj); 
}


/* ----------------------------------------------------------------------*/
/* ----------------------------------------------------------------------*/


/* transdema.
 * Calculate common state transition rates. */
void transdema(double *rates, unsigned char *nz, const unsigned char d)
{
    unsigned char dp;
    
    #ifdef DEBUG
	if (!((d >= 1) && (d <= D)))
        mexErrMsgTxt("Undefined transition rates of demand states!");
	#endif
        
	#ifdef DEBUG
    memset(rates, 0, (D+2)*sizeof(double));
    memset(nz, 0, (D+2)*sizeof(unsigned char));
	#endif

    for (dp=1; dp<=D; dp++) 
        if (d != dp) {                              /* Transit to a different state. */
            rates[dp] = *(trans+(d-1)+(dp-1)*D);    /* Transition rate (pointer to element of rates[0]). */
            if (rates[dp] > 0.0)                    /* Assign states with non-zero rates to nz[0]. */
                *(++nz) = dp;
    }
	*(++nz) = 0;
}


/* ----------------------------------------------------------------------*/
/* ----------------------------------------------------------------------*/


/* transincu.
 * Calculate transition rates for incumbent firms. */
void transincu(double *rates, unsigned char *nz, const unsigned char i, const double x, const double q, const double y)
{
	#ifdef DEBUG
	if (!((i >= 1) && (i <= M)))
        mexErrMsgTxt("Undefined transition rates of incumbent firm!");
	#endif

	#ifdef DEBUG
	memset(rates, 0, (M+1+2)*sizeof(double));
    memset(nz, 0, (M+1+2)*sizeof(unsigned char));
	#endif
    
    /* Lower productivity level: can't go down, only up to i=2. */
    if (i == 1) {
        rates[2] = hazrd(x) + hazlbd(q);
        *(++nz) = 2;
    }
    /* Higher productivity level: can't go up, only down to i=M-1. */
    else if (i == M) {
        rates[M-1] = delta;
        *(++nz) = M-1; 
    }
    /* All other productivity levels. */
    else {
        rates[i-1] = delta;                     /* Hazard of going down. */
        rates[i+1] = hazrd(x) + hazlbd(q);      /* Hazard of going up. */
        *(++nz) = i-1; *(++nz) = i+1;
    } 
    
    /* Exit rate. */
    if (y > 0.0) {
        rates[M+1] = hazex(y);
        *(++nz) = M+1;
    }

	*(++nz) = 0;
}


/* ----------------------------------------------------------------------*/
/* ----------------------------------------------------------------------*/


/* transentr.
 * Calculate transition rates for potential entrants. */
void transentr(double *rates, unsigned char *nz, const double xe, const double ye, const unsigned char we0)
{
	#ifdef DEBUG
	memset(rates, 0, (M+1+2)*sizeof(double));
    memset(nz, 0, (M+1+2)*sizeof(unsigned char));
	#endif

    if (ye > 0.0) {
        rates[we0] = hazen(ye);
        *(++nz) = we0;
    }

	*(++nz) = 0;
}

 
/* ------------------------------------------------------------------------
 * Compute sums of values and rates. 
 * ------------------------------------------------------------------------
 *
 * Note: rates[0] and nz[0] are vectors with D+2 elements and describe the 
 * transition rates for the demand state. rates[n] and nz[n] are vectors
 * with M+1+2 elements and describe the transition rates for firm n's state.
 * nz indexes the nonzero rates: nz[n][1]>0 indexes the first nonzero 
 * rate, nz[n][2]>0 the second one, and so on until nz[n][m]=0. The sum of 
 * values uses Vspec instead of the value function to properly account for the 
 * value of entry and exit. */


/* summation_V1.
 * Compute sums of values (SV) and hazard rates (SH). */
void summation_V1(double *SV, double *SH, const unsigned char *di, double *rates[], unsigned char *nz[], const double *Vspec)
{
    unsigned char dip[NMAX+1], in;
    int n, m, k, pos[NMAX+1];
	double *V;

	#ifdef DEBUG
	checkstate(di, "summation_V1");
	#endif
    
    /* Pre-allocate with zeros. */
	memset(SV, 0, (NMAX+1)*sizeof(double));
	memset(SH, 0, (NMAX+1)*sizeof(double));

	/* Change 'd' and 'i_n' in (d,i1,i2,...,iN) and lookup value function.
     * 'n' starts at 0 to include common state 'd'. */
    for (n=0; n<=N; n++) {
        /* Only states with non-zero rates (in=nz[n][m])>0. 
         * nz[0] stops after D+1 indices since there's a zero in D+2. */
        for (m=1; (in=nz[n][m])>0; m++) {                  
			memcpy(dip, di, (NMAX+1)*sizeof(unsigned char));
			change(dip, n, in, pos);
			V = lookup_V1(encode(dip));

			/* Build up sums for incumbent firms. */
			for (k=1; (k<=N) && (di[k]<M+1); k++) {
				if ((n==k) && (in==M+1))                    /* Exit. */
					SV[k] += Vspec[k] * rates[n][in];       
                else                                        /* Continue. */
					SV[k] += V[pos[k]] * rates[n][in];
				SH[k] += rates[n][in];                      /* Sum hazards. */
			}
        }
    }

	/* Sums for potential entrants (short-lived). */
    /* Loop backwards since entrants are the last in the array (it's sorted!). */
	for (n=N; (n>=1) && (di[n]==M+1); n--)
        for (m=1; (in=nz[n][m])>0; m++) {
			SV[n] += Vspec[n] * rates[n][in];
			SH[n] += rates[n][in];
		}
}


/* ------------------------------------------------------------------------
 * Encoding, decoding, and sorting industry states. 
 * ----------------------------------------------------------------------*/

/* encode.
 * Encode state (d,i1,i2,...,iN). 
 * 'd' is common state (index); 'in' are firm-specific states. */
unsigned long encode(const unsigned char *di)
{
    unsigned long code, j = 0;
    int n;

	#ifdef DEBUG
	checkstate(di, "encode");
	#endif

	/* Map (i1,i2,...,iN) into j=0,...,J-1. */
	for (n=N; n>0; n--)
	    j += *(binom+di[n]+n-2+(di[n]-1)*(N+M+1));
	
	/* Compute linear index from subscript (d,j) in a Matlab array of size (D,J). */
	code = (di[0]-1)+D*j;

	return(code);
}


/* decode.
 * Decode state j into (d,i1,i2,...,iN). */
void decode(unsigned char *di, const unsigned long code)
{
    unsigned long j;
    int n;

	#ifdef DEBUG
	if (!((code >= 0) && (code < S)))
        mexErrMsgTxt("Undefined state (decode)!");
	#endif

	/* Compute subscript (d,j) in a Matlab array of size (D,J) from linear index. */
	memset(di, 1, (NMAX+1)*sizeof(unsigned char));
	j = code/D; 
	di[0] = code-j*D+1;

	/* Map j=0,...,J-1 into (i1,i2,...,iN). */
	for (n=N; n>0; n--) {
		while (*(binom+di[n]+n-1+di[n]*(N+M+1)) <= j)
			di[n]++;
		j -= *(binom+di[n]+n-2+(di[n]-1)*(N+M+1));
	}

	return;
}


/* sort.
 * Sort state (d,i1,i2,...,iN) such that i1<=...<=iN. */

/* Bubble sort. */
#ifdef BUBBLE
void sort(unsigned char *di, int *pos)
{
	int k, n, ind[NMAX+1], indn;
	unsigned char in;

	for (n=0; n<N+1; n++)
		ind[n] = n;

    for (k=2; k<N+1; k++)
        for (n=k; n>1; n--)
            if (di[n-1]>di[n]) {
                in = di[n]; indn = ind[n];
                di[n] = di[n-1]; ind[n] = ind[n-1];
                di[n-1] = in; ind[n-1] = indn;
            }

	for (n=0; n<N+1; n++)
		pos[ind[n]] = n;
}
/* Insertion sort. */
#else
void sort(unsigned char *di, int *pos)
{
	int k, n, ind[NMAX+1], indk;
	unsigned char ik;

	for (n=0; n<N+1; n++)
		ind[n] = n;

	for (k=2; k<N+1; k++) {
		ik = di[k]; indk = ind[k];
		for (n=k; (n>1) && (di[n-1]>ik); n--) {
			di[n] = di[n-1]; ind[n] = ind[n-1];
		}
		di[n] = ik; ind[n] = indk;
	}

	for (n=0; n<N+1; n++)
		pos[ind[n]] = n;
}
#endif


/* change.
 * Change nth coordinate of state (d,i1,i2,...,iN) to 'in' such that i1<=...<=iN. */
void change(unsigned char *di, const int n, const unsigned char in, int *pos)
{
	int k, ind[NMAX+1];

	#ifdef DEBUG
	checkstate(di, "change");
	#endif

	for (k=0; k<N+1; k++)
		ind[k] = k;

	if (n<1)
		di[n] = in;
	else {
		if (in>di[n]) {
			for (k=n; (k<N) && (in>di[k+1]); k++) {
				di[k] = di[k+1]; ind[k] = ind[k+1];
			}
			di[k] = in; ind[k] = n;
		}
		if (in<di[n]) {
			for (k=n; (k>1) && (in<di[k-1]); k--) {
				di[k] = di[k-1]; ind[k] = ind[k-1];
			}
			di[k] = in; ind[k] = n;
		}
	}

	for (k=0; k<N+1; k++)
		pos[ind[k]] = k;
}


/* checkstate.
 * Check if state is sorted and within range. */
void checkstate(const unsigned char *di, const char *caller)
{
	int n;

	for (n=1; n<N; n++)
		if (!(di[n] <= di[n+1]))
	        mexErrMsgIdAndTxt("", "State (d,i1,...,iN) must be sorted (%s)!", caller);
	if (!((di[0] >= 1) && (di[0] <= D)))
        mexErrMsgIdAndTxt("", "State (d,i1,...,iN) must be within range (%s)!", caller);
	for (n=1; n<=N; n++)
		if (!((di[n] >= 1) && (di[n] <= M+1)))
	        mexErrMsgIdAndTxt("", "State (d,i1,...,iN) must be within range (%s)!", caller);
}


/* checkparm.
 * Check if parameter is within range. */
double checkparm(double val, double from, double to, const char *parm, const char *caller)
{   
    double temp;
	if (!((val >= from) && (val <= to)))
        mexErrMsgIdAndTxt("", "%s must be between %g and %g (%s)!", parm, from, to, caller);
    temp = MAX(val,from);
	return (MIN(temp,to));
}



/* ------------------------------------------------------------------------
 * Helper functions. 
 * ----------------------------------------------------------------------*/

/* match.
 * Search a list for the first item to match a key (ignoring case). Return 
 * the index of the matching item or -1 if none exists. The list must be 
 * terminated by a pointer to an empty string. 
 * 'strnicmp' in original code replaced by 'strncasecmp' since the former is for Windows only. */
int match(const char *list[], const char *key)
{
    int len, ind;
    
    len = strlen(key);

    for (ind=0; list[ind]; ind++)
        if (strncasecmp(list[ind], key, len) == 0)          
            return ind;
    
    return -1;
}

/* read.
 * Read a matlab variable (i.e., an mxArray) from a MAT-file and return a pointer to it. 
 * Careful with const qualifiers and what types mx functions return. 
 * E.g.: *from and *to need to be const because mxGetDimensions 
 * returns a const; otherwise there will be a warning and potential error. */
mxArray *read(MATFile *mfp, const char *var_name)
{
    mxArray *array_ptr;
    const int *from, *to;
    
    if ((array_ptr = matGetVariable(mfp, var_name)) == NULL)
        mexErrMsgIdAndTxt("", "Cannot read matlab variable %s.", var_name); 

    if ((mxGetNumberOfElements(array_ptr) == 1) && mxIsNumeric(array_ptr)) {
        /* mexPrintf("Reading matlab variable %s=%g.\n", var_name, mxGetScalar(array_ptr)); */
        mxGetScalar(array_ptr);
    }
    else {        
        /* mexPrintf("Reading matlab variable %s (", var_name);*/
        from = mxGetDimensions(array_ptr);
        to = from+mxGetNumberOfDimensions(array_ptr)-1;
        /* for (; from<to; from++)
            mexPrintf("%dx", *from);
        mexPrintf("%d).\n", *from);
         */
    }

    return array_ptr;
}

/* write.
 * Write a matlab variable (i.e., an mxArray) to a MAT-file.
 * Careful with const qualifiers and what mx functions return. 
 * E.g.: *from and *to need to be const because mxGetDimensions 
 * returns a const; otherwise there will be a warning and potential error. */
void write(MATFile *mfp, const char *var_name, const mxArray *array_ptr)
{
    const int *from, *to;               

    if (matPutVariable(mfp, var_name, array_ptr) != 0)
        mexErrMsgIdAndTxt("", "Cannot write matlab variable %s.", var_name); 

    if ((mxGetNumberOfElements(array_ptr) == 1) && mxIsNumeric(array_ptr)) {
        /*mexPrintf("Writing matlab variable %s=%f.\n", var_name, mxGetScalar(array_ptr));*/
        mxGetScalar(array_ptr);
    }
    else {        
        /* mexPrintf("Writing matlab variable %s (", var_name); */
        from = mxGetDimensions(array_ptr);
        to = from+mxGetNumberOfDimensions(array_ptr)-1;
        /* for (; from<to; from++)
            mexPrintf("%dx", *from);
        mexPrintf("%d).\n", *from);
         */
    }

    return;
}

/* create.
 * Allocate dynamic memory to a data structure and return a pointer to it. */
double *create(void)
{
    double *data_ptr;
    
    if ((data_ptr = mxCalloc(N*S, sizeof(double))) == NULL)
        mexErrMsgTxt("Cannot create data structure.");

    return data_ptr;
}

/* mxFree.
 * Free dynamic memory allocated to a data structure. */
void destroy(double *data_ptr)
{
    mxFree(data_ptr);
}

/* copy.
 * Copy a data structure. */
void copy(double *to_data_ptr, const double *from_data_ptr)
{
    memcpy(to_data_ptr, from_data_ptr, N*S*sizeof(double));
}



/**************************************************************************
 * EOF *
 *************************************************************************/

/* ------------------------------------------------------------------------
 * Continuous-time dynamic stochastic game with R&D and learning by doing.
 *
 * Check first-order conditions and second order (existence) conditions for
 * a given converged solution. Also, examine non decreasing/increasing 
 * nature of the value function in own/competitor's productivity. 
 * 
 * Written by Bernardo Diaz de Astarloa, May 2015 @ Penn State. 
 *
 * --------------------------------------------------------------------- */


#include "mex.h"
#include "matrix.h"
#include "mat.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdarg.h>


/* Compiler directives. To be modified by user. */
#define NMAX 10                             /* Max. # of firms allowed in the market. */
#define DMAX 20                             /* Max. # of commons states allowed. */

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

double alpha;                               /* Demand: Price coefficient. */
double alpha1;                              /* Transition rates: Incumbent productivity. */
double alpha2;                              /* Transition rates: Incumbent productivity. */
double eta1;                                /* Transition rates: Incumbent lbd. */
double eta2;                                /* Transition rates: Incumbent lbd. */

double *V;                                  /* Inputs: Value function. */
double *p;                                  /* Inputs: Investment policy. */

mxArray *foc_ptr; double *foc;              /* Output: foc check. */
mxArray *ex_ptr; double *ex;                /* Output: 2nd order check (existence). */
mxArray *dvown_ptr; double *dvown;          /* Output: non decreasing/increasing property. */
mxArray *dvoth_ptr; double *dvoth;          /* Output: non decreasing/increasing property. */

/* Function prototypes. */
void init(const mxArray *prhs[]);
void checks(void);

/* Hazard rates functional form definitions. */
double hazlbd(const double);
double hazlbdprime(const double);
double hazlbdprime2(const double);

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
#define	lookup_V(s) (V+N*(s)-1)
#define	lookup_p(s) (p+N*(s)-1)

/* Need to define 'max' and 'min' which are not defined in C. 
 * Comment this out if in Windows (maybe change to lower case too?). */
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)


/********************************************************************************
 * Main routine.
 ********************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
	/* Check for proper number (and type) of arguments. */	
    if (nrhs != 3)
	    mexErrMsgIdAndTxt( "MATLAB:compMPE:invalidNumInputs","Three inputs required.");
    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:compMPE:inputNotStruct","First input must be a structure.");  
    
    /* Initialize. */
    init(prhs);

    /* Check conditions. */
    checks();
    
    /* Allocate output.
     * This is better than writing to a file since it won't create problems
     * when parallelizing. */  
    plhs[0] = mxDuplicateArray(foc_ptr); 
    plhs[1] = mxDuplicateArray(ex_ptr);
    plhs[2] = mxDuplicateArray(dvown_ptr);
    plhs[3] = mxDuplicateArray(dvoth_ptr);
    
	return;
}


/* Initialize. 
 * Reads and assigns parameter values and structures from MATLAB. */
void init(const mxArray *prhs[])
{    
    /* Get parameters and objects from input structure. 
     * In the future include warnings for potential errors when reading 
     * the structure (empty, types, etc.). */
    
    /* Scalars. */
    D = (unsigned char) checkparm(*mxGetPr(mxGetField(prhs[0],0,"D")), 1, DMAX, "D", "init");
    N = (int) checkparm(*mxGetPr(mxGetField(prhs[0],0,"N")), 2, NMAX, "N", "init");
    M = (unsigned char) checkparm(*mxGetPr(mxGetField(prhs[0],0,"M")), 1, UCHAR_MAX-1, "M", "init");
	S = (unsigned long) checkparm(*mxGetPr(mxGetField(prhs[0],0,"S")), 1, ULONG_MAX, "S", "init");        
    alpha = *mxGetPr(mxGetField(prhs[0],0,"alpha"));
    eta1 = *mxGetPr(mxGetField(prhs[0],0,"eta1"));
    eta2 = *mxGetPr(mxGetField(prhs[0],0,"eta2"));
    dsize = (unsigned char) *mxGetPr(mxGetField(prhs[0],0,"dsize"));
    p0size = (unsigned char) *mxGetPr(mxGetField(prhs[0],0,"p0size"));

    /* Arrays (matrices). */
    binom_ptr = mxGetField(prhs[0],0,"binom"); binom = (unsigned long*) mxGetData(binom_ptr);
    state_ptr = mxGetField(prhs[0],0,"state"); state = (unsigned char*) mxGetData(state_ptr);
    
    mgc_ptr = mxGetField(prhs[0],0,"mgc"); mgc = mxGetPr(mgc_ptr);
    pfor_ptr = mxGetField(prhs[0],0,"pfor"); pfor = mxGetPr(pfor_ptr);
    mkt_ptr = mxGetField(prhs[0],0,"mkt"); mkt = mxGetPr(mkt_ptr);

    /* Get pointers to value function and prices. */
    V = mxGetPr(prhs[1]); 
    p = mxGetPr(prhs[2]); 

    /* Create output matrices. */ 
    foc_ptr = mxCreateDoubleMatrix(N,S,mxREAL); foc = mxGetPr(foc_ptr);
    ex_ptr = mxCreateDoubleMatrix(N,S,mxREAL); ex = mxGetPr(ex_ptr);
    dvown_ptr = mxCreateDoubleMatrix(N,S,mxREAL); dvown = mxGetPr(dvown_ptr);
    dvoth_ptr = mxCreateDoubleMatrix(N,S,mxREAL); dvoth = mxGetPr(dvoth_ptr);

	return;
}


/* Perform checks. 
 * Checks FOCs, existence condition, and non-drecreasing/increasing 
 * property of value function. */
void checks(void)
{
    
	int n, k, st, iter = 1, pos[NMAX+1];
	unsigned long s, cntr;
	double diff, alldiff[NMAX+1], *Vup;
    double *Vs, *ps, temp1, temp2;
   
    unsigned char *di, *agg[D+1], dip[NMAX+1];
    unsigned char i, j, ind;
    double p0state, dstate, cn;
    
    double d, hsV, q[NMAX+1], sh[NMAX+1], hq[NMAX+1], hqq[NMAX+1], exist;

    
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
    
    /* Start iteration over states. */
    for (cntr=0; s=cntr, cntr<S; cntr++) {
        /* Preallocate arrays with zeros. */
        memset(alldiff, 0, (NMAX+1)*sizeof(double));

        /* Obtain unique state. */
        /* 'di' is N+1 vector (d,i1,...,iN). */
        /* dstate is demand level associated to d. */
        /* p0state is import prices level associated to d. */
        di = lookup_state(s);

        dstate = mkt[agg[di[0]][0]-1];
        p0state = pfor[agg[di[0]][1]-1];

        /* Lookup value and price for state s. */
        /* These are pointers to the mxArrays. */
        Vs = lookup_V(s);
        ps = lookup_p(s);
        
    
        /* Compute market shares and hazard derivatives. */
        d = 0.0;
        /* Compute denominator of mkt shares. */
        for (i=1; i<=N && di[i]<M+1; i++) 
            d += exp((-alpha)*(ps[i] - p0state));
        /* Compute mkt shares, and hazard derivatives. */
        for (i=1; i<=N && di[i]<M+1; i++) { 
            sh[i] = exp((-alpha)*(ps[i] - p0state))/(1+d);
            q[i] = sh[i]*dstate;
            hq[i] = hazlbdprime(q[i]);
            hqq[i] = hazlbdprime2(q[i]);
        }

        /* Loop over firms. */
        for (n=1; n<=N; n++) {
            if (di[n] < M+1) {
                
                /* Marginal cost. */
                cn = mgc[di[n]-1];

                /* deltaVs = Vup - V for incumbent n. */
                for (k=1; k<=N && di[k]<M; k++) {
                    memcpy(dip, di, (NMAX+1)*sizeof(unsigned char));
                    change(dip, k, di[k]+1, pos);           
                    Vup = lookup_V(encode(dip));
                    alldiff[k] = Vup[pos[n]] - Vs[n];
                }
                
                /* Sum of deltaVs term for FOCs. */
                hsV = 0.0;
                for (i=1; i<=N && di[i]<M+1; i++) { 
                    hsV += hq[i]*(sh[i])*(alldiff[i]);
                }
                
                /* FOC. */ 
                double f = (1/alpha) - (1-sh[n])*(ps[n]-cn) - hq[n]*(alldiff[n]) + hsV;
                
                /* Existence condition (second order). */
                double hqqn = hqq[n];
                hsV = 0.0;
                hqq[n] = 0;
                for (i=1; i<=N && di[i]<M+1; i++)
                    hsV += hqq[i]*(sh[i])*(sh[i])*(alldiff[i]);
                
                exist = (1/(alpha*q[n])) - hqqn*(1-sh[n])*(1-sh[n])*(alldiff[n]) - hsV;
                
                /* Fill checks. */
                foc[N*(s)+n-1] = f;
                ex[N*(s)+n-1] = exist;
                dvown[N*(s)+n-1] = alldiff[n];
                hsV = 0.0;
                for (i=1; i<=N; i++) { 
                    if (i!=n && alldiff[i]>0) {
                        hsV += 1;
                    }
                }
                dvoth[N*(s)+n-1] = (hsV>0);
                
            } else {
                
                foc[N*(s)+n-1] = 0;
                ex[N*(s)+n-1] = 2;
                dvown[N*(s)+n-1] = 2;
                dvoth[N*(s)+n-1] = -2;
            }
        }   /* for (n=1; n<=N; n++) */
    }   /* for (cntr=0; s=cntr, cntr<S; cntr++) */
    
    /* Free memory of structures. */
    for (i=1; i<=D; i++) {
        mxFree(agg[i]);
    }
    
    return;
}


/* ------------------------------------------------------------------------
 * Hazards and derivatives. 
 * ----------------------------------------------------------------------*/


/* hazlbd.
 * Compute hazard of LBD-productivity jump. */
double hazlbd(const double q)
{
    /* double h = eta1*q/(10+eta1*q); */
    double h = eta1*pow(q,eta2);
    return h;
}
/* hazlbdprime.
 * Compute derivative of hazard of LBD-productivity jump. */
double hazlbdprime(const double q)
{
    /* double hq = eta1/pow(10+eta1*q,2); */
    double hq = eta2*eta1*pow(q,eta2-1);
    return hq;
}
/* hazlbdprime2.
 * Compute second derivative of hazard of LBD-productivity jump. */
double hazlbdprime2(const double q)
{
    double hqq = (eta2-1)*eta2*eta1*pow(q,eta2-2);
    return hqq;
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


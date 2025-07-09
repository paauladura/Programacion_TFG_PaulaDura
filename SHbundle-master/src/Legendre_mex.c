/*
 * compile by:
 *    mex Legendre_mex.c legendre_od.c x_numbers.c
 */

#include <string.h> // for "strcmp"
#include "ctype.h"  // for "tolower"
#include <math.h>   // <-- absolutely necessary!!
#include "mex.h"

#include "legendre.h"
#include "x_numbers.h"

//#define DEBUG

// behaviour
#define NONE           0 // nothing, catch errors
#define PLM            1 // copies behaviour of Nico's PLM (nearly) with two paramaters 
#define PLM3           3 // copies behaviour of Nico's PLM (nearly) with three paramaters
// speed-optimized, calculates all P_lm up to l = m = Lmax;
// output in order/degree-sorting (P_00, P10, P20, ..., P_11, P21, ..., P_22, ...)#define SPEED_OD       4
#define SPEED_OD       4 //
// w/o X-numbers?
#define AUTO          32 // automatic determination: use Fukushima's X-numbers above Lmax >= 1500 (slower but stable)?
#define XNUM          64 // force X-numbers
#define ERROR       1024 // if wrong parameter was used

// ID values for mex-parameter-classes
#define ID_CHAR        4 
#define ID_DOUBLE      6 

// some constants
#define METHOD_STRLEN 16 // including \0 as end of string, i.e. 15 chars effectively

// X-numbers
#define XNUMDEG 1500  // the function switches to X-numbers for Lmax >= XNUMDEG

void lower_string(char s[])
{
  int c = 0;
  while (s[c]) {
    if (s[c] >= 'A' && s[c] <= 'Z') 
      s[c] |= 0x20;
    c++;
  }
}

//// MAIN FUNCTION ////
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
  int behaviour = NONE | AUTO, // defaults
    plm = 0, strparams = 0;
  char method[METHOD_STRLEN], *p;
  int i;

  int Lmax, order, degreelen, degreecnt, thetalen, numlp, matsiz1, matsiz2;
  double *degree, *theta;
  double *nL = NULL, *dnL = NULL, *ddnL = NULL;
  double *Wmm, *Wlm_1, *Wlm_2;

  // check input parameters
  if (nrhs < 2) mexErrMsgTxt("Not enough input arguments! Check the help text!"); // 2 parameters are minimum!

  // start searching for string parameters at 3rd parameter (i.e. position 2)
  for (i = 2; i < nrhs; i++) {
    if (mxGetClassID(prhs[i]) == ID_CHAR) { // do we have a string-parameter here?
      if (!(mxGetString(prhs[i], method, METHOD_STRLEN - 1))) // get first parameter string
        lower_string(method); // only lowercase if mxGetString succeeded
      else
        method[0] = '\0';
      strparams++; // one additional string parameter

      if (strcmp(method, "plm") == 0) {
        switch (i) {
          case 2:
            behaviour |= PLM;
            break;
          case 3:
            behaviour |= PLM3;
            break;
          default:
            behaviour = ERROR;
        }
      }	else if (strcmp(method, "speed_od") == 0) {
        if (i == 2) {
          behaviour |= SPEED_OD;
        } else {
          behaviour = ERROR;
        }
      } else if (strcmp(method, "xnum") == 0) { // forced (i.e. switch off auto, switch on X-numbers)
        behaviour &= ~AUTO;
        behaviour |= XNUM;
      } else if (strcmp(method, "noxnum") == 0) { // forced (i.e. switch off both auto and X-numbers)
        behaviour &= ~(AUTO | XNUM);
      }
    }
  }

  if (!(behaviour & (PLM | PLM3 | SPEED_OD)) || (behaviour & (PLM | PLM3))) { // no mode selected (assume PLM), or PLM/PLM3 selected
    switch (nrhs - strparams) {
    case 2: // check if params are all double --> behaviour: PLM(l_vec, theta_vec)
      if ((mxGetClassID(prhs[0]) == ID_DOUBLE) & (mxGetClassID(prhs[1]) == ID_DOUBLE))
	behaviour |= PLM; // assuming m = 0
      else
	behaviour = ERROR; // wrong kind of parameters!
      break;
    case 3: // check if params are all double --> behaviour: PLM(l_vec, m, theta_vec)
      if ((mxGetClassID(prhs[0]) == ID_DOUBLE) & (mxGetClassID(prhs[1]) == ID_DOUBLE) & (mxGetClassID(prhs[2]) == ID_DOUBLE)) 
	behaviour |= PLM3;
      else
	behaviour = ERROR; // wrong kind of parameters!
      break;
    }
  } else if (behaviour & SPEED_OD) {
    if (nrhs == 2) {
      if ((mxGetClassID(prhs[0]) == ID_DOUBLE) & (mxGetClassID(prhs[1]) == ID_DOUBLE))
	behaviour |= SPEED_OD; // assuming m = 0
      else
	behaviour = ERROR; // wrong kind of parameters!
    }
  } else
    behaviour = ERROR;
      
  if (behaviour & ERROR)
    mexErrMsgTxt("Something is wrong with your parameters! Check the help text!");

  // set up computation method and read parameters
  switch (behaviour & (PLM | PLM3 | SPEED_OD)) {
  case PLM:
    degree    = mxGetPr(prhs[0]);
    degreelen = (int)mxGetN(prhs[0]); // only length of lying vector!
    order     = 0;                    // order assumed to be 0
    theta     = mxGetPr(prhs[1]);
    thetalen  = (int)mxGetN(prhs[1]); // only length of lying vector!
    // searches for maximum > 0 in vector l, returns the maximum as integer
    Lmax = 0;
    for (i = 0; i < degreelen; i++)
      if (Lmax < (int)degree[i]) Lmax = (int)degree[i];
    break;

  case PLM3:
    degree    = mxGetPr(prhs[0]);
    degreelen = ((int)mxGetN(prhs[0]) > (int)mxGetM(prhs[0]) ? (int)mxGetN(prhs[0]) : (int)mxGetM(prhs[0])); // take bigger length as vector length! 
    order     = (int)mxGetScalar(prhs[1]);
    theta     = mxGetPr(prhs[2]);
    thetalen  = ((int)mxGetN(prhs[2]) > (int)mxGetM(prhs[2]) ? (int)mxGetN(prhs[2]) : (int)mxGetM(prhs[2])); // take bigger length as vector length! 
    // searches for maximum > 0 in vector l, returns the maximum as integer
    Lmax = 0;
    for (i = 0; i < degreelen; i++)
      if (Lmax < (int)degree[i]) Lmax = (int)degree[i];
    if (order > Lmax) // if order too big --> cut it down
      order = Lmax;
    if (order < 0) // disallow negative order
      order = 0;
    break;
    
  case SPEED_OD:
    Lmax     = (int)mxGetScalar(prhs[0]); 
    order    = Lmax;
    theta    = mxGetPr(prhs[1]);
    thetalen = (int)mxGetN(prhs[1]); // only length of lying vector!
    break;
    
  default:
    mexErrMsgTxt("Some unknown calculation method made it through to here! Please tell the developers!");
  }

  // switch on high degree/order for Lmax >= XNUMDEG
  if ((behaviour & AUTO) && (Lmax >= XNUMDEG)) {
    behaviour |= XNUM;
#ifdef DEBUG
    printf("INFO: Lmax = %d >= %d --> automatically switched to X-numbers\n", Lmax, XNUMDEG);
#endif
  }

  // set up output variables, initialize some variables to proper values
  if (nlhs > 3)
    mexErrMsgTxt("Wrong number of output arguments! Check the help text!");

  switch (behaviour & (PLM | PLM3 | SPEED_OD)) {
  case PLM:
  case PLM3:
    numlp  = (Lmax - order + 1);
    matsiz1 = thetalen;
    matsiz2 = degreelen;
    break;
  case SPEED_OD:
    matsiz1 = numlp = (int)(((double)Lmax + 1.0) * ((double)Lmax / 2.0 + 1.0));
    matsiz2 = thetalen;
    break;
  }

  plm = behaviour & PLM; // are we in one of the PLM modes?
  
#ifdef DEBUG
  if (plm) printf("PLM mode\n");
  else printf("SPEED_OD mode\n");
#endif

  plhs[0]  = mxCreateDoubleMatrix(matsiz1, matsiz2, mxREAL);
  nL       = mxGetPr(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(matsiz1, matsiz2, mxREAL);
    dnL     = mxGetPr(plhs[1]);
    if (nlhs > 2) { 
      plhs[2] = mxCreateDoubleMatrix(matsiz1, matsiz2, mxREAL);
      ddnL    = mxGetPr(plhs[2]);
    }
  }

  // allocate memory
  if (plm) {
    Wmm    = (double*)malloc(sizeof(double) * (order + 1));
    Wlm_1  = (double*)malloc(sizeof(double) * (Lmax - order + 1));
    Wlm_2  = (double*)malloc(sizeof(double) * (Lmax - order));
  } else {
    Wmm    = (double*)malloc(sizeof(double) * Lmax);
    Wlm_1  = (double*)malloc(sizeof(double) * Lmax * (Lmax + 1) / 2);
    Wlm_2  = (double*)malloc(sizeof(double) * Lmax * (Lmax - 1) / 2);
  }

  // prefactors for normalisation
  if (plm) {
    init_legendre_1st_kind_plm(Lmax, order, Wmm, Wlm_1, Wlm_2);
  } else {
    init_legendre_1st_kind_od_speed(Lmax, Wmm, Wlm_1, Wlm_2);
  }
  
  // calculation of fully normalised associated Legendre functions
  if ((behaviour & XNUM) == XNUM) { // SWITCH TO X-NUMBERS!!
    if (plm) {
      legendre_1st_kind_xnum_plm(Lmax, order, numlp, degreelen, degree, thetalen, theta, nL, dnL, ddnL, Wmm, Wlm_1, Wlm_2);
    } else {
      legendre_1st_kind_od_xnum_speed(Lmax, thetalen, theta, nL, dnL, ddnL, Wmm, Wlm_1, Wlm_2);
    }
  } else { // USE NORMAL (UNSTABLE) WAY
    if (plm) 
      legendre_1st_kind_plm(Lmax, order, numlp, degreelen, degree, thetalen, theta, nL, dnL, ddnL, Wmm, Wlm_1, Wlm_2);
    else // all degrees/orders (no PLM behaviour, no X-numbers)
      legendre_1st_kind_od_speed(Lmax, thetalen, theta, nL, dnL, ddnL, Wmm, Wlm_1, Wlm_2);
  }
 
  // don't forget to free the allocated memory:
  free(Wmm); free(Wlm_1); free(Wlm_2);
}


#include <R.h>

static double parms[11];

/* define parameters as macros */

/* AdpA */
#define beta_AdpA      parms[0]
#define gamma_AdpA     parms[1]
#define k1_AdpA        parms[2]
#define k2_AdpA        parms[3]
#define sigma_AdpA     parms[4]
#define n1             parms[5]
#define n2             parms[6]

/* BldA */
#define gamma_BldA      parms[7]
#define k1_BldA         parms[8]
#define sigma_BldA      parms[9]
#define p               parms[10]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
    int N=11;
    odeparms(&N, parms);
}

/* derivatives */
void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
    // if (ip[0] < 0) error("nout should be zero");

    ydot[0] = (beta_AdpA * pow(y[0],n1)) / (pow(k1_AdpA,n1) + pow(y[0],n1)) + (gamma_AdpA * pow(y[1],n2)) / (pow(k2_AdpA,n2) + pow(y[1],n2)) - sigma_AdpA * y[0];
    ydot[1] = (gamma_BldA * pow(y[0],p)) / (pow(k1_BldA,p) + pow(y[0],p)) - sigma_BldA * y[1];

    yout[0] = ydot[0];
    yout[1] = ydot[1];
}

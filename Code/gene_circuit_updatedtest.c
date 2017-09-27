#include <R.h>

static double parms[10];

/* define parameters as macros */

/* AdpA */
#define beta_AdpA      parms[0]
#define k1_AdpA        parms[1]
#define k2_AdpA        parms[2]
#define sigma_AdpA     parms[3]

/* BldA */
#define gamma_BldA      parms[4]
#define k1_BldA         parms[5]
#define sigma_BldA      parms[6]
#define p               parms[7]
#define sigma_adpAchange parms[8]
#define sigma_bldAchange parms[9]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
    int N=10;
    odeparms(&N, parms);
}

void event(int *n, double *t, double *y) {
	sigma_AdpA= sigma_adpAchange;
	sigma_BldA= sigma_bldAchange;
}

/* derivatives */
void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
	double n1=1.2;
    double n2= 1.68;
    double gamma_AdpA= 320; 

    // if (ip[0] < 0) error("nout should be zero");
    ydot[0] = (beta_AdpA * pow(y[0],n1)) / (pow(k1_AdpA,n1) + pow(y[0],n1)) + (gamma_AdpA * pow(y[1],n2)) / (pow(k2_AdpA,n2) + pow(y[1],n2)) - sigma_AdpA * y[0];
    ydot[1] = (gamma_BldA * pow(y[0],p)) / (pow(k1_BldA,p) + pow(y[0],p)) - sigma_BldA * y[1];

    yout[0] = ydot[0];
    yout[1] = ydot[1];
}
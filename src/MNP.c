/******************************************************************
 This file contains Cpp functions for fitting symmetric multinomial probit models.
 It is adapted from the codebase the codebase of MNP (JSS 2005), now modified to run the algorithm preseneted in
 the sMNP paper of Burgette, Puelz, Hahn (2020)
 *******************************************************************/

// To do:
// load in W from the R side
// n_dim refers to the number of columns in Sigma_{-b}
// figure out prior for beta

#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cMNPgibbs(int *piNDim, int *piNCovCat, int *piNCovResp, int *piNSamp, int *piNGen, 
	       double *pdA0, int *piNu0, double *pdS, double *pdX, double *pdW,
			   int *bStart, // Initial b value
	       int *y,        /* response variable: -1 for missing */
	       double *pdbeta, double *pdSigma, int *piImp, 
	       int *invcdf,   /* use inverse cdf for TruncNorm? */
	       int *piBurnin, /* the number of burnin */
	       int *piKeep,
	       int *verbose,  /* 1 if extra print is needed */ 
	       int *latent, /* 1 if W is stored */
           int *corrected, /*record beta at correct step*/
	       double *pdStore){
  
  /* parameters from R */
    int n_samp = *piNSamp; /* sample size */
    int n_gen = *piNGen;   /* number of gibbs draws */
	int n_cov_resp_spec = *piNCovResp;
	int n_cov_cat_spec_par = *piNCovCat;
    int n_dim = *piNDim - 1;   /* The number of indentifiable dimension (p-1) */
	int n_dim_long = n_dim + 1;  // number of choices	
	int n_cov_short = n_cov_cat_spec_par * n_dim + n_cov_resp_spec;
	int n_cov_long = n_cov_cat_spec_par * n_dim_long + n_cov_resp_spec;
    int nu0 = *piNu0;      /* nu0: degree of freedom, S: scale matrix */
    int imp = *piImp;      /* 0: improper prior for beta - algorithms 1 & 2
			    1: proper prior for beta */
    int keep = 1;          /* counter for keep */
    int progress = 1;      /* counter for progress report */
    double *beta;	         /* The model parameters */
    double **Sigma;        /* The unidentifiable variance matrix */
    double **X;            /* The covariates and outcome var */
    double **A0;           /* prior precision for beta */
    double **S;            /* The prior parameters for Sigma */
	double *betaLong;

  /* model parameters */
    double **W;                /* The missing data! */
	double **XminusB;
    double cmean;  	     /* The conditional mean for Wij */
    double cvar;               /* The conditional variance for Wij */
    double **Xbeta;
    double **Xbeta_new;
    double ***PerSig;          /* Permuted Sigma for sampling W */
    int kindx, lindx;          /* Indexes for the permutation of Sigma */
    double *mbeta;             /* Posterior mean of beta */
    double **Vbeta;   	     /* The var of beta, given Yaug */
    double **Xstar;            /* L*X where L is the Chol factor of SigInv */
    double **SigInv;           /* The Inverse of Sigma */
    double *epsilon;           /* The error term */
    double *epsilon_new;           /* The error term */
	double *epsilonMinJ;     // with one dimension removed
    double **R;                /* Sum of squares matrix e'e */
    double **SS;	             /* Sum of squares for sweep operator */
    double alpha2 = 1.0;             /* alpha^2: Inv-chisq df=nu0 */
    double ss;                 /* scale parameter for alpha^2 */
	int *hldMult;            // Hold outcome of multinomial draw
    int i, j, k, l, main_loop; /* used for loops */
	int bCat=0; // faux base param

	double wSum;
	double currWMax;
	double sampWMax;
	double sampWMin;
	double *logDet;

  /* temporay storages */
    int itemp, itemp2, itempMax, itempMin, *ivtemp, itempS = 0, itempP=ftrunc((double) n_gen/10);
	double dtemp;
    double *vtemp;
    double **mtemp,**mtemp1,**mtemp2;
	double maxDet;

	Rprintf("%i dimensions.\n", n_dim);
	
    /** get random seed **/
    GetRNGstate();

    /** defining vectors and matricies **/
    W = doubleMatrix(n_samp, n_dim_long);
	
    X = doubleMatrix(n_samp*n_dim_long+n_cov_long, n_cov_long+1);
	XminusB = doubleMatrix(n_samp*n_dim_long + n_cov_long, n_cov_short + 1); //need long rows; go long then cut down
    Xbeta = doubleMatrix(n_samp, n_dim_long);
    Xbeta_new = doubleMatrix(n_samp, n_dim_long);
    SS = doubleMatrix(n_cov_short+1, n_cov_short+1);
    Sigma = doubleMatrix(n_dim, n_dim);
    SigInv = doubleMatrix(n_dim, n_dim);
    epsilon = doubleArray(n_samp*n_dim_long);
    epsilon_new = doubleArray(n_samp*n_dim_long);
	epsilonMinJ = doubleArray(n_samp * n_dim);
    R = doubleMatrix(n_dim, n_dim);
    beta = doubleArray(n_cov_short);
	betaLong = doubleArray(n_cov_long);
    mbeta = doubleArray(n_cov_short);
    Vbeta = doubleMatrix(n_cov_short, n_cov_short);
    A0 = doubleMatrix(n_cov_short, n_cov_short);
    S = doubleMatrix(n_dim, n_dim);
    vtemp = doubleArray(n_dim-1);
	hldMult = intArray(n_dim_long);
    ivtemp = intArray(n_dim_long);
    Xstar = doubleMatrix(n_samp*n_dim +n_cov_short, n_cov_short+1);
    mtemp = doubleMatrix(n_cov_short, n_cov_short);
    mtemp1 = doubleMatrix(n_dim, n_dim);
    mtemp2 = doubleMatrix(n_dim, n_dim);
    PerSig = doubleMatrix3D(n_dim, n_dim, n_dim);
	logDet = doubleArray(n_dim_long);

  /** Packing X, A0, S, beta, Sigma  **/
    itemp = 0;
        for (k = 0; k < n_cov_long; k++)
            for (j = 0; j < n_dim_long; j++)
                for (i = 0; i < n_samp; i++)
                    X[i*n_dim_long+j][k] = pdX[itemp++];
    
    itemp = 0;
        for (k = 0; k < n_cov_short; k++)
            for (j = 0; j < n_cov_short; j++)
                A0[j][k] = pdA0[itemp++];
	

    itemp = 0;
        for (k = 0; k < n_dim; k++)
            for (j = 0; j < n_dim; j++)
                S[j][k] = pdS[itemp++];
	
	PdoubleMatrix(S, n_dim,n_dim);

    itemp = 0;
        for (j=0; j<n_cov_short;j++)
            beta[j]=pdbeta[itemp++];
	
	for(j=0; j<n_cov_long;j++)
		betaLong[j] = 0.0;

    itemp = 0;
        for (k = 0; k < n_dim; k++)
            for (j = 0; j < n_dim; j++)
                Sigma[j][k] = pdSigma[itemp++];
	
  dinv(Sigma,n_dim,SigInv); 
	
  /* add prior information as additional data points
   here, this is where big B is stored*/
    if(imp){
        for(j=0;j<n_cov_short;j++)
            for(k=0;k<=n_cov_short;k++)
                Xstar[n_samp*n_dim+j][k]=0;
            }
        else{
            dcholdc(A0,n_cov_short,mtemp); /* Cholesky decomposition */
            for(j=0;j<n_cov_short;j++){
                Xstar[n_samp*n_dim+j][n_cov_short]= 0;
            for(k=0;k<n_cov_short;k++){
                Xstar[n_samp*n_dim+j][n_cov_short]+=mtemp[j][k]*0;
                Xstar[n_samp*n_dim+j][k]=mtemp[j][k];
                    }
                }
            }


	
	itemp = 0;
	for(k = 0; k < n_dim_long; k++)
		for(i=0; i<n_samp; i++){
			W[i][k] = pdW[itemp];
			itemp++;
		}
    /*** GIBBS SAMPLER! ***/
    for(main_loop=1; main_loop<=n_gen; main_loop++){
        for(j=0;j<n_dim;j++)
            for(k=0;k<n_dim;k++) mtemp1[j][k]=0;
        for(i=0;i<n_dim;i++)
            for(j=0;j<n_dim;j++)
                for(k=0;k<n_dim;k++)
                    mtemp1[j][k]+=S[j][i]*SigInv[i][k];
        
        
        // sigma --> sigmatilde
        for(j=0;j<n_dim;j++){
            for(k=0;k<n_dim;k++){
                Sigma[j][k]*=alpha2;
                SigInv[j][k]/=alpha2;
            }
        }
        
        
        // beta --> betatilde (puelz -- 12/10/2019)
        for(j=0;j<n_cov_short;j++)
            beta[j] *= sqrt(alpha2);
        for(j=0;j<n_cov_long;j++)
            betaLong[j] *= sqrt(alpha2);
        
        // ----------------------------------------------------------
        /** permutation of Sigma  **/
        
        
        for(j=0;j<n_dim;j++){
            kindx = 0;
            for(k=0;k<n_dim;k++){
                lindx=0;
                for(l=0;l<n_dim;l++){
                    if(j!=k)
                        if(j!=l)
                            PerSig[j][k+kindx][l+lindx]=Sigma[k][l];
                        else{
                            lindx=-1;
                            PerSig[j][k+kindx][n_dim-1]=Sigma[k][l];
                        }
                        else{
                            kindx=-1;
                            if(j==l){
                                lindx=-1;
                                PerSig[j][n_dim-1][n_dim-1]=Sigma[j][j];
                            }
                            else
                                PerSig[j][n_dim-1][l+lindx]=Sigma[k][l];
                        }
                }
            }
            dinv(PerSig[j],n_dim,PerSig[j]);
        }
        
        
        /** Truncated Multinomial Normal Sampling for W **/
        //  PintArray(y, 10);
        
        for(i=0;i<n_samp;i++){
            for(j=0;j < n_dim_long;j++) {
                Xbeta[i][j]=0;
                for(k=0;k<n_cov_long;k++) Xbeta[i][j]+=X[i*n_dim_long+j][k]*betaLong[k];
            }
            
            for(j=0;j<n_dim_long;j++){
                if(j != bCat){
                    wSum = 0.0;
                    for(k = 0; k<n_dim_long; k++)
                        wSum += W[i][k];
                    wSum = wSum - (W[i][bCat] + W[i][j]);
                    
                    if(y[i] == bCat){
                        currWMax = -INFINITY;
                        sampWMin = -10000;
                        for(k=0;k<n_dim_long; k++)
                            if((k != bCat) && (k != j))
                                currWMax = fmax2(currWMax, W[i][k]);
                        sampWMax = fmin2(-.5 * wSum, -1*(currWMax + wSum));
                    }
                    else if(y[i] != j){
                        itemp = y[i];
                        sampWMax = W[i][itemp];
                        sampWMin = -1*wSum - W[i][itemp];
                    }
                    else {
                        currWMax = -INFINITY;
                        for(k =0; k<n_dim_long; k++)
                            if ((k != bCat) && (k != j))
                                currWMax = fmax2(currWMax, W[i][k]);
                        sampWMin = fmax2(-.5*wSum, currWMax);
                        sampWMax = 10000;
                    }
                    
                    /* conditional mean and variance */
                    itemp=0;
                    for (k=0;k<n_dim_long;k++)
                        if((j!=k) && (k != bCat))
                            vtemp[itemp++]=W[i][k]-Xbeta[i][k];
                    
                    cmean=Xbeta[i][j];
                    itemp = j - (j > bCat);
                    cvar=1/PerSig[itemp][n_dim-1][n_dim-1];
                    for(k=0;k<(n_dim-1);k++)
                        cmean-=PerSig[itemp][n_dim-1][k]*vtemp[k]*cvar;
                    /* sampling each W[i][j] conditionally on the other elements */
                    
                    if(sampWMax < sampWMin){
                        Rprintf("%i bCat \n", bCat);
                        Rprintf("%i j\n", j);
                        Rprintf("%i y[i] \n", y[i]);
                        Rprintf("%f sampWMax\n", sampWMax);
                        Rprintf("%f sampWMix\n", sampWMin);
                        for(k=0; k<n_dim_long; k++)
                            Rprintf("%f \n", W[i][k]);
                        PdoubleArray(W[i], n_dim_long);
                    }
                    
                    if(y[i] == -1)
                        W[i][j]=cmean+norm_rand()*sqrt(cvar);
                    else
                        W[i][j]=TruncNorm(sampWMin, sampWMax, cmean, cvar,*invcdf);
                    wSum = 0.0;
                    for(k=0; k<n_dim_long;k++)
                        wSum += W[i][k];
                    wSum = wSum - W[i][bCat];
                    W[i][bCat] = -1*wSum;
                }
            }
        }
        
        for(i =0; i<n_samp; i++)
            for(j=0; j<n_dim_long; j++)
                X[i*n_dim_long+j][n_cov_long]=W[i][j];
        
        // Get rid of the columns of X that relate to bCat-th category:
        itemp = 0;
        for(k=0;k<(n_cov_cat_spec_par * n_dim_long); k++) // n_cov_cat_spec_par is the number of class-specific parameters
            if((k % n_dim_long) != bCat){
                for(i=0; i<(n_samp*n_dim_long + n_cov_long); i++){
                    XminusB[i][itemp]=X[i][k];
                }
                itemp++;
            }
        
        itemp = n_cov_cat_spec_par * n_dim;
        for(k=(n_cov_cat_spec_par * n_dim_long); k < (n_cov_long +1); k++){
            for(i=0; i<(n_samp*n_dim_long + n_cov_long); i++)
                XminusB[i][itemp] = X[i][k];
            itemp++;
        }
        
        // Get rid of rows that that correspond to bCat-th category
        itemp = 0;
        for(i=0; i<(n_dim_long * n_samp + n_cov_long); i++)
            if ((i % n_dim_long) != bCat){
                for(j = 0; j<(n_cov_short + 1); j++)
                    XminusB[itemp][j] = XminusB[i][j];
                itemp++;
            }
        
        
        // ----------------------------------------------------------
        // sigmatilde --> sigma
        for(j=0;j<n_dim;j++){
            for(k=0;k<n_dim;k++){
                Sigma[j][k]/=alpha2;
                SigInv[j][k]*=alpha2;
            }
        }
        
        // ----------------------------------------------------------
        /* construct matrix Xstar  .. multiply X and W
         by the Inverse of the Cholesky factor, need this
         to get XSigX and XSigW */
        dcholdc(SigInv,n_dim,mtemp1);
        for(i=0;i<n_samp*n_dim;i++)
            for(j=0;j<=n_cov_short;j++) Xstar[i][j]=0;
        for(i=0;i<n_samp;i++)
            for(j=0;j<n_dim;j++)
                for(k=0;k<n_dim;k++)
                    for(l=0;l<=n_cov_short;l++)
                        Xstar[i*n_dim+k][l]+=mtemp1[j][k]*XminusB[i*n_dim+j][l];
        
        /* construct SS matrix for SWEEP */
        for(j=0;j<=n_cov_short;j++)
            for(k=0;k<=n_cov_short;k++) SS[j][k]=0;
        for(i=0;i<n_samp;i++)
            for(j=0;j<n_dim;j++)
                for(k=0;k<=n_cov_short;k++)
                    for(l=0;l<=n_cov_short;l++)
                        SS[k][l]+=Xstar[i*n_dim+j][k]*Xstar[i*n_dim+j][l];
        
        /* SWEEP to get posterior mean and variance for beta
         the SWEEP operation is an iterative way compute the
         inverse of a matrix */
        for(j=0;j<n_cov_short;j++)
            SWP(SS,j,n_cov_short+1);
        
        /* draw beta given Sigma and W */
        for(j=0;j<n_cov_short;j++){
            mbeta[j]=SS[j][n_cov_short];
            for(k=0;k<n_cov_short;k++) Vbeta[j][k]=-SS[j][k]*alpha2;
        }
        
        rMVN(beta,mbeta,Vbeta,n_cov_short);
        
        // expand beta
        for(i = 0; i<n_cov_cat_spec_par; i++){
            dtemp = 0.0;
            for(j = 0; j<n_dim;j++)
                dtemp += beta[i*n_dim + j];
            for(j=0; j<n_dim_long;j++){
                if(j != bCat)
                    betaLong[i * n_dim_long + j] = beta[i * n_dim + j - (j > bCat)];
                else betaLong[i * n_dim_long + j] = -1 * dtemp;
                
            }
        }
        for(i = n_cov_cat_spec_par *n_dim; i<n_cov_short; i++)
            betaLong[i + n_cov_cat_spec_par] = beta[i];
        
        
        /* draw Sigma and bCat given beta and W */
        for(i=0;i<n_samp;i++)
            for(j=0;j<n_dim_long;j++){
                epsilon[i*n_dim_long+j]=X[i*n_dim_long+j][n_cov_long];
                for(k=0;k<n_cov_long;k++)
                    epsilon[i*n_dim_long+j]-=X[i*n_dim_long+j][k]*betaLong[k];
            }
        
        for(i = 0; i< n_dim_long; i++){ // calculate OP_b for each b
            itemp = 0;
            for(j = 0; j<(n_dim_long * n_samp); j++)
                if((j % n_dim_long) != i){
                    epsilonMinJ[itemp] = epsilon[j];
                    itemp++;
                }
            for(j=0;j<n_dim;j++)
                for(k=0;k<n_dim;k++)
                    R[j][k]=0;
            for(l=0;l<n_samp;l++)
                for(j=0;j<n_dim;j++)
                    for(k=0;k<n_dim;k++)
                        R[j][k]+=epsilonMinJ[l*n_dim+j]*epsilonMinJ[l*n_dim+k];
            for(j=0;j<n_dim;j++)
                for(k=0;k<n_dim;k++)
                    mtemp1[j][k]=S[j][k]+R[j][k];
            logDet[i] = (n_samp + nu0)/-2 * ddet(mtemp1, n_dim, 1);
        }
        
        maxDet = logDet[0];
        for(i = 1; i < n_dim_long; i++)
            maxDet = fmax2(maxDet, logDet[i]);
        dtemp = 0.0;
        for(i=0; i< n_dim_long; i++){
            logDet[i] = logDet[i] - maxDet;
            logDet[i] = exp(logDet[i]);
            dtemp += logDet[i];
        }
        for(i=0;i<n_dim_long; i++)
            logDet[i] = logDet[i]/dtemp;
        
        // sample bCat here
        rmultinom(1, logDet, n_dim_long, hldMult);
        
        bCat = 0;
        for(i=0; i<n_dim_long; i++)
            bCat += i*hldMult[i];
        
        itemp = 0;
        for(j = 0; j<(n_dim_long * n_samp); j++)
            if((j % n_dim_long) != bCat){
                epsilonMinJ[itemp] = epsilon[j];
                itemp++;
            }
        
        for(j=0;j<n_dim;j++)
            for(k=0;k<n_dim;k++)
                R[j][k]=0;
        for(l=0;l<n_samp;l++)
            for(j=0;j<n_dim;j++)
                for(k=0;k<n_dim;k++)
                    R[j][k]+=epsilonMinJ[l*n_dim+j]*epsilonMinJ[l*n_dim+k];
        for(j=0;j<n_dim;j++)
            for(k=0;k<n_dim;k++)
                mtemp1[j][k]=S[j][k]+R[j][k];
        
        dinv(mtemp1,n_dim,mtemp2);
        rWish(SigInv,mtemp2,nu0+n_samp,n_dim);
        dinv(SigInv,n_dim,Sigma);
        
        // update alpha2 here!! (and only here)
        alpha2=0;
        for(k=0;k<n_dim; k++)
            alpha2+=Sigma[k][k];
        alpha2 = alpha2/n_dim;
        // sigmatilde --> sigma
        for(j=0;j<n_dim;j++){
            for(k=0;k<n_dim;k++){
                Sigma[j][k]/=alpha2;
                SigInv[j][k]*=alpha2;
            }
        }
        
        // update beta
        for(j=0;j<n_cov_short;j++)
            beta[j] /= sqrt(alpha2);
        for(j=0;j<n_cov_long;j++)
            betaLong[j] /= sqrt(alpha2);
        
        
        // ISSUE #1 (puelz) is X=tilde{W}?
        for(i=0;i<n_samp;i++){
            for(j=0;j<n_dim_long;j++){
                W[i][j]=X[i*n_dim_long+j][n_cov_long]/sqrt(alpha2);
            }
        }
        
        
        /* print Gibbs draws for all the schmes! */
        R_CheckUserInterrupt();
        if(main_loop > *piBurnin) {
            if(keep==*piKeep) {
                for(j=0;j<n_cov_long;j++)
                    pdStore[itempS++]=betaLong[j];
                for(j=0;j<n_dim;j++)
                    for(k=0;k<n_dim;k++)
                        if(j<=k)
                            pdStore[itempS++]=Sigma[j][k];
                pdStore[itempS++]= (double)bCat;
                keep=1;
            }
            else
                keep++;
        }
        
        if(*verbose) {
            if(main_loop == itempP) {
                Rprintf("%3d percent done (corrected sMNP - 12-11-2019).\n", progress*10);
                itempP+=ftrunc((double) n_gen/10); progress++;
                R_FlushConsole();
            }
        }
        
        
    }
   /* end of Gibbs sampler */
  
  /** write out the random seed **/
  
    PutRNGstate();
  
  /** freeing memory **/
	  

    FreeMatrix(W, n_samp);
    FreeMatrix(X, n_samp*n_dim_long+n_cov_long);
    FreeMatrix(XminusB, n_samp*n_dim_long + n_cov_long);
    FreeMatrix(Xbeta, n_samp);
    FreeMatrix(SS, n_cov_short+1);
    FreeMatrix(Sigma, n_dim);
    FreeMatrix(SigInv, n_dim);
    free(epsilon);
    free(epsilonMinJ);
    FreeMatrix(R, n_dim);
    free(beta);
    free(betaLong);
    free(mbeta);
    FreeMatrix(Vbeta, n_cov_short);
    FreeMatrix(A0, n_cov_short);
    FreeMatrix(S, n_dim);
    free(vtemp);
    free(ivtemp);
    free(hldMult);
    FreeMatrix(Xstar, n_samp*n_dim +n_cov_short);
    FreeMatrix(mtemp, n_cov_short);
    FreeMatrix(mtemp1, n_dim);
    FreeMatrix(mtemp2, n_dim);
    Free3DMatrix(PerSig, n_dim, n_dim);
    free(logDet);
/* main */
}

/* unitility function */
void R_max_col2(double **matrix, 
                int nr,
                int nc,
                int *maxes,
                int ties_meth) {
    
    int *ncol = intArray(1);
    int *nrow = intArray(1);
    int *ties = intArray(1);
    int *itmp = intArray(1);
    double *tmp = doubleArray(nr*nc);
    int i, j, k;
    
    ncol[0] = nc; nrow[0] = nr; ties[0] = ties_meth;
    i = 0;
    for (k = 0; k < nc; k++)
        for (j = 0; j < nr; j++)
            tmp[i++] = matrix[j][k];
    
    R_max_col(tmp, nrow, ncol, maxes, ties);
    
    free(ncol);
    free(nrow);
    free(itmp);
    free(tmp);
    
}
 
void predict(double *dX,     /* X matrix */
             int *nobs,      /* number of observations */
             double *dcoef,  /* coefficients */
             double *dSigma, /* covariances */
             int *dBase,
             int *ndims,     /* number of dimensions */
             int *ncovs,     /* number of covariates */
             int *ndraws,    /* number of MCMC draws */
             int *moredraws, /* number of extra draws */
             int *verbose,
             double *prob,   /* probability output */
             double *choice, /* choice output */
             double *order  /* order output */
) {
    
    int n_samp = *nobs;
    int n_draw = *ndraws;
    int n_dim = *ndims;
    int n_cov = *ncovs;
    int n_extra = *moredraws;
    
    double **X = doubleMatrix(n_samp*n_dim, n_cov);
    double *Xbeta = doubleArray(n_cov);
    double *zeros =doubleArray(n_dim-1);
    double *vtemp = doubleArray(n_dim-1);
    double *vtemp2 = doubleArray(n_dim);
    int *baseSamp = intArray(n_draw);
    double **W = doubleMatrix(n_extra, n_dim);
    double **beta = doubleMatrix(n_draw, n_cov);
    double **mtemp = doubleMatrix(n_dim, n_dim);
    double ***Sigma = doubleMatrix3D(n_draw, n_dim-1, n_dim-1);
    
    int i, j, k, main_loop, itemp, itempP, itempO, itempC;
    int total = n_extra*n_samp*n_draw;
    int count, progress = 1;
    int itempQ = ftrunc((double) total/10);
    int *maxdim = intArray(n_extra);
    int *ind = intArray(n_dim);
    int *sumorder = intArray(n_dim);
    int *probTemp = intArray(n_dim);
    double dtemp;
    
    for(k=0; k<(n_dim -1); k++)
        zeros[k] = 0.0;
    
    /** reading the data */
    itemp = 0;
    for (k = 0; k < n_cov; k++)
        for (j = 0; j < n_dim; j++)
            for (i = 0; i < n_samp; i++)
                X[i*n_dim+j][k] = dX[itemp++];
    
    
    PdoubleMatrix(X, 10, n_cov);
    
    /** reading the MCMC draws **/
    itemp = 0;
    for (k = 0; k < n_cov; k++)
        for (j = 0; j < n_draw; j++)
            beta[j][k] = dcoef[itemp++];
    
    PdoubleMatrix(beta, 10, n_cov);
    
    itemp = 0;
    for (k = 0; k < n_dim-1; k++)
        for (j = k; j < n_dim-1; j++)
            for (i = 0; i < n_draw; i++) {
                Sigma[i][j][k] = dSigma[itemp++];
                Sigma[i][k][j] = Sigma[i][j][k];
            }
    
    itemp =0;
    for (k=0; k<n_draw; k++)
        baseSamp[k] =dBase[itemp++];
    
    
    /** get random seed **/
    GetRNGstate();
    
    /** Posterior predictive simulations **/
    itempC = 0; itempO = 0; itempP = 0; count = 0;
    itempQ = 0;
    
    /* loop for observations */
    for (i = 0; i < n_samp; i++) {
        if (n_extra == 1) {
            for (j = 0; j < n_dim; j++)
                probTemp[j] = 0;
        }
        /* loop for MCMC draws */
        for (main_loop = 0; main_loop < n_draw; main_loop++) {
            if (n_extra > 1) {
                for (j = 0; j < n_dim; j++)
                    probTemp[j] = 0;
            }
            /* compute the mean for each dimension */
            for (j = 0; j < n_dim; j++) {
                Xbeta[j] = 0;
                for (k = 0; k < n_cov; k++)
                    Xbeta[j] += X[i*n_dim+j][k] * beta[main_loop][k];
            }
            /* PdoubleArray(Xbeta, n_dim); */
            /* sample W */
            for (j = 0; j < n_extra; j++) {
                /*dinv(Sigma[main_loop], n_dim, mtemp);*/
                rMVN(vtemp, zeros, Sigma[main_loop], n_dim-1);
                for(k=0;k<n_dim; k++)
                    vtemp2[k] = 0;
                for(k=0;k<n_dim; k++)
                    if(k != baseSamp[main_loop])
                        vtemp2[k] = vtemp[k - (baseSamp[main_loop] < k)];
                dtemp = 0;
                for(k=0;k<n_dim; k++)
                    dtemp += vtemp2[k];
                vtemp2[baseSamp[main_loop]] = -1*dtemp;
                for (k = 0; k < n_dim; k++)
                    W[j][k] = vtemp2[k] + Xbeta[k];
                
            }
            /* which dimension is max for each of n_extra W draws? */
            R_max_col2(W, n_extra, n_dim, maxdim, 1);
            /* order */
            for (j = 0; j < n_extra; j++) {
                for (k = 0; k < n_dim; k++) {
                    ind[k] = k; sumorder[k] = 0;
                }
                revsort(W[j], ind, n_dim);
                for (k = 0; k < n_dim; k++)
                    sumorder[ind[k]] += (k+1);
                if(*verbose) {
                    if(count == itempQ) {
                        Rprintf("%3d percent done.\n", progress*10);
                        itempQ += ftrunc((double) total/10); progress++;
                        R_FlushConsole();
                    }
                    count++;
                }
            }
            if (n_extra > 1) {
                /* store probability and mean order */
                for (j = 0; j < n_dim; j++) {
                    itemp = 0;
                    for (k = 0; k < n_extra; k++)
                        if (maxdim[k] == (j+1))  // this should be j
                            itemp++;
                    prob[itempP++] = ((double) itemp / (double) n_extra);
                    order[itempO++] = ((double) sumorder[j] / (double) n_extra);
                }
            } else {
                /* store choice */
                for (j = 0; j < n_dim; j++) {
                    if (maxdim[0] == (j+1)) {
                        choice[itempC++] = j;
                        probTemp[j]++;
                    }
                    order[itempO++] = sumorder[j];
                }
            }
        }
        if (n_extra == 1)
            for (j = 0; j < n_dim; j++)
                prob[itempP++] = ((double) probTemp[j] / (double) n_draw);
    }
    
    /** write out the random seed **/
    PutRNGstate();
    
    /* freeing memory */
    FreeMatrix(X, n_samp*n_dim);
    free(vtemp);
    free(Xbeta);
    FreeMatrix(W, n_extra);
    FreeMatrix(beta, n_draw);
    FreeMatrix(mtemp, n_dim);
    Free3DMatrix(Sigma, n_draw, n_dim-1);
    free(maxdim);
    free(ind);
    free(sumorder);
    free(probTemp);
}

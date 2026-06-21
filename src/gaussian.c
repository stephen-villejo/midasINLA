#include <assert.h>
#include <strings.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cgeneric.h"

#define Calloc(n_, type_) (type_ *)calloc((n_), sizeof(type_))

double *inla_cgeneric_gaussian_midas(inla_cgeneric_cmd_tp cmd,
                                    double *theta,
                                    inla_cgeneric_data_tp *data)
{

    double *ret = NULL;

	assert(!strcasecmp(data->ints[0]->name, "n"));
	int N = data->ints[0]->ints[0];

	assert(!strcasecmp(data->ints[1]->name, "debug"));
	int debug = data->ints[1]->ints[0];

	assert(!strcasecmp(data->ints[2]->name, "K"));
	int K = data->ints[2]->ints[0];

	double prec_high = exp(15.0);
    double mu_val = (theta ? K * (1.0 / (1 + exp(-theta[0]))) : NAN);
    double sigma_val = (theta ? exp(theta[1]) : NAN);
	double beta1 = (theta ? theta[2] : NAN);
    
    double *w = NULL;
	if (theta) {
		w = Calloc(K+1, double);
		double sumw = 0.0;
		for (int k = 0; k <= K; k++) {
            w[k] = exp(-((k - mu_val)*(k - mu_val)) / (2 * sigma_val * sigma_val));
			sumw += w[k];
		}
		sumw = 1.0 / sumw;
		for (int k = 0; k <= K; k++) {
			w[k] *= sumw;
		}
	}

    switch(cmd) {
        case INLA_CGENERIC_GRAPH:
            ret = Calloc(2 + 2*N, double);
            ret[0] = N;
            ret[1] = N;
            for (int i = 0; i < N; i++) {
                ret[2 + i] = i;
                ret[2 + N + i] = i;
            }
            break;
        
        case INLA_CGENERIC_Q:
            ret = Calloc(2 + N, double);
            ret[0] = -1;
            ret[1] = N;
            for (int i = 0; i < N; i++) {
                ret[2 + i] = prec_high;
            }
            break;

        case INLA_CGENERIC_MU:
            ret = Calloc(N + 1, double);
            ret[0] = N;
            double *x = data->mats[0]->x;

            if (debug) {
                static int first =  1;
                if (first) {
                    printf("(nrow ncol) = (%d %d)\n", data->mats[0]->nrow,  data->mats[0]->ncol);
                }
                first = 0;
            }

            for (int i = 0; i < N; i++) {
                double agg = 0.0;
                for (int k = 0; k <= K; k++) {
                    agg += w[k] * x[i + k*N];  
                }
                ret[1 + i] = beta1 * agg;
            }
            break;

        case INLA_CGENERIC_LOG_PRIOR:
            ret = Calloc(1, double);
            ret[0] = -0.5 * theta[0]*theta[0] + -0.5 * theta[1]*theta[1] + -0.5 * theta[2]*theta[2];
            break;

        case INLA_CGENERIC_INITIAL:
            ret = Calloc(4, double);
            ret[0] = 3;
            ret[1] = 0.0;
            ret[2] = 0.0;
            ret[3] = 0.0;
            break;

        case INLA_CGENERIC_LOG_NORM_CONST:
            ret = Calloc(1, double);
            ret[0] = 0.0;
            break;

        case INLA_CGENERIC_QUIT:
	    case INLA_CGENERIC_VOID:
        default:
            break;
        }

        if (w) {
            free(w);
        }

    return ret;

}
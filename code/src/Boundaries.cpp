#include <iostream>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <cstring>

#include "Boundaries.h"

void ApplySpongeLayer(State& S1)
{
    unsigned int N = S1.N;

    // left sponge layer
    #pragma omp parallel for
    for(unsigned int i=0; i<nSP; i++)
    {
        float Cs = pow(alphaSP, pow(gammaSP, i));

        S1.Eta[i] = S1.Eta[i] / Cs;
        S1.Hn[i] = S1.Eta[i] - S1.Zn[i];

        S1.Pn[i] = S1.Pn[i] / Cs;
        S1.Un[i] = S1.Un[i] / Cs;
    }

    // right sponge layer
    if(S1.Grid == grid::staggered)
    {
        #pragma omp parallel for
        for(unsigned int i=N-nSP; i<N; i++)
        {
            float Cs = pow(alphaSP, pow(gammaSP, N-1-i));

            S1.Eta[i] = S1.Eta[i] / Cs;
            S1.Hn[i] = S1.Eta[i] - S1.Zn[i];

            S1.Pn[i+1] = S1.Pn[i+1] / Cs;
            S1.Un[i+1] = S1.Un[i+1] / Cs;
        }
    }
    else
    {
        #pragma omp parallel for
        for(unsigned int i=N-nSP; i<N; i++)
        {
            float Cs = pow(alphaSP, pow(gammaSP, N-1-i));

            S1.Eta[i] = S1.Eta[i] / Cs;
            S1.Hn[i] = S1.Eta[i] - S1.Zn[i];

            S1.Pn[i] = S1.Pn[i] / Cs;
            S1.Un[i] = S1.Un[i] / Cs;
        }
    }
    
}
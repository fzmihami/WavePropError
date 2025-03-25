#include <iostream>
#include <omp.h>
#include <math.h>
#include <cstring>

#include "Dispersion.h"

Dispersion::Dispersion(const State& S0)
{
    this->N = S0.N;
    this->PhiC = new float[N]; memset(this->PhiC, 0, N * sizeof(float));

    if(S0.Grid == grid::staggered)
    {
        this->PhiM = new float[N+1]; memset(this->PhiM, 0, (N+1) * sizeof(float));
        this->Uxx = new float[N+1]; memset(this->Uxx, 0, (N+1) * sizeof(float));

        this->D1 = new float[N+1]; memset(this->D1, 0, (N+1) * sizeof(float));
        this->D2 = new float[N+1]; memset(this->D2, 0, (N+1) * sizeof(float));
        this->D3 = new float[N+1]; memset(this->D3, 0, (N+1) * sizeof(float));
    }
    else
    {
        this->PhiM = new float[N]; memset(this->PhiM, 0, N * sizeof(float));
        this->Uxx = new float[N]; memset(this->Uxx, 0, N * sizeof(float));

        this->D1 = new float[N]; memset(this->D1, 0, N * sizeof(float));
        this->D2 = new float[N]; memset(this->D2, 0, N * sizeof(float));
        this->D3 = new float[N]; memset(this->D3, 0, N * sizeof(float));
    }

    // the diagonal terms of the diseprion matrix can be precomputed
    // since they only depend on the still water level, which assumed to be constant in time
    float dn = S0.depth;
    unsigned int i_end = S0.Grid == grid::staggered ? N+1 : N;

    for(unsigned int i=0; i<i_end; i++) { D2[i] = 1.0 - 2*alpha1*dn*dn * S0.dx_1*S0.dx_1;}
	for(unsigned int i=0; i<i_end-1; i++) { D1[i] = alpha1*dn*dn*S0.dx_1*S0.dx_1; }
	for(unsigned int i=1; i<i_end; i++) { D3[i] = alpha1*dn*dn*S0.dx_1*S0.dx_1; }
}

Dispersion::~Dispersion()
{
    delete[] PhiC;
    delete[] PhiM;
    delete[] Uxx;
}

void Dispersion::ComputeDispersion(State& S0)
{
    // compute second-order derivative of the velocity field
    unsigned int i_end = S0.Grid == grid::staggered ? N+1 : N;
    #pragma omp parallel for
    for (unsigned int i = 1; i < i_end-1; i++)
    {
        Uxx[i] = (S0.Un[i+1] - 2.0 * S0.Un[i] + S0.Un[i-1]) * S0.dx_1 * S0.dx_1;
    }

    // compute dispersive terms
    if(S0.Grid == grid::staggered)
    {
        #pragma omp parallel for
        for (unsigned int i = 1; i < N; i++)
        {
            float h = S0.depth;
            float Za = S0.depth*gamma1;

            PhiC[i] = (pow(Za, 2)/2.0f - pow(h, 2)/6.0f) * h * (Uxx[i+1] - Uxx[i]) * S0.dx_1
                        + (Za + h/2.0f) *h*h*  (Uxx[i+1] - Uxx[i]) * S0.dx_1;

            PhiM[i] = pow(Za, 2)/2.0f * Uxx[i] + Za * h * Uxx[i];
        }

        PhiC[N-1] = 0;
    }

    else {
        #pragma omp parallel for
        for (unsigned int i = 1; i < N-1; i++)
        {
            float h = S0.depth;
            float Za = S0.depth*gamma1;

            PhiC[i] = (pow(Za, 2)/2.0f - pow(h, 2)/6.0f) * h * (Uxx[i+1] - Uxx[i-1]) * (0.5f * S0.dx_1)
                        + (Za + h/2.0f) *h*h*  (Uxx[i+1] - Uxx[i-1]) * (0.5f * S0.dx_1);

            PhiM[i] = pow(Za, 2)/2.0f * Uxx[i] + Za * h * Uxx[i];

        }
    }

}

void Dispersion::SolveTridiagonalMatrix(State& S1)
{
    // solve the tridiagonal matrix
    unsigned int i_end = S1.Grid == grid::staggered ? N+1 : N;

    // compute Un
    if(S1.Grid == grid::staggered)
    {
        #pragma omp parallel for
        for (unsigned int i = 0; i < N+1; i++)
        {
            float H1;
            if(i ==0 ) H1 = S1.Hn[i];
            else if(i == N) H1 = S1.Hn[i-1];
            else H1 = 0.5*(S1.Hn[i] + S1.Hn[i-1]);

            if(H1 > hmin)
            {
                S1.Un[i] = S1.Pn[i] / H1;
            }
            else
            {
                S1.Un[i] = 0.0f;
            }
        }
    }
    else
    {
        #pragma omp parallel for
        for (unsigned int i = 0; i < N; i++)
        {
            if(S1.Hn[i] > hmin)
            {
                S1.Un[i] = S1.Pn[i] / S1.Hn[i];
            }
            else
            {
                S1.Un[i] = 0.0f;
            }
        }
    }

    // solve the tridiagonal matrix
    ThomasAlgorithm(D1, D2, D3, S1.Un, i_end);


}

void ThomasAlgorithm(const float* a, const float* b, const float* c, float* d, int n) {

    // create local arrays for the forward sweep
    float* c_prime = new float[n];

    // forward sweep
    c_prime[0] = c[0] / b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; ++i) {
        float denom = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = c[i] / denom;
        d[i] = (d[i] - a[i] * d[i - 1]) / denom;
    }

    // back substitution
    d[n - 1] = d[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        d[i] = d[i] - c_prime[i] * d[i + 1];
    }

    delete[] c_prime;
}
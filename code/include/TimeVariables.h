#pragma once

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include "ReadUserInput.h"
#include "GlobalVariables.h"

class State
{
    public:
        unsigned int N;
        float dx;
        float dx_1;
        float depth;
        grid Grid;

    public:
        float *Hn, *Eta, *Zn;
        float *Pn, *Un;

    public:
        State(const Input& In);
        ~State();

        void Update(State& other);
};


void PrintVariables(State& S0, State& S2, const Input &In, const float& t1, const float& dt1, unsigned int& Kprint, float& Tprint);
void printINFO(const Input& In);

void RecordTimeSeries(State& S0, State& S2, const Input &In, std::vector<float> *Hgauges, const float& t1, const float& dt1, unsigned int& Kgauges, float& Tgauges);
void PrintTimeSeries(std::vector<float> *Hgauges, const Input &In);

void linear_interpolation(float *H1, float *H2, unsigned int N, const float& t,const float& dt, const float& Tprint);
float linear_interpolation(float value1, float value2, const float& t,const float& dt, const float& Tprint);
std::string intToStringWithLeadingZeros(int number, int totalLength);


void printProgress(double percentage);
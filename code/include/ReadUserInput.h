#pragma once

#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <fstream>

#include "GlobalVariables.h"

class Input
{
    public:
        scheme Scheme;
        unsigned int OrderReconstruction;
        unsigned int OrderTime;
        float CourantNumber;
        float dx;
        float LengthDomain;
        float WaterDepth;
        float Time;
        float dtOutput;
        float dtGauges;
        std::vector<unsigned int> IndexGauges;
        std::vector<std::vector<float>> WaveInput;
        std::string ResFolder;

    public:
        Input(std::string nametest);
        void LOGOUT();
};

scheme convertScheme(std::string string1);

std::vector<unsigned int> convertIndexGauges(std::string string1, float dx);

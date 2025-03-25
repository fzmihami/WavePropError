#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <omp.h>
#include <cmath>

#include "ReadUserInput.h"
#include "GlobalVariables.h"

Input::Input(std::string nametest)
{
    std::ifstream file;
    std::vector<std::string> parameters;  // Array to store all parameters

    // read all lines from file
    file.open(nametest);
    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            line = line.substr(line.find(":") + 1);
            parameters.push_back(line);
        }
        file.close();
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
    }

    // Assign values to variables
    Scheme = convertScheme(parameters[0]);
    OrderReconstruction = std::stoi(parameters[1]);
    OrderTime = std::stoi(parameters[2]);
    CourantNumber = std::stof(parameters[3]);
    dx = std::stof(parameters[4]);
    LengthDomain = std::stof(parameters[5]);
    WaterDepth = std::stof(parameters[6]);
    Time = std::stof(parameters[7]);
    dtOutput = std::stof(parameters[8]);
    IndexGauges = convertIndexGauges(parameters[9], dx);
    dtGauges = std::stof(parameters[10]);

    // read wave input data as 2D vector of floats from string array
    for (int i = 12; i < parameters.size(); i++)
    {
        std::vector<float> row;
        std::string temp = parameters[i];
        std::string delimiter = " ";
        size_t pos = 0;
        std::string token;
        while ((pos = temp.find(delimiter)) != std::string::npos)
        {
            token = temp.substr(0, pos);
            row.push_back(std::stof(token));
            temp.erase(0, pos + delimiter.length());
        }
        row.push_back(std::stof(temp));
        WaveInput.push_back(row);
    }

    // define and create results folder
    std::string resultsFolder = nametest;
    resultsFolder.replace(0, 7, "results/");
    resultsFolder = resultsFolder.substr(0, resultsFolder.find("."));   //remove file extension

    std::string command = "rm -rf " + resultsFolder;
    system(command.c_str());

    command = "mkdir " + resultsFolder;
    system(command.c_str());

    command = "mkdir " + resultsFolder + "/FreeSurface";
    system(command.c_str());

    command = "mkdir " + resultsFolder + "/TimeSeries";
    system(command.c_str());

    ResFolder = resultsFolder;

}

void Input::LOGOUT()
{
    std::cout << std::endl;
    std::cout << "Input parameters" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout << "Numerical scheme :" << ListeScheme[Scheme] << std::endl;
    std::cout << "Order of reconstruction :" << OrderReconstruction << std::endl;
    std::cout << "Order of RK time integration :" << OrderTime << std::endl;
    std::cout << "Courant number :" << CourantNumber << std::endl;
    std::cout << "Grid spacing :" << dx << " m" << std::endl;
    std::cout << "Length of the domain :" << LengthDomain << " m" << std::endl;
    std::cout << "Water depth :" << WaterDepth << " m" << std::endl;
    std::cout << "Total computation time :" << Time << " s" << std::endl;
    std::cout << "Time step for the output :" << dtOutput << " s" << std::endl;
    std::cout << "Time step for the gauges :" << dtGauges << " s" << std::endl;

    std::cout << std::left << std::setw(25) << "frequency [Hz]" << std::setw(25) << "amplitude [m]" << std::setw(25) << "phase [rad]" << std::endl;
    for (int i = 0; i < WaveInput.size(); i++)
    {
        for (int j = 0; j < WaveInput[i].size(); j++)
        {
            std::cout << std::setw(25) << WaveInput[i][j];
        }
        std::cout << std::setw(25) << std::endl;

        if (i == 3)
        {
            std::cout << "..." << std::endl;
            break;
        }
    }
}

// function to convert string to enum scheme
scheme convertScheme(std::string string1)
{
    // remove leading and trailing whitespaces
    string1.erase(std::remove(string1.begin(), string1.end(), ' '), string1.end());

    int index=-1;
    for(int i=0; i<4; i++){
        if(string1 == ListeScheme[i]){
            index = i;
            break;
        }
    }

    if(index < 0){
        std::cout << "Scheme not found" << std::endl;
        exit(1);
    }

    return static_cast<scheme>(index);
}

// function to convert string to vector of unsigned int gauges index
std::vector<unsigned int> convertIndexGauges(std::string string1, float dx)
{
    // remove leading and trailing whitespaces
    string1.erase(std::remove(string1.begin(), string1.end(), ' '), string1.end());

    std::vector<float> gauges;
    std::string temp = string1.substr(1, string1.size() - 2);  // remove brackets
    if (temp.empty())
    {
        std::vector<unsigned int> gaugesIndex;
        return gaugesIndex;
    }

    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    while ((pos = temp.find(delimiter)) != std::string::npos)
    {
        token = temp.substr(0, pos);
        gauges.push_back(std::stoi(token));
        temp.erase(0, pos + delimiter.length());
    }
    gauges.push_back(std::stoi(temp));

    // convert float to unsigned int
    std::vector<unsigned int> gaugesIndex;
    for (int i = 0; i < gauges.size(); i++)
    {
        gaugesIndex.push_back(std::round(gauges[i] / dx));
    }

    return gaugesIndex;
}


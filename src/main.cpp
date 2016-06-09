#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include "Individual.h"
#include "diversity_functions.h"
#include "initialisation_functions.h"
#include "test_functions.h"
#include <math.h>
#include "RandomVectorGenerator.h"
#include <sstream>

#define PROBLEM_SIZE 2
#define POPULATION_SIZE 100
#define MIN_SPACE -10
#define MAX_SPACE 10

std::vector<double> addVectors(std::vector<double> v1, std::vector<double> v2) {
    std::vector<double> newV;

    if (v1.size() == 0) {
        v1.resize(v2.size());
    }

    for (unsigned long i = 0; i < v1.size(); i++) {
        newV.push_back(v1.at(i) + v2.at(i));
    }

    return newV;
}

std::vector<double> averageVector(std::vector<double> v, int n) {
    std::vector<double> newV;

    for (unsigned long i = 0; i < v.size(); i++) {
        newV.push_back(v.at(i) /= n);
    }

    return newV;
}

double initialiseAndAverageDiversity(std::function<std::vector<Individual>(int, int, double, double,
    std::function<double(Individual)>, RandomVectorGenerator)> init_function,
 std::function<double(std::vector<Individual>, int, float, float)> diversity_func) {

    double measure;

    for (int rI = 0; rI < 30; rI++) {
         RandomVectorGenerator rvca(rI);
        std::vector<Individual> population = init_function(POPULATION_SIZE, PROBLEM_SIZE, MIN_SPACE, MAX_SPACE, sphereFunction,rvca);

        double m = diversity_func(population, PROBLEM_SIZE, MIN_SPACE, MAX_SPACE);

        measure += m;
    }

    return (measure / 30);
}

void printAndCalculatePopulationDiversity(std::vector<Individual> population,
    std::function<double(std::vector<Individual>, int)> diversityFunc) {

    double diversity = diversityFunc(population, PROBLEM_SIZE);

    std::cout << diversity << "\n";
}

// void printCSVLine() {
//     //Print
//     ofstream csvFile;
//     csvFile.open(filename);
//     csvFile << "Generation" << ",";
//     //Print Scales
//     for (double x = MIN_SCALE;  x <= MAX_SCALE; x += SCALE_INC) {
//         csvFile << x << ",";
//     }
//     csvFile << "\n";

//     int gI = 0;
//     for (auto &vars : variances) {
//         csvFile << gI << ",";
//         for (auto &var : vars) {
//             csvFile << var << ",";
//         }
//         csvFile << "\n";
//         gI++;
//     }
//     csvFile.close();
// }

void printMethodLine(std::string title, double variance, double max, double mean, double radius, double averageDistance, double dgea, double entropy) {
    // cout << title << "          " << diversity << "\n";
     // printf("%20s %15f %15f %15f %15f \n", title.c_str(), log(variance), log(max), log(mean), log(radius));
     printf("%20s %15f %15f %15f %15f %15f %15f %15f \n", title.c_str(), variance, max, mean, radius, averageDistance,  dgea, entropy);
}

void writeVizFile(std::vector<Individual> population, std::string label) {
    std::ofstream csvFile;
    std::string filename = "./viz_output/data/"+label+".csv";
    csvFile.open (filename, std::ios::out);
    if(!csvFile.is_open())
    {
        std::cout << "Unable to create file " << filename << std::endl;
        //print error;
        return;
    }
     csvFile << "d1,d2" << "\n";
    for (size_t i = 0; i < population.size(); i++) {

        std::vector<double> atts = population.at(i).atts;
        for (int x = 0; x < PROBLEM_SIZE; x++) {
            if (x < PROBLEM_SIZE-1) {
                csvFile << atts.at(x) << ",";
            } else {
                csvFile << atts.at(x);
            }
        }
        csvFile.flush();
        csvFile << "\n";
    }
    csvFile.close();
}

void findMinMax(std::function<std::vector<Individual>(int, int, double, double,
    std::function<double(Individual)>, RandomVectorGenerator)> init_function,
 std::function<double(std::vector<Individual>, int, float, float)> diversity_func, std::string label) {

    double low = 0;
    double high = 0;
    std::vector<Individual> lowPopulation;
    std::vector<Individual> highPopulation;
    for (int i = 0; i < 1000; i++) {

        RandomVectorGenerator rvca(i);
        std::vector<Individual> population = init_function(POPULATION_SIZE, PROBLEM_SIZE, MIN_SPACE, MAX_SPACE, sphereFunction,rvca);
        double diversity = diversity_func(population, PROBLEM_SIZE, MIN_SPACE, MAX_SPACE);

        if (i == 0) {
            low = diversity;
            high = diversity;
            lowPopulation = population;
            highPopulation = population;
        } else {
            if (diversity < low) {
                low = diversity;
                lowPopulation = population;
            }
            if (diversity > high) {
                high = diversity;
                highPopulation = population;
            }
        }
    }


    //Print CSV
    std::ostringstream lowstr;
    std::ostringstream highstr;
    lowstr << low;
    highstr << high;

    writeVizFile(lowPopulation, label+"_low_"+lowstr.str());
    writeVizFile(highPopulation, label+"_high_"+highstr.str());

    std::cout << label << ":" << "Low: " << low << ", High: " << high << "\n";
}


void viz() {
    std::cout << "QRS Initialisation" << std::endl;
    findMinMax(sobolQRSInitialisation, varianceDiversity, "QRSVariance");
    findMinMax(sobolQRSInitialisation, maxDiversity, "QRSMax");
    findMinMax(sobolQRSInitialisation, meanDiversity, "QRSMean");
    findMinMax(sobolQRSInitialisation, swarmRadiusDiversity, "QRSRadius");
    findMinMax(sobolQRSInitialisation, averageDistanceAroundCentreDiversity, "QRSAverageDistance");
    findMinMax(sobolQRSInitialisation, dgeaDiversity, "QRSDGEA");
    findMinMax(sobolQRSInitialisation, entropyDiversity, "QRSEntropy");
    std::cout << "Random Initialisation" << std::endl;
    findMinMax(randomInitialisation, varianceDiversity, "RandomVariance");
    findMinMax(randomInitialisation, maxDiversity, "RandomMax");
    findMinMax(randomInitialisation, meanDiversity, "RandomMean");
    findMinMax(randomInitialisation, swarmRadiusDiversity, "RandomRadius");
    findMinMax(randomInitialisation, averageDistanceAroundCentreDiversity, "RandomAverageDistance");
    findMinMax(randomInitialisation, dgeaDiversity, "RandomDGEA");
    findMinMax(randomInitialisation, entropyDiversity, "RandomEntropy");
    std::cout << "Chaotic Initialisation" << std::endl;
    findMinMax(chaoticBasedInitialisation, varianceDiversity, "ChaoticVariance");
    findMinMax(chaoticBasedInitialisation, maxDiversity, "ChaoticMax");
    findMinMax(chaoticBasedInitialisation, meanDiversity, "ChaoticMean");
    findMinMax(chaoticBasedInitialisation, swarmRadiusDiversity, "ChaoticRadius");
    findMinMax(chaoticBasedInitialisation, averageDistanceAroundCentreDiversity, "ChaoticAverageDistance");
    findMinMax(chaoticBasedInitialisation, dgeaDiversity, "ChaoticDGEA");
    findMinMax(chaoticBasedInitialisation, entropyDiversity, "ChaoticEntrop");
    std::cout << "Opposition Initialisation" << std::endl;
    findMinMax(oppositionBasedInitialisation, varianceDiversity, "OppositionVariance");
    findMinMax(oppositionBasedInitialisation, maxDiversity, "OppositionMax");
    findMinMax(oppositionBasedInitialisation, meanDiversity, "OppositionMean");
    findMinMax(oppositionBasedInitialisation, swarmRadiusDiversity, "OppositionRadius");
    findMinMax(oppositionBasedInitialisation, averageDistanceAroundCentreDiversity, "OppositionAverageDistance");
    findMinMax(oppositionBasedInitialisation, dgeaDiversity, "OppositionDGEA");
    findMinMax(oppositionBasedInitialisation, entropyDiversity, "OppositionEntropy");
    std::cout << "Chaotic-Opposition Initialisation" << std::endl;
    findMinMax(chaoticOppositionBasedInitialisation, varianceDiversity, "ChaoticOppositionVariance");
    findMinMax(chaoticOppositionBasedInitialisation, maxDiversity, "ChaoticOppositionMax");
    findMinMax(chaoticOppositionBasedInitialisation, meanDiversity, "ChaoticOppositionMean");
    findMinMax(chaoticOppositionBasedInitialisation, swarmRadiusDiversity, "ChaoticOppositionRadius");
    findMinMax(chaoticOppositionBasedInitialisation, averageDistanceAroundCentreDiversity, "ChaoticOppositionAverageDistance");
    findMinMax(chaoticOppositionBasedInitialisation, dgeaDiversity, "ChaoticOppositionDGEA");
    findMinMax(chaoticOppositionBasedInitialisation, entropyDiversity, "ChaoticOppositionEntropy");
    std::cout << "Orthogonal Design Initialisation" << std::endl;
    findMinMax(orthogonalArray, varianceDiversity, "ODVariance");
    findMinMax(orthogonalArray, maxDiversity, "ODMax");
    findMinMax(orthogonalArray, meanDiversity, "ODMean");
    findMinMax(orthogonalArray, swarmRadiusDiversity, "ODRadius");
    findMinMax(orthogonalArray, averageDistanceAroundCentreDiversity, "ODAverageDistance");
    findMinMax(orthogonalArray, dgeaDiversity, "ODDGEA");
    findMinMax(orthogonalArray, entropyDiversity, "ODEntropy");
}

int main(int argc, char ** argv)
{

    if (argc > 1 && strcmp(argv[1], "viz")==0) {
        viz();
        return 0;
    }

    printf("%20s %15s %15s %15s %15s %15s %15s %15s \n", "Method", "Variance", "Max-Diversity", "Mean-Diversity", "Swarm Radius", "ADAC", "DGEA", "Entropy");
    printf("-------------------------------------------------------------------------------------------------------------------------------------------\n");



    printMethodLine("QRS",
        initialiseAndAverageDiversity(sobolQRSInitialisation, varianceDiversity),
        initialiseAndAverageDiversity(sobolQRSInitialisation, maxDiversity),
        initialiseAndAverageDiversity(sobolQRSInitialisation, meanDiversity),
        initialiseAndAverageDiversity(sobolQRSInitialisation, swarmRadiusDiversity),
        initialiseAndAverageDiversity(sobolQRSInitialisation, averageDistanceAroundCentreDiversity),
        initialiseAndAverageDiversity(sobolQRSInitialisation, dgeaDiversity),
        initialiseAndAverageDiversity(sobolQRSInitialisation, entropyDiversity));

    printMethodLine("Uniform-Random",
        initialiseAndAverageDiversity(randomInitialisation, varianceDiversity),
        initialiseAndAverageDiversity(randomInitialisation, maxDiversity),
        initialiseAndAverageDiversity(randomInitialisation, meanDiversity),
        initialiseAndAverageDiversity(randomInitialisation, swarmRadiusDiversity),
        initialiseAndAverageDiversity(randomInitialisation, averageDistanceAroundCentreDiversity),
        initialiseAndAverageDiversity(randomInitialisation, dgeaDiversity),
        initialiseAndAverageDiversity(randomInitialisation, entropyDiversity));

    printMethodLine("Chaotic",
        initialiseAndAverageDiversity(chaoticBasedInitialisation, varianceDiversity),
        initialiseAndAverageDiversity(chaoticBasedInitialisation, maxDiversity),
        initialiseAndAverageDiversity(chaoticBasedInitialisation, meanDiversity),
        initialiseAndAverageDiversity(chaoticBasedInitialisation, swarmRadiusDiversity),
        initialiseAndAverageDiversity(chaoticBasedInitialisation, averageDistanceAroundCentreDiversity),
        initialiseAndAverageDiversity(chaoticBasedInitialisation, dgeaDiversity),
        initialiseAndAverageDiversity(chaoticBasedInitialisation, entropyDiversity));

    printMethodLine("Opposition",
        initialiseAndAverageDiversity(oppositionBasedInitialisation, varianceDiversity),
        initialiseAndAverageDiversity(oppositionBasedInitialisation, maxDiversity),
        initialiseAndAverageDiversity(oppositionBasedInitialisation, meanDiversity),
        initialiseAndAverageDiversity(oppositionBasedInitialisation, swarmRadiusDiversity),
        initialiseAndAverageDiversity(oppositionBasedInitialisation, averageDistanceAroundCentreDiversity),
        initialiseAndAverageDiversity(oppositionBasedInitialisation, dgeaDiversity),
        initialiseAndAverageDiversity(oppositionBasedInitialisation, entropyDiversity));

    printMethodLine("Chaotic-Opposition",
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, varianceDiversity),
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, maxDiversity),
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, meanDiversity),
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, swarmRadiusDiversity),
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, averageDistanceAroundCentreDiversity),
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, dgeaDiversity),
        initialiseAndAverageDiversity(chaoticOppositionBasedInitialisation, entropyDiversity));

    // printMethodLine("Uniform Design",
    //     initialiseAndAverageDiversity(uniformDesignInitialisation, varianceDiversity, rvc),
    //     initialiseAndAverageDiversity(uniformDesignInitialisation, maxDiversity, rvc),
    //     initialiseAndAverageDiversity(uniformDesignInitialisation, meanDiversity, rvc),
    //     initialiseAndAverageDiversity(uniformDesignInitialisation, swarmRadiusDiversity, rvc),
    //     initialiseAndAverageDiversity(uniformDesignInitialisation, dgeaDiversity, rvc),
    //     initialiseAndAverageDiversity(uniformDesignInitialisation, entropyDiversity, rvc));

    printMethodLine("Orthogonal Design",
        initialiseAndAverageDiversity(orthogonalArray, varianceDiversity),
        initialiseAndAverageDiversity(orthogonalArray, maxDiversity),
        initialiseAndAverageDiversity(orthogonalArray, meanDiversity),
        initialiseAndAverageDiversity(orthogonalArray, swarmRadiusDiversity),
        initialiseAndAverageDiversity(orthogonalArray, averageDistanceAroundCentreDiversity),
        initialiseAndAverageDiversity(orthogonalArray, dgeaDiversity),
        initialiseAndAverageDiversity(orthogonalArray, entropyDiversity));


    return 0;
}

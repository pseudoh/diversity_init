#include <math.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "Individual.h"
#include "diversity_functions.h"

using namespace std;

double calculateVectorMagnitude(std::vector<double> v) {
    double length = 0;

    for (auto &x : v) {
        length += pow(x, 2);
    }

    return sqrt(length);
}

double entropyDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space) {
    // double sum = 0.0;
    // for (int i = 0; i < problem_size - population.size(); i++) {

    //     sum += log((problem_size/population.size()) * (population.at(i+population.size())))
    // }
    // return sum;
}


double swarmRadiusDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space){

    std::vector<double> centroid(problem_size,0);

    for(auto& individual : population){
        for(int d = 0;d<problem_size;++d) {
            centroid.at(d) += individual.atts.at(d);
        }
    }

    for (int i = 0; i < problem_size; ++i) {
        centroid.at(i) /= population.size();
    }

    double sum = 0.0;

    for (auto &individual : population) {
        for (int d = 0; d < problem_size; d++) {
            sum += pow(individual.atts.at(d) - centroid.at(d), 2);
        }
    }

    return sqrt(sum);
};


double meanDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space){
    double result = 0;

    std::vector<double> centroid(problem_size,0);

    for(auto& individual : population){
        for(int d = 0;d<problem_size;++d) {
            centroid.at(d) += individual.atts.at(d);
        }
    }

    for (int i = 0; i < problem_size; ++i) {
        centroid.at(i) /= population.size();
    }

    double jSum = 0.0;

    for (auto &individual : population) {
        for (int d = 0; d < problem_size; d++) {
            jSum += pow(individual.atts.at(d) - centroid.at(d), 2);
        }
    }
    return jSum / (population.size() * problem_size);
};

double maxDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space){



    double maxd[problem_size];

    for (auto &individual : population) {
       for (size_t i = 0; i < individual.atts.size(); i++) {
            double v = abs(individual.atts.at(i));
            if (v > maxd[i]) {
                maxd[i] = v;
            }


       }
    }

    //Centroid
    std::vector<double> xNorj(problem_size);
    for (auto &individual : population) {
       for (size_t i = 0; i < individual.atts.size(); i++) {
            xNorj.at(i) += (individual.atts.at(i) / maxd[i]) / population.size();
       }
    }

    double jSum = 0.0;

    for (auto &individual : population) {
        for (int d = 0; d < problem_size; d++) {
            jSum += abs((individual.atts.at(d) / maxd[d])  - xNorj.at(d));
        }
    }

    jSum /= population.size() * problem_size;

    return jSum;

};

double varianceDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space){

    std::vector< std::vector<double> > strategies;

    int individualN = 0;

    for (auto &individual : population) {
        std::vector<double> strategy;
        for (int d = 0; d < problem_size; d++) {
            strategy.push_back(individual.atts.at(d));
        }
        strategies.push_back(strategy);
    }



    double finalV = 0.0;
    for (int strategyDI = 0; strategyDI < problem_size; strategyDI++) {
        double sum = 0.0;

        for (auto &strategyDValue : strategies[strategyDI]) {
            sum += strategyDValue;
        }

        sum /= population.size();

        for (size_t strategyDValueI = 0; strategyDValueI < strategies[strategyDI].size(); strategyDValueI++) {
            finalV += pow(strategies[strategyDI].at(strategyDValueI) - sum, 2);

        }

    }

    finalV /= population.size();


    return finalV;
        };


double dgeaDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space){

   std::vector<double> centroid(problem_size,0);

    for(auto& individual : population){
        for(int d = 0;d<problem_size;++d) {
            centroid.at(d) += individual.atts.at(d);
        }
    }

    for (int i = 0; i < problem_size; ++i) {
        centroid.at(i) /= population.size();
    }

     double iSum = 0.0;

    for (size_t i = 0; i < population.size(); i++) {
        Individual individual = population.at(i);
        double jSum = 0;
        for (size_t j = 0; j < individual.atts.size(); j++) {
            jSum += pow(individual.atts.at(j) - centroid.at(j), 2);
        }
        iSum += sqrt(jSum);
    }

    double l = sqrt(pow(min_space, 2) + pow(max_space, 2));

    return iSum / (l  * population.size());
};

double averageDistanceAroundCentreDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space){

   std::vector<double> centroid(problem_size,0);

    for(auto& individual : population){
        for(int d = 0;d<problem_size;++d) {
            centroid.at(d) += individual.atts.at(d);
        }
    }

    for (int i = 0; i < problem_size; ++i) {
        centroid.at(i) /= population.size();
    }

     double iSum = 0.0;

    for (size_t i = 0; i < population.size(); i++) {
        Individual individual = population.at(i);
        double jSum = 0;
        for (size_t j = 0; j < individual.atts.size(); j++) {
            jSum += pow(individual.atts.at(j) - centroid.at(j), 2);
        }
        iSum += sqrt(jSum);
    }

    return iSum / population.size();
};

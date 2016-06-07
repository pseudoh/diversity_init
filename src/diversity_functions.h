#ifndef DIVERSITY_FUNCTIONS_H
#define DIVERSITY_FUNCTIONS_H

#include <vector>
#include "Individual.h"

double swarmRadiusDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

double meanDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

double maxDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

double varianceDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

double dgeaDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

double averageDistanceAroundCentreDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

double entropyDiversity(std::vector<Individual> population, int problem_size, float min_space, float max_space);

#endif

#ifndef INTIALISATION_FUNCTIONS_H
#define INTIALISATION_FUNCTIONS_H

#include <vector>
#include "Individual.h"
#include "RandomVectorGenerator.h"

std::vector<Individual> orthogonalArray(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

std::vector<Individual> randomInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

std::vector<Individual> sobolQRSInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

std::vector<Individual> uniformDesignInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

std::vector<Individual> chaoticOppositionBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

std::vector<Individual> chaoticBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

std::vector<Individual> oppositionBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc);

#endif

#include <math.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include "Individual.h"
#include "test_functions.h"

using namespace std;

double sphereFunction(Individual individual){
            double fitness = 0.0;

            for (auto &x : individual.atts) {
                fitness += pow(x, 2);
            }

            return fitness;
        };

double rosenbrockFunction(Individual individual){
            double fitness = 0.0;

            for (unsigned long i = 0; i < individual.atts.size() - 1; i++) {
                double xi = individual.atts.at(i);
                fitness += pow(100 * (individual.atts.at(i + 1) - pow(xi, 2)), 2) + pow(1 - xi, 2);
            }

            return fitness;
        };

double rastriginFunction(Individual individual){
            double fitness = 0.0;

            int a = 10;

            for (unsigned long i = 0; i < individual.atts.size(); i++) {
                double xi = individual.atts.at(i);
                fitness += pow(xi, 2) - a*cos(2*M_PI*xi);
            }

            return (individual.atts.size() * a) + fitness;
        };

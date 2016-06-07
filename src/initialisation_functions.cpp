
#include "initialisation_functions.h" // for plus
#include <vector>
#include "Individual.h"
#include "RandomVectorGenerator.h"
#include <math.h>
#include <algorithm> // for transform
#include <functional> // for plus
#include "sobol.h" // for plus
#include <iostream>

using namespace std;


bool isPrime(int num) {
    if (num <= 3) {
        return num > 1;
    } else if (num % 2 == 0 || num % 3 == 0) {
        return false;
    } else {
        for (int i = 5; i * i <= num; i += 6) {
            if (num % i == 0 || num % (i + 2) == 0) {
                return false;
            }
        }
        return true;
    }
}

std::vector<Individual> orthogonalArray(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){


    //Params
    int q = 3;
    int s = 12; //M = S * pow(q, j1)

    //Generate J1
    int j1 = 1;
    while ((pow(q, j1) - 1) / (q - 1) < problem_size) {
        j1++;
    }

    int n1 = problem_size;

    if ((pow(q, j1) - 1) / (q - 1) > problem_size) {
        n1 = (pow(q, j1) - 1) / (q - 1);
    }

    //Construct basic columns
    std::vector<std::vector<double>> aij(pow(q, j1));
    for (int k = 1; k <= j1; k++) {
        double j = ((pow(q, k - 1) - 1) / (q - 1)) + 1;
        for (int i = 1; i <= pow(q, j1); i++) {
            int x = fmod((((i - 1) / (pow(q, j1 - k)))), q);
            aij.at(i - 1).push_back(x);

        }

    }

    std::vector<std::vector<double>> aji(n1);

    for (size_t i = 0; i < aij.at(0).size();i++) {
        for (size_t x = 0; x < aij.size(); x++) {
            // cout << x << "\n";
            aji.at(i).push_back(aij.at(x).at(i));
        }
    }


    // Construct non basic columns
    for (int k = 2; k <= j1; k++) {
        double j = ((pow(q, k - 1) - 1) / (q - 1)) + 1;
        for (int s = 1; s <= j - 1; s++) {
            for (int t = 1; t <= q - 1; t++) {
                // as * t
                std::vector<double> as = aji.at(s - 1);
                std::vector<double> aj = aji.at(j - 1);
                std::vector<double> x;

                for (size_t i = 0; i < as.size(); i++) {
                    // x.push_back(fmod((as.at(i) * t) + aj.at(i), q));
                    aji.at((j + (s - 1) * (q - 1)+t) - 1).push_back(fmod((as.at(i) * t) + aj.at(i), q));

                }
            }
        }
    }

    // Increment aij
    for (int i = 0; i < n1; i++) {
        for (int x = 0; x < pow(q, j1); x++) {
            aji.at(i).at(x) += 1;
        }
    }

    //Delete the last N1 - N
    for (int i = 0; i < n1 - problem_size; i++) {
        aji.pop_back();
    }

    // Print out
    // for (int i = 0; i < aji.size(); i++) {
    //     for (int x = 0; x < pow(q, j1); x++) {
    //         cout <<  aji.at(i).at(x) << ",";
    //     }
    //     cout << "\n";
    // }


    //Initialisation
    std::vector<double> l;
    std::vector<double> u;
    std::vector<double> s1;

    for (int i = 0; i < problem_size; i++) {
        l.push_back(min_space);
        u.push_back(max_space);
        if (i == 0) {
            s1.push_back(1);
        } else {
            s1.push_back(0);
        }

    }
    // std::vector<double> l = {0.5, 3.5, 4.5};
    // std::vector<double> u = {10.5, 6.5, 7.5};
    // std::vector<double> s1 = {1,0,0};

    //Get max us-ls
    double usls = 0;
    for (int i = 0; i < problem_size; i++) {
        if (u.at(i) - l.at(i) > usls) {
            usls = u.at(i) - l.at(i) ;
        }
    }

    // cout << usls << "\n";

    //Divide subspaces
    std::vector<std::vector<std::vector<double>>> spaces;
    for (int i = 1; i <= s; i++) {
        std::vector<double> liv;
        std::vector<double> luv;
        std::vector<std::vector<double>> space;
        for (size_t x = 0; x < l.size(); x++) {
            double li = l.at(x) + ((i - 1) * (usls/s) * s1.at(x));
            double ui = u.at(x) - ((s - 1) * (usls/s) * s1.at(x));

            // cout << li << "-" << ui << ",";
            liv.push_back(li);
            luv.push_back(ui);
        }

        space.push_back(liv);
        space.push_back(luv);

        spaces.push_back(space);
    }

    // cout << "\n" << "Quantization" << "\n";
    //Quantize
    std::vector<std::vector<std::vector<double>>> quants(spaces.size());
    for (size_t si = 0; si < spaces.size(); si++) {

        std::vector<std::vector<double>> space = spaces.at(0);
        std::vector<double> ls = space.at(0);
        std::vector<double> us = space.at(1);

        // cout << "Space " << si + 1 << "\n";

        for (int i = 0; i < problem_size; i++) {

            std::vector<double> quant;
            double li = ls.at(i);
            double ui = us.at(i);

            quant.push_back(li); //li
            // cout << li << ",";
            for (int j = problem_size - 1; j <= problem_size - 1; j++) {
                double mi = li + (j - 1) * ((ui - li) / (q - 1));
                // cout << mi  << ",";
                quant.push_back(mi);
            }
            quant.push_back(ui); //li
            // cout << ui << "\n";
            quants.at(si).push_back(quant);
        }
    }


    // //Apply L3

    std::vector<Individual> population(pow(q, j1) * s);


    for (int si = 0; si < s; si++) {
        for (size_t i = 0; i < aji.size(); i++) {
            for (size_t x = 0; x < aji.at(i).size(); x++) {
                double lv = aji.at(i).at(x);
                std::vector<double> a = quants.at(si).at(i);
                // cout << lv - 1 << ",";
                // cout << a.at(lv - 1) << ",";

                population.at((aji.at(i).size() * si) + x).atts.push_back(a.at(lv - 1));
            }
            // cout << "\n";
        }
        // cout << "\n\n";
    }

    //Choose random = population_size
    std::random_shuffle ( population.begin(), population.end() );


    population.resize(population_size);

    return population;
};

std::vector<Individual> sobolQRSInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){

    int DIM_MAX = problem_size;

    std::vector<Individual> population;


  int dim_num;
  int i;
  int j;
  float r[DIM_MAX];
 int seed;
  int seed_in;
  int seed_out;

    seed = 0;

    for ( i = 0; i <= population_size - 1; i++ )
    {
      seed_in = seed;
      i4_sobol ( problem_size, &seed, r );
      seed_out = seed;

       Individual individual;
       for ( j = 0; j < problem_size; j++ )
        {
           individual.atts.push_back(min_space + (max_space - min_space) * r[j]);
        }

        population.push_back(individual);

    }

    return population;
};

std::vector<Individual> randomInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){

    std::vector<Individual> population;
    for (int i = 0; i < population_size; i++) {
        Individual individual;
        individual.atts = rvc.getUniformRandomVector(min_space, max_space, problem_size);
        population.push_back(individual);
    }
    return population;
};

std::vector<Individual> uniformDesignInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){



    double randomPrimes[population_size];

    //Find all the prime numbers which
    //are less than M , where M is the size of population.

    for (int i = 0; i < population_size; i++) {
        // double prime = std::round(rvc.getNormalRandomWithinRange(1, population_size - 1));
        // if (!isPrime(prime) || prime < population_size) {
        //     i--;
        // } else {
        //     randomPrimes[i] = prime;
        // }
        //Get every prime up to popualation size
    }

    std::vector<Individual> population;

    for (int i = 0; i < population_size; i++) {
        Individual individual;
        for (int j = 1; j <= problem_size; j++) {
            double x = fmod(i * randomPrimes[j], population_size);
            individual.atts.push_back(x);
        }
        population.push_back(individual);
    }

    return population;

};

std::vector<Individual> chaoticOppositionBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){

    std::vector<Individual> population;


    int max_k = 300;

    for (int i = 0; i < population_size; i++) {

        std::vector<std::vector<double>> ch;
        std::vector<double> chi;
        ch.push_back(chi);

        Individual individual;

        for (int j = 0; j < problem_size; j++) {

            double chi = rvc.getUniformRandom();

            for (int k = 1; k <= max_k; k++) {
                double s = sin(M_PI * chi);
                chi = s;
            }

            double att = min_space + (chi * (max_space - min_space));

            individual.atts.push_back(att);
        }


        population.push_back(individual);
    }


    for (int i = 0; i < population_size; i++) {
        Individual individual;

        for (int j = 0; j < problem_size; j++) {
            double ox = min_space + max_space - population.at(i).atts.at(j);
            individual.atts.push_back(ox);
        }


        population.push_back(individual);
    }

    std::random_shuffle ( population.begin(), population.end() );


    population.resize(population_size);



    return population;
};


std::vector<Individual> chaoticBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){

    std::vector<Individual> population;


    int max_k = 300;

    for (int i = 0; i < population_size; i++) {

        std::vector<std::vector<double>> ch;
        std::vector<double> chi;
        ch.push_back(chi);

        Individual individual;

        for (int j = 0; j < problem_size; j++) {

            double chi = rvc.getUniformRandom();

            for (int k = 1; k <= max_k; k++) {
                double s = sin(M_PI * chi);
                chi = s;
            }

            double att = min_space + (chi * (max_space - min_space));
            individual.atts.push_back(att);
        }

        population.push_back(individual);
    }



    return population;
};

std::vector<Individual> oppositionBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){


    std::vector<Individual> population = randomInitialisation(population_size, problem_size,
        min_space, max_space, fitness_func, rvc);

    for (int i = 0; i < population_size; i++) {
        Individual individual;

        for (int j = 0; j < problem_size; j++) {
            double ox = min_space + max_space - population.at(i).atts.at(j);
            individual.atts.push_back(ox);
        }

        population.push_back(individual);
    }

    std::random_shuffle ( population.begin(), population.end() );


    population.resize(population_size);


    return population;
};

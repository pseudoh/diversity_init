
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

std::vector<std::vector<double>> createBasicColumns(int j1, int q) {
    // cout << "Basic Columns" << "\n";
    std::vector<std::vector<double>> aij;
    for (int k = 1; k <= j1; k++) {
        double j = ((pow(q, k - 1) - 1) / (q - 1)) + 1;
        std::vector<double> jRow;
        for (int i = 1; i <= pow(q, j1); i++) {
            int x = fmod((((i - 1) / (pow(q, j1 - k)))), q);
            jRow.push_back(x);
            // cout << x << ",";
        }
        aij.push_back(jRow);
        // cout << "\n";
    }
    return aij;
}

std::vector<std::vector<double>> createNonBasicColumns(int j1, int q,std::vector<std::vector<double>> basicColumns) {
    // cout << "Non Basic Columns" << "\n";
    std::vector<std::vector<double>> aij;
    for (int k = 2; k <= j1; k++) {
        double j = ((pow(q, k - 1) - 1) / (q - 1)) + 1;

        for (int s = 1; s <= j - 1; s++) {
            for (int t = 1; t <= q - 1; t++) {
                // as * t
                std::vector<double> as = basicColumns.at(s - 1);
                std::vector<double> aj = basicColumns.at(j - 1);
                std::vector<double> x;

                for (size_t i = 0; i < as.size(); i++) {
                    // x.push_back(fmod((as.at(i) * t) + aj.at(i), q));
                    double nonbasicA = fmod((as.at(i) * t) + aj.at(i), q);
                    // aji.at((j + (s - 1) * (q - 1)+t) - 1).push_back(nonbasicA);
                    x.push_back(nonbasicA);
                    // cout << nonbasicA << ",";
                }
                // cout << "\n";
                aij.push_back(x);
            }

        }

    }
    return aij;
}


std::vector<std::vector<std::vector<double>>> generateSubspaces(
    int s,
    int problem_size,
    std::vector<double> minSpace,
    std::vector<double> maxSpace) {

    //Get max us-ls
    double usls = 0;
    for (int i = 0; i < problem_size; i++) {
        if (maxSpace.at(i) - minSpace.at(i) > usls) {
            usls = maxSpace.at(i) - minSpace.at(i) ;
        }
    }


    std::vector<double> s1(problem_size, 0); //1s
    s1.at(0) = 1;

    std::vector<std::vector<std::vector<double>>> spaces;
    for (int i = 1; i <= s; i++) {
        std::vector<double> liv;
        std::vector<double> luv;
        std::vector<std::vector<double>> space;
        for (size_t x = 0; x < problem_size; x++) {
            double li = minSpace.at(x) + ((i - 1) * (usls/s) * s1.at(x));
            double ui = maxSpace.at(x) - ((s - i) * (usls/s) * s1.at(x));

            liv.push_back(li);
            luv.push_back(ui);
        }

        space.push_back(liv);
        space.push_back(luv);

        spaces.push_back(space);
    }

    return spaces;
}


std::vector<std::vector<std::vector<double>>> quantize(
    int problem_size,
    int q,
    std::vector<std::vector<std::vector<double>>> spaces) {

    std::vector<std::vector<std::vector<double>>> quants(spaces.size());
    for (size_t si = 0; si < spaces.size(); si++) {

        std::vector<std::vector<double>> space = spaces.at(si);
        std::vector<double> ls = space.at(0);
        std::vector<double> us = space.at(1);

        cout << "Quants " << si + 1 << "\n";

        for (int i = 0; i < problem_size; i++) {

            std::vector<double> quant;
            double li = ls.at(i);
            double ui = us.at(i);

            quant.push_back(li); //li
            cout << li << ",";
            for (int j = problem_size - 1; j <= problem_size - 1; j++) {
                double mi = li + (j - 1) * ((ui - li) / (q - 1));
                cout << mi  << ",";
                quant.push_back(mi);
            }
            quant.push_back(ui); //li
            cout << ui << "\n";
            quants.at(si).push_back(quant);
        }
    }

    return quants;
}

std::vector<Individual> orthogonalArray(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){


    std::vector<Individual> population;


    double levels_per_dimension = ceil(pow(population_size,1.0/problem_size));

    double n = pow(levels_per_dimension, problem_size);

    // std::cout << levels_per_dimension << "," << n << "\n";

    for (int i = 1; i < population_size - 1; i++) {

        Individual individual;
        for (int j = 0; j < problem_size; j++) {

            int level = 1 + (rvc.getUniformRandom() * (n - 1));

            double f = 1+(fmod(floor(level  / pow(levels_per_dimension, j)),levels_per_dimension));
             // double randomPoint = min_space + (((f - 2) * (1/(levels_per_dimension-1))) * (max_space - min_space));
            // cout <<  f * (1/(levels_per_dimension-1)) << ",";
            individual.atts.push_back(f);
        }
        // std::cout <<  "\n";

        population.push_back(individual);
    }
    return population;

}


std::vector<Individual> orthogonalArraya(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){


    //Params
    int q = 3;
    int s = 5; //M = S * pow(q, j1)

    //***
    //*** Algorithm 2
    //***

    //Step 1: Select Smallest j1
    int j1 = 1;
    while ((pow(q, j1) - 1) / (q - 1) < problem_size) {
        j1++;
    }

    //Step 2: Calculate N1
    int n1 = problem_size;
    if ((pow(q, j1) - 1) / (q - 1) > problem_size) {
        n1 = (pow(q, j1) - 1) / (q - 1);
    }

    //Step 3: Execute Algorithm 1

    //***
    //*** Algorithm 1
    //***

    //Construct Basic Columns
    std::vector<std::vector<double>> basicColumns = createBasicColumns(j1, q);
    std::vector<std::vector<double>> nonBasicColumns = createNonBasicColumns(j1, q, basicColumns);


    //Concat Columns
    std::vector<std::vector<double>> columns;
    std::copy(basicColumns.begin(), basicColumns.end(), std::back_inserter(columns));
    std::copy(nonBasicColumns.begin(), nonBasicColumns.end(), std::back_inserter(columns));

    //Increment Columns
    for (int i = 0; i < columns.size(); i++) {
        for (int x = 0; x < columns.at(i).size(); x++) {
            columns.at(i).at(x)++;
        }
    }

    //Step 4: Delete the last N1-N
    for (int i = 0; i < n1 - problem_size; i++) {
        columns.pop_back();
    }

    cout << "Columns" << "\n";
    for (int i = 0; i < columns.size(); i++) {
        for (int x = 0; x < columns.at(i).size(); x++) {
            cout << columns.at(i).at(x) << ",";
        }
        cout << "\n";
    }

    //Columns = L32(33)


    //***
    //*** Algorithm 3
    //***

    //Step 1: Divide to Subspaces
    // std::vector<double> minSpaceVector = {0.5, 3.5, 4.5};
    // std::vector<double> maxSpaceVector = {10.5, 6.5, 7.5};

    std::vector<double> minSpaceVector(problem_size, min_space);
    std::vector<double> maxSpaceVector(problem_size, max_space);

    std::vector<std::vector<std::vector<double>>> spaces = generateSubspaces(s, problem_size, minSpaceVector, maxSpaceVector);

    cout << "\nSpaces" << "\n";
    for (int i = 0; i < spaces.size(); i++) {
        cout << "Space " << i << " \n";
        for (int x = 0; x < spaces.at(i).size(); x++) {
            for (int y = 0; y < spaces.at(i).at(x).size(); y++) {
                cout << spaces.at(i).at(x).at(y) << ",";
            }
            cout << "\n";
        }
        cout << "\n";
    }

    //Step 2: Quantize
    std::vector<std::vector<std::vector<double>>> quants = quantize(problem_size, q, spaces);

    //Step 3: Select Chromosomes - Apply L3w(33)

    std::vector<Individual> population(s * pow(q, j1));

    for (int si = 0; si < s; si++) {
        for (size_t i = 0; i < columns.size(); i++) {
            // Individual individual;
            for (size_t j = 0; j < columns.at(i).size(); j++) {
                double lv = columns.at(i).at(j);
                std::vector<double> a = quants.at(si).at(i);
                population.at((columns.at(i).size() * si) + j).atts.push_back(a.at(lv - 1));
            }
        }
    }

    // std::random_shuffle ( population.begin(), population.end());
    if (s * pow(q, j1) > population_size) {
        population.resize(population_size);
    }

    cout << "\nPopulation: " << population.size() << "\n";
    for (int i = 0; i < population.size(); i++) {
        for (int x = 0; x < population.at(i).atts.size(); x++) {
            cout << population.at(i).atts.at(x) << ",";
        }
        cout << "\n";
    }


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


std::vector<Individual> oldChaoticBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
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

std::vector<Individual> chaoticBasedInitialisation(int population_size, int problem_size, double min_space, double max_space,
    std::function<double(Individual)> fitness_func, RandomVectorGenerator rvc){

    std::vector<Individual> population;


    int max_k = 300;

    Individual individual;

    for (int d  = 0; d < problem_size; d++) {
        double att = min_space + (rvc.getUniformRandom() * (max_space - min_space));
        individual.atts.push_back(att);
    }

    population.push_back(individual);

    for (int i = 1; i < population_size - 1; i++) {



        Individual individual;
        Individual previousIndividual = population.at(i - 1);

        for (int j = 0; j < problem_size; j++) {


            double s = 1.0 * sin(M_PI * previousIndividual.atts.at(j));

            double att = min_space + (s * (max_space - min_space));
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

#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <gtest/gtest.h>
#include "Individual.h"
#include "diversity_functions.h"

using ::testing::EmptyTestEventListener;
using ::testing::InitGoogleTest;
using ::testing::Test;
using ::testing::TestCase;
using ::testing::TestEventListeners;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::UnitTest;

TEST(DiversityMeasuresTest, MeanDiversity){
    std::vector<Individual> population;

    std::vector<std::vector<double> > x = {
        {1,2,1},
        {4,3,3},
        {5,4,5},
        {6,6,9},
        {8,7,1}
    };

    for (int i = 0; i < x.size(); i++) {
        Individual new_ind; //Create constructor for individual
        new_ind.atts  = x.at(i);
        population.push_back(new_ind);
    }

    EXPECT_FLOAT_EQ(meanDiversity(population, 3, -10, 10), 5.92);
}

TEST(DiversityMeasuresTest, MaxDiversity){
    std::vector<Individual> population;

    std::vector<std::vector<double> > x = {
        {1,2,1},
        {4,3,3},
        {5,4,5},
        {6,6,9},
        {8,7,1}
    };

    for (int i = 0; i < x.size(); i++) {
        Individual new_ind; //Create constructor for individual
        new_ind.atts  = x.at(i);
        population.push_back(new_ind);
    }

    EXPECT_FLOAT_EQ(maxDiversity(population, 3, -10, 10), 0.251481481);
}

TEST(DiversityMeasuresTest, SwarmRadius){
    std::vector<Individual> population;

    std::vector<std::vector<double> > x = {
        {1,2,1},
        {4,3,3},
        {5,4,5},
        {6,6,9},
        {8,7,1}
    };

    for (int i = 0; i < x.size(); i++) {
        Individual new_ind; //Create constructor for individual
        new_ind.atts  = x.at(i);
        population.push_back(new_ind);
    }

    EXPECT_FLOAT_EQ(swarmRadiusDiversity(population, 3, -10, 10), 9.423375192);
}

TEST(DiversityMeasuresTest, DGEATest){
    std::vector<Individual> population;

    std::vector<std::vector<double> > x = {
        {1,2,1},
        {4,3,3},
        {5,4,5},
        {6,6,9},
        {8,7,1}
    };

    for (int i = 0; i < x.size(); i++) {
        Individual new_ind; //Create constructor for individual
        new_ind.atts  = x.at(i);
        population.push_back(new_ind);
    }

    EXPECT_FLOAT_EQ(dgeaDiversity(population, 3, -10, 10), 0.267728098);
}

TEST(DiversityMeasuresTest, EntropyMeasure){
    std::vector<Individual> population;

    std::vector<std::vector<double> > x = {
        {1,2,1},
        {4,3,3},
        {5,4,5},
        {6,6,9},
        {8,7,1}
    };

    for (int i = 0; i < x.size(); i++) {
        Individual new_ind; //Create constructor for individual
        new_ind.atts  = x.at(i);
        population.push_back(new_ind);
    }

    EXPECT_FLOAT_EQ(entropyDiversity(population, 3, -10, 10), 4.159146681);
}

TEST(DiversityMeasuresTest, AverageDistanceAroundCentre){
    std::vector<Individual> population;

    std::vector<std::vector<double> > x = {
        {1,2,1},
        {4,3,3},
        {5,4,5},
        {6,6,9},
        {8,7,1}
    };

    for (int i = 0; i < x.size(); i++) {
        Individual new_ind; //Create constructor for individual
        new_ind.atts  = x.at(i);
        population.push_back(new_ind);
    }

    EXPECT_FLOAT_EQ(averageDistanceAroundCentreDiversity(population, 3, -10, 10), 3.786247079);
}

int main(int argc, char **argv) {
  InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

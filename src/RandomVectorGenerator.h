#include <vector>
#include <random>

#pragma once
class RandomVectorGenerator {
    std::mt19937 randomGenerator;
    std::normal_distribution<double> normalDistribution;
    std::uniform_distribution<double> uniformDistribution;
public:
    RandomVectorGenerator(double seed);
    RandomVectorGenerator();
    std::vector<double> getRandomVector(double min, double max, int dimensions);
    std::vector<double> getUniformRandomVector(double min, double max, int dimensions);
    double getNormalRandom();
    double getUniformRandom();
    double getNormalRandomWithinRange(double min, double max);
    double getUniformRandomWithinRange(double min, double max);
};

#include <vector>
#include <random>
#include <algorithm>
#include "RandomVectorGenerator.h"

using namespace std;

RandomVectorGenerator::RandomVectorGenerator(double seed) {
    randomGenerator = std::mt19937(seed);
}

RandomVectorGenerator::RandomVectorGenerator() {
    randomGenerator = std::mt19937();
}

std::vector<double> RandomVectorGenerator::getRandomVector(double min, double max, int dimensions) {
    std::vector<double> randomVector;
    for (int i = 0; i < dimensions; i++) {
        randomVector.push_back(getNormalRandomWithinRange(min, max));
    }
    return randomVector;
}

std::vector<double> RandomVectorGenerator::getUniformRandomVector(double min, double max, int dimensions) {
    std::vector<double> randomVector;
    for (int i = 0; i < dimensions; i++) {
        randomVector.push_back(getUniformRandomWithinRange(min, max));
    }
    return randomVector;
}

double RandomVectorGenerator::getNormalRandom() {
    return normalDistribution(randomGenerator);
}

double RandomVectorGenerator::getUniformRandom() {
    return uniformDistribution(randomGenerator);
}

double RandomVectorGenerator::getUniformRandomWithinRange(double min, double max) {
    return min + (max - min) * getUniformRandom();
}


double RandomVectorGenerator::getNormalRandomWithinRange(double min, double max) {
    return min + (max - min) * getNormalRandom();
}

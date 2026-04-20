#ifndef NOVOLEITOR_INSTANCE_H
#define NOVOLEITOR_INSTANCE_H

#include <limits>
#include <string>
#include <vector>

typedef struct Instance {
    std::string name;
    int dimension;
    bool explicitCoord = false;

    std::vector<double> distMatrix;
    std::vector<double> xCoord;
    std::vector<double> yCoord;

    [[nodiscard]] double getDistance(const int i, const int j) const {
        return distMatrix[i * dimension + j];
    }

    static constexpr double INF = std::numeric_limits<double>::infinity();
} Instance;

#endif

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "instance.h"
#include "parser.h"

namespace dist {
    double euc(const double x1, const double y1, const double x2, const double y2) {
        const double dx = x2 - x1;
        const double dy = y2 - y1;
        return sqrt(dx * dx + dy * dy);
    }

    double pseudoEuc(const double x1, const double y1, const double x2, const double y2) {
        const double dx = x2 - x1;
        const double dy = y2 - y1;

        const double rij = std::sqrt((dx * dx + dy * dy) / 10.0);
        const double tij = std::round(rij);

        return (tij < rij) ? (tij + 1.0) : tij;
    }

    // Pareceu meio esoterico no original, mas eh para converter para radianos
    double toRadians(const double coord) {
        constexpr double TSP_PI = 3.141592; // Para match com a precisao historica de pi que usaram
        const int deg = static_cast<int>(coord);
        const double min = coord - deg;
        return TSP_PI * (deg + 5.0 * min / 3.0) / 180.0;
    }

    double geo(const double lat1, const double lon1, const double lat2, const double lon2) {
        constexpr double EARTH_RADIUS = 6378.388; // kilometres, idealised.

        const double q1 = std::cos(lon1 - lon2);
        const double q2 = std::cos(lat1 - lat2);
        const double q3 = std::cos(lat1 + lat2);

        const int distance = static_cast<int>(
            EARTH_RADIUS * std::acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0
        );
        return static_cast<double>(distance);
    }
}

Instance parseTSP(const std::string &filepath) {
    Instance instance;

    instance.name = std::filesystem::path(filepath)
            .stem()
            .string();

    // ReSharper disable once CppTooWideScopeInitStatement
    std::ifstream file(filepath);

    if (!file) {
        throw std::runtime_error("Erro ao abrir arquivo " + filepath);
    }

    std::string token;
    std::string ewType;
    std::string ewFormat;

    auto readValue = [&file]() {
        std::string value;
        file >> value;

        if (value == ":") {
            file >> value;
        } else if (!value.empty() && value.front() == ':') {
            value = value.substr(1);
        }

        return value;
    };

    bool isExplicit = false;

    while (file >> token) {
        if (token == "DIMENSION" || token == "DIMENSION:") {
            instance.dimension = std::stoi(readValue());
        } else if (token == "EDGE_WEIGHT_TYPE" || token == "EDGE_WEIGHT_TYPE:") {
            ewType = readValue();
        } else if (token == "EDGE_WEIGHT_FORMAT" || token == "EDGE_WEIGHT_FORMAT:") {
            ewFormat = readValue();
        } else if (token == "EDGE_WEIGHT_SECTION") {
            isExplicit = true;
            break;
        } else if (token == "NODE_COORD_SECTION") {
            isExplicit = false;
            break;
        } else if (token == "EOF") {
            break;
        }
    }

    if (instance.dimension <= 0) {
        throw std::runtime_error("Erro ao processar as dimensoes.");
    }

    const int dim = instance.dimension;

    instance.distMatrix.assign(dim * dim, Instance::INF);

    if (isExplicit) {
        if (ewFormat == "FULL_MATRIX") {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    double value;
                    file >> value;

                    if (i != j) instance.distMatrix[i * dim + j] = value;
                }
            }
        } else if (ewFormat == "UPPER_ROW") {
            for (int i = 0; i < dim; i++) {
                for (int j = i + 1; j < dim; j++) {
                    double value;
                    file >> value;
                    instance.distMatrix[i * dim + j] = value;
                    instance.distMatrix[j * dim + i] = value; // Propriedade simetrica
                }
            }
        } else if (ewFormat == "LOWER_ROW") {
            for (int i = 1; i < dim; i++) {
                for (int j = 0; j < i; j++) {
                    double value;
                    file >> value;
                    instance.distMatrix[i * dim + j] = value;
                    instance.distMatrix[j * dim + i] = value;
                }
            }
        } else if (ewFormat == "UPPER_DIAG_ROW") {
            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    double value;
                    file >> value;
                    if (i != j) {
                        instance.distMatrix[i * dim + j] = value;
                        instance.distMatrix[j * dim + i] = value;
                    }
                }
            }
        } else if (ewFormat == "LOWER_DIAG_ROW") {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j <= i; j++) {
                    double value;
                    file >> value;
                    if (i != j) {
                        instance.distMatrix[i * dim + j] = value;
                        instance.distMatrix[j * dim + i] = value;
                    }
                }
            }
        } else if (ewFormat == "UPPER_COL") {
            for (int j = 1; j < dim; j++) {
                for (int i = 0; i < j; i++) {
                    double value;
                    file >> value;
                    instance.distMatrix[i * dim + j] = value;
                    instance.distMatrix[j * dim + i] = value;
                }
            }
        } else if (ewFormat == "LOWER_COL") {
            for (int j = 0; j < dim; j++) {
                for (int i = j + 1; i < dim; i++) {
                    double value;
                    file >> value;
                    instance.distMatrix[i * dim + j] = value;
                    instance.distMatrix[j * dim + i] = value;
                }
            }
        } else if (ewFormat == "LOWER_DIAG_COL") {
            for (int j = 0; j < dim; j++) {
                for (int i = j; i < dim; i++) {
                    double value;
                    file >> value;
                    if (i != j) {
                        instance.distMatrix[i * dim + j] = value;
                        instance.distMatrix[j * dim + i] = value;
                    }
                }
            }
        } else {
            throw std::runtime_error("Tipo de edge weight format nao suportado: " + ewFormat);
        }
    } else {
        instance.explicitCoord = true;
        instance.xCoord.resize(dim);
        instance.yCoord.resize(dim);

        for (int i = 0; i < dim; i++) {
            int id;
            file >> id >> instance.xCoord[i] >> instance.yCoord[i];
        }

        if (ewType == "GEO") {
            std::vector<double> lats(dim), lons(dim);

            for (int i = 0; i < dim; i++) {
                lats[i] = dist::toRadians(instance.xCoord[i]);
                lons[i] = dist::toRadians(instance.yCoord[i]);
            }

            for (int i = 0; i < dim; i++) {
                for (int j = i + 1; j < dim; j++) {
                    double d = dist::geo(lats[i], lons[i], lats[j], lons[j]);
                    instance.distMatrix[i * dim + j] = d;
                    instance.distMatrix[j * dim + i] = d;
                }
            }
        } else {
            for (int i = 0; i < dim; i++) {
                for (int j = i + 1; j < dim; j++) {
                    double d = 0.0;
                    if (ewType == "EUC_2D") {
                        d = std::floor(
                            dist::euc(instance.xCoord[i], instance.yCoord[i], instance.xCoord[j], instance.yCoord[j])
                            + 0.5
                        );
                    } else if (ewType == "CEIL_2D") {
                        d = std::ceil(
                            dist::euc(instance.xCoord[i], instance.yCoord[i], instance.xCoord[j], instance.yCoord[j])
                        );
                    } else if (ewType == "ATT") {
                        d = dist::pseudoEuc(
                            instance.xCoord[i], instance.yCoord[i], instance.xCoord[j], instance.yCoord[j]
                        );
                    } else {
                        throw std::runtime_error("erro ao processar o edge type, nao suportado: " + ewType);
                    }

                    instance.distMatrix[i * dim + j] = d;
                    instance.distMatrix[j * dim + i] = d;
                }
            }
        }
    }

    return instance;
}

void printMatrix(const Instance &instance) {
    const int dim = instance.dimension;

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            const double dist = instance.distMatrix[i * dim + j];

            if (dist == Instance::INF) {
                std::cout << std::setw(8) << "INF";
            } else {
                std::cout << std::setw(8) << dist;
            }
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}

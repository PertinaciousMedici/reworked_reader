#ifndef NOVOLEITOR_PARSER_H
#define NOVOLEITOR_PARSER_H

#include "instance.h"

Instance parseTSP(const std::string &filepath);
void printMatrix(const Instance &instance);

#endif

#ifndef CLASSIC_OPENMP_2D_SAP_H
#define CLASSIC_OPENMP_2D_SAP_H

#include <string>
#include <vector>

int classic_openmp_2d_sap(const std::vector<std::string> &x, const std::vector<std::string> &y, const std::vector<std::string> &trans_x, const std::vector<std::string> &trans_y, double &preprocessing_time, double &computation_time);

#endif
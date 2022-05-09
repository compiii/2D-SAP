#ifndef CLASSIC_MPI_2D_SAP_H
#define CLASSIC_MPI_2D_SAP_H

#include <string>
#include <vector>
#include <cmath>
#include "../1D/datatype.h"

int classic_mpi_2d_sap(const std::vector<std::string> &x, const std::vector<std::string> &y, const std::vector<std::string> &trans_x, const std::vector<std::string> &trans_y, const int &rank, const int &max_process, double &preprocessing_time, double &pre_comm_time, double &main_time, double &main_comm_time);
bound get_bound(const int &m, const int &rank, const int &max_process);
void bcast(int **&tab, bound &bd, const int &b, const int &s, const int &m, const int &rank, const int &max_process);
#endif
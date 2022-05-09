#ifndef CLASSIC_PAR_EDIT_DISTANCE
#define CLASSIC_PAR_EDIT_DISTANCE

#include "datatype.h"
#include <mpi.h>

#include <string>
#include <vector>

int classic_par_edit_distance(const std::string &a, const std::string &b, const int &rank, const int &max_process, double &computation_time, double &communication_time);
void communicate_block_data(int **&scores, const int &n, const int type, const block_data &src, const block_data dest, const int &is_line, const int &is_column, const int &is_oblique);
void receive_block_data(const block_data &bd, const MPI_Datatype &mpi_vector, int **&scores, const unsigned long &absc, const unsigned long &ord, const int &tag);
void send_block_data(const block_data &bd, const MPI_Datatype &mpi_vector, int **&scores, const unsigned long &absc, const unsigned long &ord, const int &tag);
#endif // CLASSIC_PAR_EDIT_DISTANCE
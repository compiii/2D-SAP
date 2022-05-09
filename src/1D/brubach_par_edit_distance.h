#ifndef BRUBACH_PAR_EDIT_DISTANCE
#define BRUBACH_PAR_EDIT_DISTANCE

#include "datatype.h"
#include <mpi.h>

#include <string>
#include <vector>

int brubach_par_edit_distance(const std::string &str_1, const std::string &str_2, const int &rank, const int &max_process, double &computation_time, double &communication_time, double &table_build_duration, double &efficiency);
//void communicate_block_data(int **&scores, const int &n, const int type, const block_data &src, const block_data dest, const int &is_line, const int &is_column, const int &is_oblique);
//void receive_block_data(const block_data &bd, const MPI_Datatype &mpi_vector, int **&scores, const unsigned long &absc, const unsigned long &ord, const int &tag);
//void send_block_data(const block_data &bd, const MPI_Datatype &mpi_vector, int **&scores, const unsigned long &absc, const unsigned long &ord, const int &tag);
#endif // BRUBACH_PAR_EDIT_DISTANCE
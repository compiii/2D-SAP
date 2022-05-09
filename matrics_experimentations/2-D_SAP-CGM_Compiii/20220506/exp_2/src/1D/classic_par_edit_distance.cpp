#include "classic_par_edit_distance.h"

#include "datatype.h"
#include "partitioning_classic.h"
#include "../common/handy_chrono.h"
#include "classic_seq_edit_distance.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

int classic_par_edit_distance(const string &str_1, const string &str_2, const int &rank, const int &max_process, double &computation_time, double &communication_time)
{
    int max_diag = 0;
    vector<block> tab_block = partitioning_classic::build(str_1, str_2, rank, max_process, max_diag);
    int max_eval_block = tab_block.size();
    int m = str_1.size();
    int n = str_2.size();

    int **scores = new int *[m + 1];
    scores[0] = new int[(m + 1) * (n + 1)];
    for (int i = 1; i < m + 1; i++)
    {
        scores[i] = scores[0] + i * (n + 1);
    }
    chrono_clock start, end;

    start = now();
    for (auto path : tab_block)
    {
        if (path.data.core_data.first_bound.begin == 0)
        {
            bound b = path.data.core_data.second_bound;
            for (int j = b.begin; j <= b.end; j++)
                scores[0][j] = j;
        }

        if (path.data.core_data.second_bound.begin == 0)
        {
            bound b = path.data.core_data.first_bound;
            for (int i = b.begin; i <= b.end; i++)
                scores[i][0] = i;
        }
    }
    end = now();
    computation_time += higher_precision_duration(start, end);

    int i = 0, last = 0;
    block bk;
    for (int d = 1; (d <= max_diag) && (i != max_eval_block); d++)
    {
        bk = tab_block[i];
        while (bk.data.core_data.diag == d && (i != max_eval_block))
        {
            if (d == max_diag)
                last = 1;
            start = now();
            for (auto path : bk.need)
            {
                if (path.is_row == 1 && path.is_column == 1)
                {
                    communicate_block_data(scores, n, 0, bk.data, path, 1, 0, 0);
                    communicate_block_data(scores, n, 0, bk.data, path, 0, 1, 0);
                }
                else if (path.is_row == 1 && path.is_oblique == 1)
                {
                    communicate_block_data(scores, n, 0, bk.data, path, 1, 0, 0);
                }
                else if (path.is_column == 1 && path.is_oblique == 1)
                {
                    communicate_block_data(scores, n, 0, bk.data, path, 0, 1, 0);
                }
                else
                {
                    communicate_block_data(scores, n, 0, bk.data, path, path.is_row, path.is_column, path.is_oblique);
                }
            }
            end = now();
            communication_time += higher_precision_duration(start, end);

            start = now();
            classic_seq_edit_distance(str_1, str_2, bk, scores);
            end = now();
            computation_time += higher_precision_duration(start, end);

            start = now();
            for (auto path : bk.depend)
            {
                if (path.is_row == 1 && path.is_column == 1)
                {
                    communicate_block_data(scores, n, 1, bk.data, path, 1, 0, 0);
                    communicate_block_data(scores, n, 1, bk.data, path, 0, 1, 0);
                }
                else if (path.is_row == 1 && path.is_oblique == 1)
                {
                    communicate_block_data(scores, n, 1, bk.data, path, 1, 0, 0);
                }
                else if (path.is_column == 1 && path.is_oblique == 1)
                {
                    communicate_block_data(scores, n, 1, bk.data, path, 0, 1, 0);
                }
                else
                {
                    communicate_block_data(scores, n, 1, bk.data, path, path.is_row, path.is_column, path.is_oblique);
                }
            }
            end = now();
            communication_time += higher_precision_duration(start, end);
            i++;
            if (i != max_eval_block)
                bk = tab_block[i];
        }
    }

    if (last)
        return scores[m][n];
    else
        return -1;
}

void communicate_block_data(int **&scores, const int &n, const int type, const block_data &src, const block_data dest, const int &is_row, const int &is_column, const int &is_oblique)
{
    int count, blocklengths, stride;
    unsigned long absc, ord;
    block_data bd = (type == 1 ? src : dest);
    MPI_Datatype mpi_vector;

    if (is_column == 1)
    {
        count = 1;
        blocklengths = bd.core_data.d.column;
        stride = 1;
        MPI_Type_vector(count, blocklengths, stride, MPI_INT, &mpi_vector);
        MPI_Type_commit(&mpi_vector);
        absc = bd.core_data.first_bound.end;
        ord = bd.core_data.second_bound.begin;
    }
    else if (is_row == 1)
    {
        count = bd.core_data.d.row;
        blocklengths = 1;
        stride = n + 1;
        MPI_Type_vector(count, blocklengths, stride, MPI_INT, &mpi_vector);
        MPI_Type_commit(&mpi_vector);

        absc = bd.core_data.first_bound.begin;
        ord = bd.core_data.second_bound.end;
    }
    else if (is_oblique == 1)
    {
        count = 1;
        blocklengths = 1;
        stride = 1;
        MPI_Type_vector(count, blocklengths, stride, MPI_INT, &mpi_vector);
        MPI_Type_commit(&mpi_vector);
        absc = bd.core_data.first_bound.end;
        ord = bd.core_data.second_bound.end;
    }

    bd = dest;

    if (type == 0)
        receive_block_data(bd, mpi_vector, scores, absc, ord, dest.core_data.id);
    else
        send_block_data(bd, mpi_vector, scores, absc, ord, src.core_data.id);
    MPI_Type_free(&mpi_vector);
}

void receive_block_data(const block_data &bd, const MPI_Datatype &mpi_vector, int **&scores, const unsigned long &absc, const unsigned long &ord, const int &tag)
{
    MPI_Recv(&scores[absc][ord], 1, mpi_vector, bd.rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void send_block_data(const block_data &bd, const MPI_Datatype &mpi_vector, int **&scores, const unsigned long &absc, const unsigned long &ord, const int &tag)
{
    MPI_Request r;
    MPI_Isend(&scores[absc][ord], 1, mpi_vector, bd.rank, tag, MPI_COMM_WORLD, &r);
}
#include "brubach_par_edit_distance.h"

#include "datatype.h"
#include "partitioning_4r.h"
#include "../common/handy_chrono.h"
#include "brubach_seq_edit_distance.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

int brubach_par_edit_distance(const std::string &str_1, const std::string &str_2, const int &rank, const int &max_process, double &computation_time, double &communication_time, double &table_build_duration, double &efficiency)
{
    int max_diag = 0;
    vector<block> tab_block = partitioning_4r::build(str_1, str_2, rank, max_process, max_diag);
    int max_eval_block = tab_block.size();

    int m = str_1.size();
    int n = str_2.size();

    int t = log10(max(m, n)) * log2(max(m, n));

    brubach_lookup_table_t lookup_table;
    row_column initial_values, rs;
    string a, b;
    vector<vector<int>> reconstitute_solution;

    chrono_clock start, end;

    int i = 0;
    block bk;
    int *send_recv;
    for (int d = 1; (d <= max_diag) && (i != max_eval_block); d++)
    {
        bk = tab_block[i];
        while (bk.data.core_data.diag == d && (i != max_eval_block))
        {
            if (!bk.need.empty())
            {
                for (auto path : bk.need)
                {
                    if (path.is_row == 1)
                    {
                        // receive row
                        start = now();
                        send_recv = new int[(bk.data.core_data.first_bound.end - bk.data.core_data.first_bound.begin)];
                        initial_values.row = vector<int>((bk.data.core_data.first_bound.end - bk.data.core_data.first_bound.begin));
                        MPI_Recv(send_recv, initial_values.row.size(), MPI_INT, path.rank, path.core_data.id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        initial_values.row.assign(send_recv, send_recv + initial_values.row.size());
                        end = now();
                        communication_time += higher_precision_duration(start, end);

                        if (2 * bk.data.core_data.c.i == bk.data.core_data.c.j)
                        {
                            initial_values.column = vector<int>(bk.data.core_data.second_bound.end - bk.data.core_data.second_bound.begin + 1, 1);
                        }
                    }
                    else if (path.is_column == 1)
                    {
                        // receive column
                        start = now();
                        send_recv = new int[(bk.data.core_data.second_bound.end - bk.data.core_data.second_bound.begin)];
                        initial_values.column = vector<int>((bk.data.core_data.second_bound.end - bk.data.core_data.second_bound.begin));
                        MPI_Recv(send_recv, initial_values.column.size(), MPI_INT, path.rank, path.core_data.id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        initial_values.column.assign(send_recv, send_recv + initial_values.column.size());
                        end = now();
                        communication_time += higher_precision_duration(start, end);
                        if (bk.data.core_data.c.i == 1 && bk.data.core_data.c.j - 1 <= max_process)
                        {
                            initial_values.row = vector<int>(bk.data.core_data.first_bound.end - bk.data.core_data.first_bound.begin + 1, 1);
                        }
                    }
                }
            }
            else
            {
                if (bk.data.core_data.id == 0)
                {
                    initial_values = {vector<int>(bk.data.core_data.first_bound.end, 1), vector<int>(bk.data.core_data.second_bound.end, 1)};
                }
                else
                {
                    initial_values.column = vector<int>(bk.data.core_data.second_bound.end - bk.data.core_data.second_bound.begin + 1, 1);
                }
            }

            start = now();
            a = str_1.substr(bk.data.core_data.first_bound.begin, bk.data.core_data.first_bound.end - bk.data.core_data.first_bound.begin);
            b = str_2.substr(bk.data.core_data.second_bound.begin, bk.data.core_data.second_bound.end - bk.data.core_data.second_bound.begin);
            rs = brubach_seq_edit_distance(a, b, initial_values, lookup_table, t, table_build_duration, efficiency);
            end = now();
            computation_time += higher_precision_duration(start, end);

            if (bk.data.core_data.c.j - 1 == max_diag)
            {
                reconstitute_solution.push_back(rs.row);
            }

            if (!bk.depend.empty())
            {
                int sent_row = 0;
                for (auto path : bk.depend)
                {
                    if (path.is_row == 1)
                    {
                        sent_row = 1;
                        start = now();
                        MPI_Request r;
                        send_recv = new int[rs.row.size()];
                        copy(rs.row.begin(), rs.row.end(), send_recv);
                        MPI_Isend(send_recv, rs.row.size(), MPI_INT, path.rank, bk.data.core_data.id, MPI_COMM_WORLD, &r);
                        end = now();
                        communication_time += higher_precision_duration(start, end);
                    }

                    else if (path.is_column == 1)
                    {
                        // send column
                        start = now();
                        MPI_Request r;
                        send_recv = new int[rs.column.size()];
                        copy(rs.column.begin(), rs.column.end(), send_recv);
                        MPI_Isend(send_recv, rs.column.size(), MPI_INT, path.rank, bk.data.core_data.id, MPI_COMM_WORLD, &r);
                        end = now();
                        communication_time += higher_precision_duration(start, end);
                    }
                }

                if (sent_row == 0)
                {
                    initial_values.row = rs.row;
                }
            }
            else
            {
                initial_values.row = rs.row;
            }

            i++;
            if (i != max_eval_block)
                bk = tab_block[i];
        }
    }

    // Reconstitute solution
    int sum = 0;
    start = now();
    for (auto row : reconstitute_solution)
    {
        for (unsigned int i = 1; i <= row.size(); ++i)
            sum += row[i - 1];
    }
    end = now();
    computation_time += higher_precision_duration(start, end);

    if (rank != 0)
    {
        start = now();
        MPI_Request r;
        MPI_Isend(&sum, 1, MPI_INT, 0, 237, MPI_COMM_WORLD, &r);
        end = now();
        communication_time += higher_precision_duration(start, end);
        return -1;
    }
    else
    {
        int tmp = 0;
        for (int processor = 1; processor < max_process; processor++)
        {
            start = now();
            MPI_Recv(&tmp, 1, MPI_INT, processor, 237, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            end = now();
            communication_time += higher_precision_duration(start, end);
            sum += tmp;
        }
        return n + sum;
    }
}
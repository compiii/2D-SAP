#include "classic_mpi_2d_sap.h"
#include "indice.h"
#include "../common/editing_op.h"
#include "../common/handy_chrono.h"
#include "../1D/classic_seq_edit_distance.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include <mpi.h>

using namespace std;

int classic_mpi_2d_sap(const vector<string> &x, const vector<string> &y, const vector<string> &trans_x, const vector<string> &trans_y, const int &rank, const int &max_process, double &preprocessing_time, double &pre_comm_time, double &main_time, double &main_comm_time)
{
    int m1 = x.size();
    int n1 = x[1].size();
    int m2 = y.size();
    int n2 = y[1].size();
    indice ind(m1, n1, m2, n2);
    // indice ind_bis(n1, m1, n2, m2);

    // vector<vector<int>> dr(m1 + 1, vector<int>(n1 + 1));
    // vector<vector<int>> dc(m1 + 1, vector<int>(n1 + 1));

    // vector<vector<int>> ir(m2 + 1, vector<int>(n2 + 1));
    // vector<vector<int>> ic(m2 + 1, vector<int>(n2 + 1));

    // vector<vector<vector<vector<int>>>> R(m1 + 1, vector<vector<vector<int>>>(n1 + 1, vector<vector<int>>(m2 + 1, vector<int>(n2 + 1))));
    // vector<vector<vector<vector<int>>>> C(n1 + 1, vector<vector<vector<int>>>(m1 + 1, vector<vector<int>>(n2 + 1, vector<int>(m2 + 1))));

    // vector<vector<vector<vector<int>>>> T(m1 + 1, vector<vector<vector<int>>>(n1 + 1, vector<vector<int>>(m2 + 1, vector<int>(n2 + 1))));

    int **dr = new int *[m1 + 1];
    dr[0] = new int[(m1 + 1) * (n1 + 1)];
    for (int i = 1; i < m1 + 1; i++)
    {
        dr[i] = dr[0] + i * (n1 + 1);
    }

    int **dc = new int *[m1 + 1];
    dc[0] = new int[(m1 + 1) * (n1 + 1)];
    for (int i = 1; i < m1 + 1; i++)
    {
        dc[i] = dc[0] + i * (n1 + 1);
    }

    int **ir = new int *[m2 + 1];
    ir[0] = new int[(m2 + 1) * (n2 + 1)];
    for (int i = 1; i < m2 + 1; i++)
    {
        ir[i] = ir[0] + i * (n2 + 1);
    }

    int **ic = new int *[m2 + 1];
    ic[0] = new int[(m2 + 1) * (n2 + 1)];
    for (int i = 1; i < m2 + 1; i++)
    {
        ic[i] = ic[0] + i * (n2 + 1);
    }

    int **R = new int *[m1 + 1];
    R[0] = new int[(m1 + 1) * (n1 + 1) * (m2 + 1) * (n2 + 1)];
    for (int i = 1; i < m1 + 1; i++)
    {
        R[i] = R[0] + i * (n1 + 1) * (m2 + 1) * (n2 + 1);
    }

    /*int **C = new int *[n1 + 1];
    C[0] = new int[(n1 + 1) * (m1 + 1) * (n2 + 1) * (m2 + 1)];
    for (int i = 1; i < n1 + 1; i++)
    {
        C[i] = C[0] + i * (m1 + 1) * (n2 + 1) * (m2 + 1);
    }*/
    int **C = new int *[m1 + 1];
    C[0] = new int[(m1 + 1) * (n1 + 1) * (m2 + 1) * (n2 + 1)];
    for (int i = 1; i < m1 + 1; i++)
    {
        C[i] = C[0] + i * (n1 + 1) * (m2 + 1) * (n2 + 1);
    }

    int **T = new int *[m1 + 1];
    T[0] = new int[(m1 + 1) * (n1 + 1) * (m2 + 1) * (n2 + 1)];
    for (int i = 1; i < m1 + 1; i++)
    {
        T[i] = T[0] + i * (n1 + 1) * (m2 + 1) * (n2 + 1);
    }

    // vector<vector<int>> scores1(m1 + 1, vector<int>(n1 + 1));
    int z = max({m1, m2, n1, n2});
    int **scores = new int *[z + 1];
    scores[0] = new int[(z + 1) * (z + 1)];
    for (int i = 1; i < z + 1; i++)
    {
        scores[i] = scores[0] + i * (z + 1);
    }

    scores[0][0] = 0;
    for (int i = 1; i <= z; ++i)
        scores[i][0] = scores[0][0] + del();

    for (int j = 1; j <= z; ++j)
        scores[0][j] = scores[0][0] + ins();

    // vector<vector<int>> scores2(m2 + 1, vector<int>(n2 + 1));
    /*int **scores2 = new int *[m2 + z];
    scores2[0] = new int[(m2 + z) * (n2 + z)];
    for (int i = 1; i < m2 + z; i++)
    {
        scores2[i] = scores2[0] + i * (n2 + z);
    }
    scores2[0][0] = 0;
    for (int i = 1; i <= m2 + z - 1; ++i)
        scores2[i][0] = scores2[0][0] + del();

    for (int j = 1; j <= n2 + z - 1; ++j)
        scores2[0][j] = scores2[0][0] + ins();*/

    chrono_clock start, end, tmp_clock;
    int bi, ei, bj, ej;
    bound b1 = get_bound(m1 + 1, rank, max_process);
    bound b2 = get_bound(m2 + 1, rank, max_process);
    bound b3;

    start = now();
    // if (rank == 0)
    //     cout << m1 + 1 << endl;
    // cout << rank << "\t" << b1.begin << " " << b1.end << endl;
    for (int i = (b1.begin == 0 ? 1 : b1.begin); i <= b1.end; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            dr[i][j] = j * del();
        }
    }

    for (int i = (b1.begin == 0 ? 1 : b1.begin); i <= b1.end; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            dc[i][j] = i * del();
        }
    }

    for (int i = (b2.begin == 0 ? 1 : b2.begin); i <= b2.end; i++)
    {
        for (int j = 1; j <= n2; j++)
        {
            ir[i][j] = j * ins();
        }
    }

    for (int i = (b2.begin == 0 ? 1 : b2.begin); i <= b2.end; i++)
    {
        for (int j = 1; j <= n2; j++)
        {
            ic[i][j] = i * ins();
        }
    }

    bi = (rank == 0 ? 1 : b1.begin);
    ei = (rank == max_process - 1 ? b1.end : b1.end + 1);

    for (int i = bi; i < ei; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    if (rank == 0 && 1 == 2)
                        cout << m1 + 1 << "\t" << x[i].substr(1, j) << "-" << y[k].substr(1, l) << endl;
                    R[i][ind.get_indice_cube(j, k, l)] = classic_seq_edit_distance(x[i].substr(1, j), y[k].substr(1, l), scores);
                }
            }
        }
    }

    /*for (int i = 1; i < m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    cout << R[i][ind.get_indice_cube(j, k, l)] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }*/

    for (int i = bi; i < ei; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    // cout << trans_x[i].substr(1, j) << "-" << trans_y[k].substr(1, l) << endl;
                    C[i][ind.get_indice_cube(j, k, l)] = classic_seq_edit_distance(trans_x[j].substr(1, i), trans_y[l].substr(1, k), scores);
                }
            }
        }
    }
    /*for (int i = 1; i < n1; i++)
    {
        for (int j = 1; j <= m1; j++)
        {
            for (int k = 1; k < n2; k++)
            {
                for (int l = 1; l <= m2; l++)
                {
                    cout << C[i][ind_bis.get_indice_cube(j, k, l)] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }*/
    tmp_clock = now();
    bcast(dr, b1, n1 + 1, n1 + 1, m1 + 1, rank, max_process);
    bcast(dc, b1, n1 + 1, n1 + 1, m1 + 1, rank, max_process);
    bcast(ir, b2, n2 + 1, n2 + 1, m2 + 1, rank, max_process);
    bcast(ic, b2, n2 + 1, n2 + 1, m2 + 1, rank, max_process);
    // bcast(R, b1, (n1 + 1) * (m2 + 1) * (n2 + 1), (n1 + 1) * (m2 + 1) * (n2 + 1), m1 + 1, rank, max_process);
    // bcast(C, b1, (n1 + 1) * (m2 + 1) * (n2 + 1), (n1 + 1) * (m2 + 1) * (n2 + 1), m1 + 1, rank, max_process);
    end = now();
    pre_comm_time += higher_precision_duration(tmp_clock, end);
    end = now();
    preprocessing_time += higher_precision_duration(start, end);

    start = now();
    bi = (rank == 0 ? 1 : b1.begin);
    ei = (rank == max_process - 1 ? b1.end : b1.end + 1);
    b3 = get_bound(n1 + 1, rank, max_process);
    bj = (rank == 0 ? 1 : b3.begin);
    ej = b3.end; //(rank == max_process - 1 ? b3.end : b3.end + 1);

    MPI_Datatype mpi_vector;
    int count = 1;
    int blocklength = (m2 + 1) * (n2 + 1);
    int stride = 1;
    MPI_Type_vector(count, blocklength, stride, MPI_INT, &mpi_vector);
    MPI_Type_commit(&mpi_vector);
    MPI_Request r;
    bound tmp_bd;
    if (rank > 0)
        tmp_bd = get_bound(n1 + 1, rank - 1, max_process);

    T[0][0] = 0;
    for (int i = 0; i < m1; i++)
    {
        for (int j = 0; j <= ej; j++)
        {
            for (int k = 0; k < m2; k++)
            {
                for (int l = 0; l <= n2; l++)
                {
                    if ((i == 0 || j == 0) && (k > 0 && l > 0))
                    {
                        T[0][ind.get_indice_cube(j, k, l)] = (k + 1) * (l + 1);
                        T[i][ind.get_indice_cube(0, k, l)] = (k + 1) * (l + 1);
                    }
                    else if ((i > 0 && j > 0) && (k == 0 || l == 0))
                    {
                        T[i][ind.get_indice_cube(j, 0, l)] = (i + 1) * (j + 1);
                        T[i][ind.get_indice_cube(j, k, 0)] = (i + 1) * (j + 1);
                    }
                    else if ((i == 0 && k == 0) || (j == 0 && l == 0))
                    {
                        T[i][ind.get_indice_cube(j, k, l)] = 0;
                    }
                    else
                    {
                        T[i][ind.get_indice_cube(j, k, l)] = -1;
                    }
                }
            }
        }
    }

    for (int i = 1; i < m1; i++)
    {
        // receive data
        tmp_clock = now();
        if (rank > 0)
        {
            /*MPI_Datatype mpi_vector;
            int count = 1;
            int blocklength = (m2 + 1) * (n2 + 1);
            int stride = 1;
            bound tmp_bd = get_bound(n1 + 1, rank - 1, max_process);
            MPI_Type_vector(count, blocklength, stride, MPI_INT, &mpi_vector);
            MPI_Type_commit(&mpi_vector);
            MPI_Recv(&T[i][tmp_bd.end], 1, mpi_vector, rank - 1, rank - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Type_free(&mpi_vector);*/
            MPI_Recv(&T[i][ind.get_indice_cube(tmp_bd.end, 0, 0)], 1, mpi_vector, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&C[i][ind.get_indice_cube(tmp_bd.end, 0, 0)], 1, mpi_vector, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // cout << rank << "\treceived : \tcube = " << i << " plan = " << bj << "-" << ej << "\tmy address = " << ind.get_indice_cube(b3.end, 0, 0) << "\trecv address = " << ind.get_indice_cube(tmp_bd.end, 0, 0) << endl;
        }
        end = now();
        main_comm_time += higher_precision_duration(tmp_clock, end);
        for (int j = bj; j <= ej; j++)
        {
            // computation
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    int a = T[i - 1][ind.get_indice_cube(j, k, l)] + dr[i][j];
                    int b = T[i][ind.get_indice_cube(j - 1, k, l)] + dc[i][j];
                    int c = T[i][ind.get_indice_cube(j, k - 1, l)] + ir[k][l];
                    int d = T[i][ind.get_indice_cube(j, k, l - 1)] + ic[k][l];
                    int e = T[i - 1][ind.get_indice_cube(j, k - 1, l)] + R[i][ind.get_indice_cube(j, k, l)];
                    int f = T[i][ind.get_indice_cube(j - 1, k, l - 1)] + C[i][ind.get_indice_cube(j, k, l)];
                    int g = T[i - 1][ind.get_indice_cube(j - 1, k - 1, l - 1)] + ((i == 1 || k == 1) ? 0 : C[i - 1][ind.get_indice_cube(j, k - 1, l)]) + R[i][ind.get_indice_cube(j, k, l)];
                    int h = T[i - 1][ind.get_indice_cube(j - 1, k - 1, l - 1)] + C[i][ind.get_indice_cube(j, k, l)] + ((j == 1 || l == 1) ? 0 : R[i][ind.get_indice_cube(j - 1, k, l - 1)]);
                    T[i][ind.get_indice_cube(j, k, l)] = max({a, b, c, d, e, f, g, h});
                }
            }
        }
        // send data
        tmp_clock = now();
        if (rank != (max_process - 1))
        {
            /*MPI_Datatype mpi_vector;
            int count = 1;
            int blocklength = (m2 + 1) * (n2 + 1);
            int stride = 1;
            MPI_Request r;
            MPI_Type_vector(count, blocklength, stride, MPI_INT, &mpi_vector);
            MPI_Type_commit(&mpi_vector);
            MPI_Isend(&T[i][b3.end], 1, mpi_vector, rank + 1, rank, MPI_COMM_WORLD, &r);
            MPI_Type_free(&mpi_vector);*/
            MPI_Isend(&T[i][ind.get_indice_cube(b3.end, 0, 0)], 1, mpi_vector, rank + 1, 0, MPI_COMM_WORLD, &r);
            MPI_Isend(&C[i][ind.get_indice_cube(b3.end, 0, 0)], 1, mpi_vector, rank + 1, 1, MPI_COMM_WORLD, &r);
            // cout << rank << "\tsend : \t\tcube = " << i << " plan = " << bj << "-" << ej << "\tmy address = " << ind.get_indice_cube(b3.end, 0, 0) << endl;
        }
        end = now();
        main_comm_time += higher_precision_duration(tmp_clock, end);
    }
    end = now();
    main_time += higher_precision_duration(start, end);
    MPI_Type_free(&mpi_vector);
    if (rank == (max_process - 1)) //&& 1 == 2
        return T[m1 - 1][ind.get_indice_cube(n1, m2 - 1, n2)];
    else
        return 0;
}

bound get_bound(const int &m, const int &rank, const int &max_process)
{
    bound b;
    int option = 1;
    int count = (int)floor(m / (float)max_process);
    if (m % max_process == 0 || ((count * max_process + rank + 1) > m))
        option = 0;
    else
        count++;
    if (rank == 0)
    {
        b.begin = 0;
        b.end = count - 1;
    }
    else
    {
        if (option == 1)
        {
            b.begin = (rank * count);
            b.end = b.begin + count - 1;
        }
        else
        {
            b.begin = 0;
            int id = rank - 1;
            while ((count * max_process + id + 1) > m)
            {
                b.begin += count;
                id--;
            }
            b.begin += (count + 1) * (id + 1);
            b.end = b.begin + count - 1;
        }
    }
    return b;
}

void bcast(int **&tab, bound &bd, const int &b, const int &s, const int &m, const int &rank, const int &max_process)
{
    MPI_Datatype mpi_vector;
    int id;
    int count = bd.length();
    int blocklength = b;
    int stride = s;
    MPI_Request r;
    MPI_Type_vector(count, blocklength, stride, MPI_INT, &mpi_vector);
    MPI_Type_commit(&mpi_vector);
    for (id = 0; id < max_process; id++)
    {
        if (rank != id)
        {
            // cout << rank << " s -> \t" << count << " " << blocklength << " " << s << " " << bd.begin << endl;
            MPI_Isend(&tab[bd.begin][0], 1, mpi_vector, id, 0, MPI_COMM_WORLD, &r);
        }
    }
    for (id = 0; id < max_process; id++)
    {
        if (rank != id)
        {
            bound tmp_bd = get_bound(m, id, max_process);
            // cout << rank << "\t\t\t" << m << " | " << tmp_bd.begin << "-" << tmp_bd.end << endl;
            count = tmp_bd.length();
            // cout << rank << " r -> \t" << count << " " << blocklength << " " << s << " " << tmp_bd.begin << endl;
            MPI_Type_vector(count, blocklength, stride, MPI_INT, &mpi_vector);
            MPI_Type_commit(&mpi_vector);
            MPI_Recv(&tab[tmp_bd.begin][0], 1, mpi_vector, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    MPI_Type_free(&mpi_vector);
}
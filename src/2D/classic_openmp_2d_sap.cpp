#include "classic_seq_2d_sap.h"
#include "indice.h"
#include "../common/editing_op.h"
#include "../common/handy_chrono.h"
#include "../1D/classic_seq_edit_distance.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;

int classic_openmp_2d_sap(const vector<string> &x, const vector<string> &y, const vector<string> &trans_x, const vector<string> &trans_y, double &preprocessing_time, double &computation_time)
{
    int m1 = x.size();
    int n1 = x[1].size();
    int m2 = y.size();
    int n2 = y[1].size();
    indice ind(m1, n1, m2, n2);

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

    chrono_clock start, end, tmp_clock, tmp_thread_clock;
    double tmp_time = 0;
    double tmp_thread_time = 0;
    start = now();

    tmp_clock = now();
#pragma omp parallel for shared(m1, n1, dr) schedule(static)
    for (int c = 0; c < m1 * n1; c++)
    {
        const int i = c / n1 + 1;
        const int j = c % n1 + 1;
        dr[i][j] = j * del();
    }
    /*for (int i = 1; i <= m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            // cout << i << " - " << j << " -> " << omp_get_thread_num() << " tid " << endl;
            dr[i][j] = j * del();
        }
    }*/
    end = now();
    tmp_time = higher_precision_duration(tmp_clock, end);
    // cout << "[dr] = " << tmp_time << "s" << endl;

    tmp_clock = now();
#pragma omp parallel for shared(m1, n1, dc) schedule(static)
    for (int c = 0; c < m1 * n1; c++)
    {
        const int i = c / n1 + 1;
        const int j = c % n1 + 1;
        // printf("(%d, %d)\n", i, j);
        dc[i][j] = i * del();
    }
    /*for (int i = 1; i <= m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            dc[i][j] = i * del();
        }
    }*/
    end = now();
    tmp_time = higher_precision_duration(tmp_clock, end);
    // cout << "[dc] = " << tmp_time << "s" << endl;

    tmp_clock = now();
#pragma omp parallel for shared(m2, n2, ir) schedule(static)
    for (int c = 0; c < m2 * n2; c++)
    {
        const int i = c / n2 + 1;
        const int j = c % n2 + 1;
        // printf("(%d, %d)\n", i, j);
        ir[i][j] = j * ins();
    }
    /*for (int i = 1; i <= m2; i++)
    {
        for (int j = 1; j <= n2; j++)
        {
            ir[i][j] = j * ins();
        }
    }*/
    end = now();
    tmp_time = higher_precision_duration(tmp_clock, end);
    // cout << "[ir] = " << tmp_time << "s" << endl;

    tmp_clock = now();
#pragma omp parallel for shared(m2, n2, ic) schedule(static)
    for (int c = 0; c < m2 * n2; c++)
    {
        const int i = c / n2 + 1;
        const int j = c % n2 + 1;
        // printf("(%d, %d)\n", i, j);
        ic[i][j] = i * ins();
    }
    /*for (int i = 1; i <= m2; i++)
    {
        for (int j = 1; j <= n2; j++)
        {
            ic[i][j] = i * ins();
        }
    }*/
    end = now();
    tmp_time = higher_precision_duration(tmp_clock, end);
    // cout << "[ic] = " << tmp_time << "s" << endl;

    vector<double> tmp_threads_time(omp_get_max_threads());
    vector<int> tmp_threads_load(omp_get_max_threads());
    tmp_clock = now();
#pragma omp parallel for shared(m1, n1, m2, n2, x, y, R, tmp_threads_time, tmp_threads_load) firstprivate(ind, scores) private(tmp_thread_time, tmp_thread_clock, end) schedule(static)
    for (int i = 1; i < m1; i++)
    {
        tmp_thread_clock = now();
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    R[i][ind.get_indice_cube(j, k, l)] = classic_seq_edit_distance(x[i].substr(1, j), y[k].substr(1, l), scores);
                    // printf("%d -> (%d, %d, %d, %d) = %d\n", omp_get_thread_num(), i, j, k, l, R[i][ind.get_indice_cube(j, k, l)]);
                }
            }
        }
        end = now();
        tmp_thread_time = 0;
        tmp_thread_time = higher_precision_duration(tmp_thread_clock, end);
        tmp_threads_time[omp_get_thread_num()] += tmp_thread_time;
        tmp_threads_load[omp_get_thread_num()]++;
    }
    for (int i = 0; i < (int)tmp_threads_time.size(); i++)
    {
        cout << i << " tid -> " << tmp_threads_time[i] << "s and " << tmp_threads_load[i] << " op" << endl;
    }

    /*for (int c = 0; c < (m1 - 1) * n1 * (m2 - 1) * n2; c++)
    {
        const int i = (c / (n1 * (m2 - 1) * n2)) + 1;
        const int j = (c / ((m2 - 1) * n2)) % n1 + 1;
        const int k = (c / (n2)) % (m2 - 1) + 1;
        const int l = (c) % n2 + 1;
        // printf("%d -> (%d, %d, %d, %d)\n", omp_get_thread_num(), i, j, k, l);
        R[i][ind.get_indice_cube(j, k, l)] = classic_seq_edit_distance(x[i].substr(1, j), y[k].substr(1, l), scores);
    }*/
    end = now();
    tmp_time = higher_precision_duration(tmp_clock, end);
    cout << "[R] = " << tmp_time << "s" << endl;

    /*for (int i = 0; i < m1; i++)
    {
        for (int j = 0; j <= n1; j++)
        {
            for (int k = 0; k < m2; k++)
            {
                for (int l = 0; l <= n2; l++)
                {
                    cout << R[i][j][k][l] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }*/

    tmp_clock = now();
#pragma omp parallel for shared(m1, n1, m2, n2, trans_x, trans_y, C) firstprivate(ind, scores) schedule(static)
    for (int i = 1; i < m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    C[i][ind.get_indice_cube(j, k, l)] = classic_seq_edit_distance(trans_x[j].substr(1, i), trans_y[l].substr(1, k), scores);
                }
            }
        }
    }
    /*for (int c = 0; c < (m1 - 1) * n1 * (m2 - 1) * n2; c++)
    {
        const int i = (c / (n1 * (m2 - 1) * n2)) + 1;
        const int j = (c / ((m2 - 1) * n2)) % n1 + 1;
        const int k = (c / (n2)) % (m2 - 1) + 1;
        const int l = (c) % n2 + 1;
        // printf("(%d, %d, %d, %d)\n", i, j, k, l);
        C[i][ind.get_indice_cube(j, k, l)] = classic_seq_edit_distance(trans_x[j].substr(1, i), trans_y[l].substr(1, k), scores);
    }*/
    end = now();
    tmp_time = higher_precision_duration(tmp_clock, end);
    cout << "[C] = " << tmp_time << "s" << endl;
    /*for (int i = 1; i < n1; i++)
    {
        for (int j = 1; j <= m1; j++)
        {
            for (int k = 1; k < n2; k++)
            {
                for (int l = 1; l <= m2; l++)
                {
                    cout << C[i][j][k][l] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }*/
    end = now();
    preprocessing_time += higher_precision_duration(start, end);

    start = now();
    T[0][0] = 0;
    for (int i = 0; i < m1; i++)
    {
        for (int j = 0; j <= n1; j++)
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
        for (int j = 1; j <= n1; j++)
        {
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
    }
    end = now();
    computation_time += higher_precision_duration(start, end);

    return T[m1 - 1][ind.get_indice_cube(n1, m2 - 1, n2)];
}
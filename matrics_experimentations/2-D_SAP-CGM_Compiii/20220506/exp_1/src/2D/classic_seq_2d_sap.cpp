#include "classic_seq_2d_sap.h"
#include "../common/editing_op.h"
#include "../common/handy_chrono.h"
#include "../1D/classic_seq_edit_distance.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

int classic_seq_2d_sap(const vector<string> &x, const vector<string> &y, const vector<string> &trans_x, const vector<string> &trans_y, double &preprocessing_time, double &computation_time)
{
    int m1 = x.size();
    int n1 = x[1].size();
    int m2 = y.size();
    int n2 = y[1].size();

    vector<vector<int>> dr(m1 + 1, vector<int>(n1 + 1));
    vector<vector<int>> dc(m1 + 1, vector<int>(n1 + 1));

    vector<vector<int>> ir(m2 + 1, vector<int>(n2 + 1));
    vector<vector<int>> ic(m2 + 1, vector<int>(n2 + 1));

    vector<vector<vector<vector<int>>>> R(m1 + 1, vector<vector<vector<int>>>(n1 + 1, vector<vector<int>>(m2 + 1, vector<int>(n2 + 1))));
    vector<vector<vector<vector<int>>>> C(m1 + 1, vector<vector<vector<int>>>(n1 + 1, vector<vector<int>>(m2 + 1, vector<int>(n2 + 1))));
    // vector<vector<vector<vector<int>>>> C(n1 + 1, vector<vector<vector<int>>>(m1 + 1, vector<vector<int>>(n2 + 1, vector<int>(m2 + 1))));

    vector<vector<vector<vector<int>>>> T(m1 + 1, vector<vector<vector<int>>>(n1 + 1, vector<vector<int>>(m2 + 1, vector<int>(n2 + 1))));

    chrono_clock start, end;
    start = now();
    for (int i = 1; i <= m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            dr[i][j] = j * del();
        }
    }

    for (int i = 1; i <= m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            dc[i][j] = i * del();
        }
    }

    for (int i = 1; i <= m2; i++)
    {
        for (int j = 1; j <= n2; j++)
        {
            ir[i][j] = j * ins();
        }
    }

    for (int i = 1; i <= m2; i++)
    {
        for (int j = 1; j <= n2; j++)
        {
            ic[i][j] = i * ins();
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
                    // cout << x[i].substr(1, j) << "-" << y[k].substr(1, l) << endl;
                    R[i][j][k][l] = classic_seq_edit_distance(x[i].substr(1, j), y[k].substr(1, l));
                }
            }
        }
    }

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

    for (int i = 1; i < m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    // cout << trans_x[i].substr(1, j) << "-" << trans_y[k].substr(1, l) << endl;
                    C[i][j][k][l] = classic_seq_edit_distance(trans_x[j].substr(1, i), trans_y[l].substr(1, k));
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
    T[0][0][0][0] = 0;
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
                        T[0][j][k][l] = (k + 1) * (l + 1);
                        T[i][0][k][l] = (k + 1) * (l + 1);
                    }
                    else if ((i > 0 && j > 0) && (k == 0 || l == 0))
                    {
                        T[i][j][0][l] = (i + 1) * (j + 1);
                        T[i][j][k][0] = (i + 1) * (j + 1);
                    }
                    else if ((i == 0 && k == 0) || (j == 0 && l == 0))
                    {
                        T[i][j][k][l] = 0;
                    }
                    else
                    {
                        T[i][j][k][l] = -1;
                    }
                }
            }
        }
    }
    /*for (int i = 0; i < m1; i++)
    {
        for (int j = 0; j <= n1; j++)
        {
            for (int k = 0; k < m2; k++)
            {
                for (int l = 0; l <= n2; l++)
                {
                    cout << T[i][j][k][l] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }*/

    for (int i = 1; i < m1; i++)
    {
        for (int j = 1; j <= n1; j++)
        {
            for (int k = 1; k < m2; k++)
            {
                for (int l = 1; l <= n2; l++)
                {
                    int a = T[i - 1][j][k][l] + dr[i][j];
                    int b = T[i][j - 1][k][l] + dc[i][j];
                    int c = T[i][j][k - 1][l] + ir[k][l];
                    int d = T[i][j][k][l - 1] + ic[k][l];
                    int e = T[i - 1][j][k - 1][l] + R[i][j][k][l];
                    int f = T[i][j - 1][k][l - 1] + C[i][j][k][l];
                    int g = T[i - 1][j - 1][k - 1][l - 1] + C[i - 1][j][k - 1][l] + R[i][j][k][l];
                    int h = T[i - 1][j - 1][k - 1][l - 1] + C[i][j][k][l] + R[i][j - 1][k][l - 1];
                    T[i][j][k][l] = max({a, b, c, d, e, f, g, h});
                }
            }
        }
    }
    end = now();
    computation_time += higher_precision_duration(start, end);

    return T[m1 - 1][n1][m2 - 1][n2];
}
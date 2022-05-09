#include "classic_seq_edit_distance.h"
#include "../common/editing_op.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

int classic_seq_edit_distance(const string &a, const string &b)
{
    int m = a.size();
    int n = b.size();

    vector<vector<int>> scores(m + 1, vector<int>(n + 1));
    /*int **scores = new int *[m + 1];
    scores[0] = new int[(m + 1) * (n + 1)];
    for (int i = 1; i < m + 1; i++)
    {
        scores[i] = scores[0] + i * (n + 1);
    }*/

    scores[0][0] = 0;
    for (int i = 1; i <= m; ++i)
        scores[i][0] = scores[0][0] + del();

    for (int j = 1; j <= n; ++j)
        scores[0][j] = scores[0][0] + ins();

    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            int x = scores[i - 1][j] + del();               // Deletion
            int y = scores[i][j - 1] + ins();               // Insertion
            int z = scores[i - 1][j - 1] + sub(a[i], b[j]); // Substitution

            scores[i][j] = min({x, y, z});
        }
    }

    // display(scores);
    int score = scores[m][n];
    scores.clear();
    return score;
}

int classic_seq_edit_distance(const string &a, const string &b, int **&scores)
{
    int m = a.size();
    int n = b.size();

    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            int x = scores[i - 1][j] + del();               // Deletion
            int y = scores[i][j - 1] + ins();               // Insertion
            int z = scores[i - 1][j - 1] + sub(a[i], b[j]); // Substitution

            scores[i][j] = min({x, y, z});
        }
    }

    // cout << m << " " << n << " " << scores[m][n] << endl;
    int score = scores[m][n];
    return score;
}

void display(const vector<vector<int>> &m)
{
    for (unsigned i = 0; i < m.size(); ++i)
    {
        for (unsigned j = 0; j < m[i].size(); ++j)
        {
            cout << right << setw(4) << m[i][j];
        }

        cout << endl;
    }
}

void classic_seq_edit_distance(const string &a, const string &b, const block &bk, int **&scores)
{
    bound first_bound = bk.data.core_data.first_bound;
    bound second_bound = bk.data.core_data.second_bound;

    for (int i = (first_bound.begin == 0 ? 1 : first_bound.begin); i <= first_bound.end; ++i)
    {
        for (int j = (second_bound.begin == 0 ? 1 : second_bound.begin); j <= second_bound.end; ++j)
        {
            if (a[i - 1] == b[j - 1])
            {
                scores[i][j] = scores[i - 1][j - 1];
                continue;
            }

            int x = scores[i - 1][j] + 1;     // Deletion
            int y = scores[i][j - 1] + 1;     // Insertion
            int z = scores[i - 1][j - 1] + 1; // Substitution

            scores[i][j] = min({x, y, z});
        }
    }
}

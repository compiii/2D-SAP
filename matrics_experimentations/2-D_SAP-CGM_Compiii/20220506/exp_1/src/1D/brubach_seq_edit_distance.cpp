#include "brubach_seq_edit_distance.h"

#include "../common/handy_chrono.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

void display(const vector<vector<duo>> &m);

int brubach_seq_edit_distance(
    const string &a, const string &b,
    double &table_build_duration, double &efficiency)
{
    int m = a.size();
    int n = b.size();

    row_column initial_values = {
        vector<int>(m, 1),
        vector<int>(n, 1)};

    brubach_lookup_table_t lookup_table;

    int t = 0;
    row_column rs = brubach_seq_edit_distance(
        a, b,
        initial_values,
        lookup_table, t,
        table_build_duration, efficiency);

    /*std::cout << "rank " << 0 << " : " << rs.row.size() << " => ";
    for (auto var : rs.row)
    {
        std::cout << right << setw(2) << var;
    }
    cout << std::endl;*/

    // Reconstitute solution
    int sum = n;

    for (int i = 1; i <= m; ++i)
        sum += rs.row[i - 1];

    return sum;
}

row_column brubach_seq_edit_distance(
    const std::string &a, const std::string &b,
    const row_column &initial_values,
    brubach_lookup_table_t &lookup_table, int &t,
    double &table_build_duration, double &efficiency)
{
    int m = a.size();
    int n = b.size();

    if (t == 0)
        t = log10(max(m, n)) * log2(min(m, n));
    int gm = ceil(static_cast<double>(m) / (t - 1)); // number of row blocks
    int gn = ceil(static_cast<double>(n) / (t - 1)); // number of col blocks

    vector<vector<duo>> scores(m + 1, vector<duo>(n + 1));

    int avoided = 0, total = 0;
    chrono_clock start, end;

    // Initialization

    scores[0][0].xy(0, 0);

    int pos = 0;
    for (int i = 1; i <= m; ++i)
    {
        scores[i][0].y(initial_values.row[pos++]);
        scores[i][0].x(0);
    }

    pos = 0;
    for (int j = 1; j <= n; ++j)
    {
        scores[0][j].x(initial_values.column[pos++]);
        scores[0][j].y(0);
    }

    // Precompute substrings

    vector<string> sa(gm);
    vector<string> sb(gn);

    for (int gi = 0; gi < gm; ++gi)
    {
        int row_low = (t - 1) * gi;
        sa[gi] = a.substr(row_low, t - 1);
    }

    for (int gj = 0; gj < gn; ++gj)
    {
        int col_low = (t - 1) * gj;
        sb[gj] = b.substr(col_low, t - 1);
    }

    // Core Computing

    for (int gi = 0; gi < gm; ++gi)
    {
        // indices
        int row_low = (t - 1) * gi;
        int row_high = min(row_low + (t - 1), m);

        for (int gj = 0; gj < gn; ++gj)
        {
            // indices
            int col_low = (t - 1) * gj;
            int col_high = min(col_low + (t - 1), n);

            // Four Russians routine

            // Collect table entries

            // offset vector
            vlong offsets;
            pos = 0;

            for (int i = row_low + 1; i <= row_high; ++i)
            {
                bit_set(offsets, pos++, scores[i][col_low].y());
            }

            for (int j = col_low + 1; j <= col_high; ++j)
            {
                bit_set(offsets, pos++, scores[row_low][j].x());
            }

            // Prevent useless computation
            total++;
            auto &entry = lookup_table[sa[gi]][sb[gj]];

            if (entry.find(offsets) != entry.cend())
            {
                vlong rs = entry[offsets];
                int pos = 0;

                for (int i = row_low + 1; i <= row_high; ++i)
                {
                    scores[i][col_high].y(bit_get(rs, pos++));
                }

                for (int j = col_low + 1; j <= col_high; ++j)
                {
                    scores[row_high][j].x(bit_get(rs, pos++));
                }

                // Skip computation
                avoided++;
                continue;
            }

            //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // Compute blocks entries
            //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            start = now();

            for (int i = row_low + 1; i <= row_high; ++i)
            {
                for (int j = col_low + 1; j <= col_high; ++j)
                {
                    int left = scores[i][j - 1].y();
                    int top = scores[i - 1][j].x();

                    if (a[i - 1] == b[j - 1])
                    {
                        scores[i][j].xy(-left, -top);
                        continue;
                    }

                    int cell = min({0, left, top}) + 1;

                    scores[i][j].xy(cell - left, cell - top);
                }
            }

            //-<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            // Save to prevent further recomputation

            // offset vector
            vlong rs;
            pos = 0;

            for (int i = row_low + 1; i <= row_high; ++i)
            {
                bit_set(rs, pos++, scores[i][col_high].y());
            }

            for (int j = col_low + 1; j <= col_high; ++j)
            {
                bit_set(rs, pos++, scores[row_high][j].x());
            }

            // Save
            entry[offsets] = rs;

            // Account for duration
            end = now();
            table_build_duration += higher_precision_duration(start, end);
        }
    }

    efficiency = static_cast<double>(avoided) / total;

    // Final result
    row_column rs = {
        vector<int>(a.size()),
        vector<int>(b.size())};

    pos = 0;
    for (int i = 1; i <= m; ++i)
    {
        rs.row[pos++] = scores[i][n].y();
    }

    pos = 0;
    for (int j = 1; j <= n; ++j)
    {
        rs.column[pos++] = scores[m][j].x();
    }

    // display(scores);

    // Return
    return rs;
}

int solution(const vector<vector<duo>> &scores)
{
    int m = scores.size() - 1;
    int n = m > 0 ? scores[0].size() - 1 : 0;

    int sum = 0;

    // first row
    for (int j = 1; j <= n; ++j)
        sum += scores[0][j].x();

    // last column
    for (int i = 1; i <= m; ++i)
        sum += scores[i][n].y();

    return sum;
}

void display(const vector<vector<duo>> &m)
{
    std::cout << std::endl;
    for (unsigned i = 0; i < m.size(); ++i)
    {
        for (unsigned j = 0; j < m[i].size(); ++j)
        {
            cout << right << setw(8) << to_string(m[i][j].x()) + "|" + to_string(m[i][j].y());
        }

        cout << endl;
    }
}
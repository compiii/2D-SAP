#include "common/handy_chrono.h"
#include "2D/classic_seq_2d_sap.h"
#include "2D/classic_mpi_2d_sap.h"
#include "2D/classic_openmp_2d_sap.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <mpi.h>
#include <omp.h>

using namespace std;

int main(int argc, char *argv[])
{
    int max_process, rank, algorithm = 0;

    unsigned int size_str_1 = 0, size_str_2 = 0, nb_str_1 = 0, nb_str_2 = 0;
    bool run_all_algorithms = true;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &max_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // omp_set_num_threads(max_process);
    if (argc == 6)
    {
        algorithm = atoi(argv[1]);
        if (algorithm <= 3)
            run_all_algorithms = false;
        nb_str_1 = atoi(argv[2]);
        size_str_1 = atoi(argv[3]);
        nb_str_2 = atoi(argv[4]);
        size_str_2 = atoi(argv[5]);
    }
    else if (argc == 4)
    {
        algorithm = atoi(argv[1]);
        if (algorithm <= 3)
            run_all_algorithms = false;
        nb_str_1 = atoi(argv[2]);
        size_str_1 = atoi(argv[3]);
        nb_str_2 = static_cast<int>(0.90 * nb_str_1);
        size_str_2 = static_cast<int>(0.75 * size_str_1);
    }
    else
    {
        nb_str_1 = 32;
        size_str_1 = 64;
        nb_str_2 = 16;
        size_str_2 = 32;
    }

    vector<string> str_1(nb_str_1 + 1);
    vector<string> str_2(nb_str_2 + 1);

    if constexpr (1)
    {
        ifstream datasets("input/datasets.txt");
        if (datasets)
        {
            string dataset;
            datasets >> dataset;
            ifstream input(dataset);
            if (input)
            {
                string tmp_str_1;
                string tmp_str_2;
                input >> tmp_str_1;
                input >> tmp_str_2;
                if (nb_str_1 * size_str_1 > tmp_str_1.size() || nb_str_2 * size_str_2 > tmp_str_2.size())
                {
                    if (nb_str_1 * size_str_1 > str_1.size())
                    {
                        std::cout << "[ERROR] the max size of the first string is " << tmp_str_1.size() << " in the dataset " << dataset << std::endl;
                    }
                    if (nb_str_2 * size_str_2 > tmp_str_2.size())
                    {
                        std::cout << "[ERROR] the max size of the second string is " << tmp_str_2.size() << " in the dataset " << dataset << std::endl;
                    }
                    exit(0);
                }
                string s = " ";
                for (unsigned int i = 0; i <= size_str_1; i++)
                {
                    s += " ";
                }
                str_1[0] = s;
                s = " ";
                for (unsigned int i = 0; i <= size_str_2; i++)
                {
                    s += " ";
                }
                str_2[0] = s;
                for (unsigned int i = 1; i <= nb_str_1; i++)
                {
                    str_1[i] = " " + tmp_str_1.substr(i * size_str_1, size_str_1);
                }
                for (unsigned int i = 1; i <= nb_str_2; i++)
                {
                    str_2[i] = " " + tmp_str_2.substr(i * size_str_2, size_str_2);
                }
                input.close();
            }
            datasets.close();
        }
    }

    int m1 = str_1.size();
    int n1 = str_1[1].size();
    int m2 = str_2.size();
    int n2 = str_2[1].size();
    vector<vector<char>> tmp_trans_x(n1 + 1, vector<char>(m1 + 1));
    vector<vector<char>> tmp_trans_y(n2 + 1, vector<char>(m2 + 1));
    vector<string> trans_x(n1 + 1);
    vector<string> trans_y(n2 + 1);

    string s = "";
    for (int i = 0; i < m1; ++i)
        for (int j = 0; j <= n1; ++j)
        {
            tmp_trans_x[j][i] = str_1[i][j];
        }

    for (int i = 0; i <= n1; ++i)
    {
        for (int j = 0; j < m1; ++j)
        {
            s += tmp_trans_x[i][j];
        }
        trans_x[i] = s;
        s = "";
    }

    for (int i = 0; i < m2; ++i)
        for (int j = 0; j <= n2; ++j)
        {
            tmp_trans_y[j][i] = str_2[i][j];
        }

    s = "";
    for (int i = 0; i <= n2; ++i)
    {
        for (int j = 0; j < m2; ++j)
        {
            s += tmp_trans_y[i][j];
        }
        trans_y[i] = s;
        s = "";
    }
    tmp_trans_x.clear();
    tmp_trans_y.clear();
    /*cout << "trans string 1" << endl;
    for (auto s : trans_x)
    {
        cout << s << endl;
    }*/
    /*for (int i = 0; i <= n1; ++i)
    {
        for (int j = 0; j < m1; ++j)
        {
            cout << trans_x[i][j];
        }
        cout << endl;
    }*/

    /*cout << "trans string 2" << endl;
    for (auto s : trans_y)
    {
        cout << s << endl;
    }*/
    /*for (int i = 0; i <= n2; ++i)
    {
        for (int j = 0; j < m2; ++j)
        {
            cout << trans_y[i][j];
        }
        cout << endl;
    }*/

    int score = 0;
    chrono_clock start, end;
    double total_duration = 0;
    double computation_time = 0;
    // double communication_time = 0;
    // double table_build_duration = 0;
    // double pure_duration = 0;
    double idle_time = 0;
    double preprocessing_time = 0;
    double pre_comm_time = 0;
    double main_time = 0;
    double main_comm_time = 0;

    if (rank == 0)
    {
        // Print sizes
        cout << "<< " << nb_str_1 << " x " << size_str_1 << endl;
        cout << "<< " << nb_str_2 << " x " << size_str_2 << endl;
        cout << "<< " << max_process << " processors " << endl;
        cout << "<< " << omp_get_max_threads() << " threads " << endl;

        /*cout << "string 1" << endl;
        for (auto s : str_1)
        {
            cout << s << endl;
        }
        cout << "string 2" << endl;
        for (auto s : str_2)
        {
            cout << s << endl;
        }*/
    }

    cout << fixed << setprecision(3);

    if (rank == 0 && (algorithm == 0 || run_all_algorithms))
    {
        // Classic
        total_duration = 0;
        computation_time = 0;
        idle_time = 0;
        preprocessing_time = 0;
        start = now();
        score = classic_seq_2d_sap(str_1, str_2, trans_x, trans_y, preprocessing_time, computation_time);
        end = now();
        total_duration = duration(start, end);
        idle_time = total_duration - (computation_time + preprocessing_time);
        computation_time += idle_time;

        cout << ">> Classic Seq [" << score << "] "
             << total_duration << "s \t= "
             << "(" << preprocessing_time << "s + "
             << computation_time << "s)"
             << endl;
        if constexpr (1)
        {
            ofstream output;
            output.open("output/results.csv", ios::app);
            if (output)
            {
                output << algorithm << "," << 0 << "," << rank << "," << max_process << "," << nb_str_1 << "," << size_str_1 << "," << nb_str_2 << "," << size_str_2 << "," << total_duration << "," << 0 << "," << preprocessing_time << "," << 0 << "," << computation_time << "," << score << "," << std::endl;
                output.close();
            }
        }
    }
    if ((algorithm == 1 || run_all_algorithms))
    {
        MPI_Barrier(MPI_COMM_WORLD);
        // MPI
        total_duration = 0;
        computation_time = 0;
        idle_time = 0;
        preprocessing_time = 0;
        pre_comm_time = 0;
        main_time = 0;
        main_comm_time = 0;
        start = now();
        score = classic_mpi_2d_sap(str_1, str_2, trans_x, trans_y, rank, max_process, preprocessing_time, pre_comm_time, main_time, main_comm_time);
        end = now();
        total_duration = duration(start, end);
        idle_time = total_duration - (preprocessing_time + main_time);
        main_time += idle_time;
        double pre_comp_time = preprocessing_time - pre_comm_time;
        double main_comp_time = main_time - main_comm_time;

        if (rank == (max_process - 1))
        {
            cout << ">> Classic MPI [" << score << "] "
                 << total_duration << "s \t= "
                 << "(" << pre_comp_time << "s + "
                 << pre_comm_time << "s) + ("
                 << main_comp_time << "s + "
                 << main_comm_time << "s)"
                 << endl;
        }

        if constexpr (1)
        {
            ofstream output;
            output.open("output/results.csv", ios::app);
            if (output)
            {
                output << algorithm << "," << 0 << "," << rank << "," << max_process << "," << nb_str_1 << "," << size_str_1 << "," << nb_str_2 << "," << size_str_2 << "," << total_duration << "," << pre_comm_time << "," << pre_comp_time << "," << main_comm_time << "," << main_comp_time << "," << score << "," << std::endl;
                output.close();
            }
        }
    }
    if (rank == 0 && (algorithm == 2 || run_all_algorithms))
    {
        // OpenMP
        total_duration = 0;
        computation_time = 0;
        idle_time = 0;
        preprocessing_time = 0;
        start = now();
        score = classic_openmp_2d_sap(str_1, str_2, trans_x, trans_y, preprocessing_time, computation_time);
        end = now();
        total_duration = duration(start, end);
        idle_time = total_duration - (computation_time + preprocessing_time);
        computation_time += idle_time;

        cout << ">> Classic OpenMP [" << score << "] "
             << total_duration << "s \t= "
             << "(" << preprocessing_time << "s + "
             << computation_time << "s)"
             << endl;
        /*if constexpr (1)
        {
            ofstream output;
            output.open("output/results.csv", ios::app);
            if (output)
            {
                output << algorithm << "," << 0 << "," << rank << "," << max_process << "," << nb_str_1 << "," << size_str_1 << "," << nb_str_2 << "," << size_str_2 << "," << total_duration << "," << 0 << "," << preprocessing_time << "," << 0 << "," << computation_time << "," << score << "," << std::endl;
                output.close();
            }
        }*/
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

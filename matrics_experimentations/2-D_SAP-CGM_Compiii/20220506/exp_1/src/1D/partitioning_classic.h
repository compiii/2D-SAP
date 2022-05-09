#ifndef PARTITIONING_CLASSIC_H
#define PARTITIONING_CLASSIC_H
#include "datatype.h"

#include <vector>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace std;

namespace partitioning_classic
{
    class partitioning
    {
    private:
        MPI_Datatype mpi_vector;
        MPI_Datatype coord_mpi_dt, dim_mpi_dt, bound_mpi_dt, block_core_data_mpi_dt, block_data_mpi_dt;
        long size_first_string, size_second_string;
        int max_frag, max_diag, max_block, max_block_diag, max_eval_block, rank, max_process;
        vector<int> tab_block_data_coord;
        vector<block_data> tab_block_data;
        vector<block> tab_block;

        void construct_all_block_data(const int &rank)
        {
            int id = rank;
            const_block_data_attr const_bd_attr;
            const_bd_attr.init(max_frag, max_block_diag);
            block_data bd;
            while (id < max_block)
            {
                if (max_process == 3 && (rank == 2 || rank == 0))
                {
                    const_block_data_attr const_bd_attr;
                    const_bd_attr.init(max_frag, max_block_diag);
                    bd = construct_block_data_by_id(id, &const_bd_attr);
                }
                else
                {
                    bd = construct_block_data_by_id(id, &const_bd_attr);
                }
                tab_block_data[id] = bd;
                tab_block_data_coord[get_id_block_data_by_coord(bd.core_data.c)] = bd.core_data.id;
                id += max_process;
            }

            create_block_data_as_mpi_datatype();

            bcast_block_data(max_eval_block, 1, max_process);

            update_tab_block_data_coord(rank);
            update_bound_block(rank);

            bcast_block_data(max_eval_block, 1, max_process);
        }

        block_data construct_block_data_by_id(const int &id, const_block_data_attr *const_bd_attr)
        {
            int step = const_bd_attr->step;
            int diag = const_bd_attr->diag;
            int nb_block = const_bd_attr->nb_block;
            int tmp = const_bd_attr->tmp;
            int frag_level = const_bd_attr->frag_level;
            int iter = const_bd_attr->iter;
            int init_tmp = const_bd_attr->init_tmp;
            int eval_lower = const_bd_attr->eval_lower;
            int is_ok = 0;
            int save = 0;
            block_data bd;

            while (iter <= step && is_ok == 0)
            {
                diag++;
                if ((nb_block + tmp) > id)
                {
                    coord c;
                    c.i = id - nb_block + 1;
                    c.j = diag + c.i;
                    block_core_data core_data;
                    core_data.id = id;
                    core_data.c = c;
                    core_data.diag = diag;
                    core_data.address = 0;
                    bd.core_data = core_data;
                    bd.frag_level = frag_level;
                    bd.rank = rank;
                    dim d;
                    d.row = get_nb_row_block(c, frag_level);
                    d.column = get_nb_column_block(c, frag_level);
                    bd.core_data.d = d;
                    is_ok = 1;
                    save = 1;
                }
                if (is_ok == 0 || save == 1)
                {
                    nb_block += tmp;
                    if (eval_lower == 1)
                    {
                        if (iter == (step - 1))
                            tmp = (int)ceil(max_block_diag / 2.0);
                        else
                            tmp++;
                        if (iter == step)
                        {
                            frag_level--;
                            tmp = (int)ceil(max_block_diag / 2.0) + 1;
                            if (frag_level != 0)
                                step += tmp;
                            else
                                step += 2 * floor(max_block_diag / 2.0) - 1;
                        }
                    }
                    else
                    {
                        if (init_tmp == 1)
                        {
                            tmp = max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0);
                            init_tmp = 0;
                        }
                        else
                            tmp--;
                        if (iter == step)
                        {
                            frag_level++;
                            tmp = (int)ceil(max_block_diag / 2.0);
                            init_tmp = 1;
                            if (frag_level == max_frag)
                            {
                                step += max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0);
                            }
                            else
                            {
                                step += (max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1;
                            }
                        }
                    }
                    iter++;
                    if (eval_lower == 1 && iter == (max_block_diag + max_frag * (ceil(max_block_diag / 2.0) + 1)))
                        eval_lower = 0;
                    if (save)
                    {
                        const_bd_attr->save(step, diag, nb_block, tmp, frag_level, iter, init_tmp, eval_lower);
                    }
                }
            }

            return bd;
        }

        int get_size_block_order_to_cyclic(const int &n, const int &p, coord &c, const int &type)
        {
            int diag = c.get_diag(), id;
            int size = get_max_size_block(n, p, 0);
            if ((n + 1) % p == 0)
                return size;
            if (type == 1)
            { //row
                if (diag >= max_block_diag)
                    id = c.i;
                else
                {
                    id = c.i + (max_block_diag - diag);
                }
                if (((size * p + id) > (n + 1)))
                    return size;
                else
                    return size + 1;
            }
            else if (type == 2)
            { //column
                if (diag <= max_block_diag)
                    id = c.i;
                else
                {
                    id = c.i + (diag - max_block_diag);
                }
                if (((size * p + (max_block_diag - id + 1)) > (n + 1)))
                    return size;
                else
                    return size + 1;
            }
            else
            {
                return -1;
            }
        }

        int get_nb_row_block(coord &_c, int &frag_level)
        {
            if (frag_level == 0)
            {
                coord c;
                c.i = _c.i;
                c.j = _c.j - max_frag * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
                return get_size_block_order_to_cyclic(size_first_string, max_process, c, 1);
            }
            else
            {
                int diag = _c.get_diag();
                if (diag < (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                {
                    int limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                    if (diag < limit)
                    {
                        _c.i = _c.i + (limit - diag);
                        _c.j = _c.j + 2 * (limit - diag);
                    }
                    if (diag == (limit + 1))
                    {
                        _c.i = _c.i + 1;
                        _c.j = _c.j + 2;
                        return (int)floor(get_nb_row_block(_c, --frag_level) / 2.0);
                    }
                    else
                    {
                        int a, b;
                        if (_c.i % 2 != 0)
                        {
                            a = (_c.i + 3) / 2.0;
                            b = (2 * _c.j - _c.i + 7) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)floor(get_nb_row_block(_c, --frag_level) / 2.0);
                        }
                        else
                        {
                            a = (_c.i + 2) / 2.0;
                            b = (2 * _c.j - _c.i + 6) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)ceil(get_nb_row_block(_c, --frag_level) / 2.0);
                        }
                    }
                }
                else
                {
                    int limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                    limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                    if (diag > limit)
                        _c.j = _c.i + limit;
                    if (diag == (limit - 1))
                    {
                        _c.j = _c.j - 1;
                        return (int)ceil(get_nb_row_block(_c, --frag_level) / 2.0);
                    }
                    else
                    {
                        int a, b;
                        if (_c.i % 2 != 0)
                        {
                            a = (_c.i + 1) / 2.0;
                            b = (2 * _c.j - _c.i - 3) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)floor(get_nb_row_block(_c, --frag_level) / 2.0);
                        }
                        else
                        {
                            a = _c.i / 2.0;
                            b = (2 * _c.j - _c.i - 4) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)ceil(get_nb_row_block(_c, --frag_level) / 2.0);
                        }
                    }
                }
            }
        }

        int get_nb_column_block(coord &_c, int &frag_level)
        {
            if (frag_level == 0)
            {
                coord c;
                c.i = _c.i;
                c.j = _c.j - max_frag * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
                return get_size_block_order_to_cyclic(size_second_string, max_process, c, 2);
            }
            else
            {
                int diag = _c.get_diag();
                if (diag < (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                {
                    int limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                    if (diag < limit)
                    {
                        _c.j = _c.i + limit;
                    }
                    if (diag == (limit + 1))
                    {
                        _c.j = _c.j - 1;
                        return (int)ceil(get_nb_column_block(_c, --frag_level) / 2.0);
                    }
                    else
                    {
                        int a, b;
                        if (_c.i % 2 != 0)
                        {
                            a = (_c.i + 1) / 2.0;
                            b = (2 * _c.j - _c.i + 5) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)floor(get_nb_column_block(_c, --frag_level) / 2.0);
                        }
                        else
                        {
                            a = _c.i / 2.0;
                            b = (2 * _c.j - _c.i + 4) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)ceil(get_nb_column_block(_c, --frag_level) / 2.0);
                        }
                    }
                }
                else
                {
                    int limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                    limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                    if (diag > limit)
                        _c.i = _c.j - limit;
                    if (diag == (limit - 1))
                    {
                        _c.i = _c.i + 1;
                        return (int)floor(get_nb_column_block(_c, --frag_level) / 2.0);
                    }
                    else
                    {
                        int a, b;
                        if (_c.i % 2 != 0)
                        {
                            a = (_c.i + 3) / 2.0;
                            b = (2 * _c.j - _c.i - 1) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)floor(get_nb_column_block(_c, --frag_level) / 2.0);
                        }
                        else
                        {
                            a = (_c.i + 2) / 2.0;
                            b = (2 * _c.j - _c.i - 2) / 2.0;
                            _c.i = a;
                            _c.j = b;
                            return (int)ceil(get_nb_column_block(_c, --frag_level) / 2.0);
                        }
                    }
                }
            }
        }

        void create_block_data_as_mpi_datatype()
        {
            MPI_Type_contiguous(2, MPI_INT, &coord_mpi_dt);
            MPI_Type_commit(&coord_mpi_dt);

            MPI_Type_contiguous(2, MPI_INT, &dim_mpi_dt);
            MPI_Type_commit(&dim_mpi_dt);

            MPI_Type_contiguous(2, MPI_INT, &bound_mpi_dt);
            MPI_Type_commit(&bound_mpi_dt);

            int count = 7;
            int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
            MPI_Datatype types[7] = {MPI_INT, coord_mpi_dt, MPI_INT, dim_mpi_dt, MPI_INT, coord_mpi_dt, coord_mpi_dt};
            MPI_Aint offsets[7];
            offsets[0] = offsetof(block_core_data, id);
            offsets[1] = offsetof(block_core_data, c);
            offsets[2] = offsetof(block_core_data, address);
            offsets[3] = offsetof(block_core_data, d);
            offsets[4] = offsetof(block_core_data, diag);
            offsets[5] = offsetof(block_core_data, first_bound);
            offsets[6] = offsetof(block_core_data, second_bound);
            MPI_Type_create_struct(count, blocklengths, offsets, types, &block_core_data_mpi_dt);
            MPI_Type_commit(&block_core_data_mpi_dt);

            count = 6;
            int blocklengths2[6] = {1, 1, 1, 1, 1, 1};
            MPI_Datatype types2[6] = {block_core_data_mpi_dt, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
            MPI_Aint offsets2[6];
            offsets2[0] = offsetof(block_data, core_data);
            offsets2[1] = offsetof(block_data, frag_level);
            offsets2[2] = offsetof(block_data, rank);
            offsets2[3] = offsetof(block_data, is_row);
            offsets2[4] = offsetof(block_data, is_column);
            offsets2[5] = offsetof(block_data, is_oblique);
            MPI_Type_create_struct(count, blocklengths2, offsets2, types2, &block_data_mpi_dt);
            MPI_Type_commit(&block_data_mpi_dt);
        }

        void bcast_block_data(const int &c, const int &b, const int &s)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            int id;
            int count = c;
            int blockLengths = b;
            int stride = s;
            MPI_Request r;
            MPI_Type_vector(count, blockLengths, stride, block_data_mpi_dt, &mpi_vector);
            MPI_Type_commit(&mpi_vector);
            for (id = 0; id < max_process; id++)
            {
                if (rank != id)
                {
                    MPI_Isend(&tab_block_data[rank], 1, mpi_vector, id, 0, MPI_COMM_WORLD, &r);
                }
            }
            for (id = 0; id < max_process; id++)
            {
                if (rank != id)
                {
                    count = get_count(id, max_block, max_process);
                    MPI_Type_vector(count, blockLengths, stride, block_data_mpi_dt, &mpi_vector);
                    MPI_Type_commit(&mpi_vector);
                    MPI_Recv(&tab_block_data[id], 1, mpi_vector, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Type_free(&mpi_vector);
        }

        void update_tab_block_data_coord(const int &rank)
        {
            int id = rank, iter = 0, i = 0;
            while (id < max_block)
            {
                for (iter = i; iter < id; iter++)
                {
                    tab_block_data_coord[get_id_block_data_by_coord(tab_block_data[iter].core_data.c)] = tab_block_data[iter].core_data.id;
                }
                i = id + 1;
                id += max_process;
            }
            while (i < max_block)
            {
                tab_block_data_coord[get_id_block_data_by_coord(tab_block_data[i].core_data.c)] = tab_block_data[i].core_data.id;
                i++;
            }
        }

        void update_bound_block(const int &rank)
        {
            int id = rank, offset;
            while (id < max_block)
            {
                offset = get_offset_first_bound(tab_block_data[id]);
                tab_block_data[id].core_data.first_bound.begin = 0 + offset;
                tab_block_data[id].core_data.first_bound.end = tab_block_data[id].core_data.d.row - 1 + offset;

                offset = get_offset_second_bound(tab_block_data[id]);
                tab_block_data[id].core_data.second_bound.begin = 0 + offset;
                tab_block_data[id].core_data.second_bound.end = tab_block_data[id].core_data.d.column - 1 + offset;

                id += max_process;
            }
        }

        block_data get_block_data_by_coord(const coord &c)
        {
            return tab_block_data[tab_block_data_coord[get_id_block_data_by_coord(c)]];
        }

        int get_offset_first_bound(block_data &root)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int column, frag_level = root.frag_level, diag = _c.get_diag(), limit = -1, exit;
            column = 1;

            c.i = _c.i;
            c.j = _c.j - (max_frag - frag_level) * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
            if (c.i == c.get_diag())
                column = 0;

            if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    ;
                else
                {
                    limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                }
            }
            else
            {
                limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
            }

            if (column)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag <= max_block_diag)
                        c.j--;
                    else
                        c.i++;
                    bd = get_block_data_by_coord(c);
                    return bd.core_data.d.row + get_offset_first_bound(bd);
                }
                else if (column)
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                            if (diag == (limit + 2))
                            {
                                c = _c;
                                c.j--;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.row + get_offset_first_bound(bd);
                            }
                            else
                            {
                                c = _c;
                                if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                    c.j--;
                                else
                                    c.i++;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.row + get_offset_first_bound(bd);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                c.i = 2 * _c.i;
                                c.j = _c.j + _c.i - 1;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.row + get_offset_first_bound(bd);
                            }
                            if (diag <= limit && max_frag != frag_level && exit == 0)
                            {
                                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                                if (diag == (limit + 2))
                                {
                                    exit = 1;
                                    c = _c;
                                    c.j--;
                                    bd = get_block_data_by_coord(c);
                                    return bd.core_data.d.row + get_offset_first_bound(bd);
                                }
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.j--;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.row + get_offset_first_bound(bd);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        if (diag == (limit - 1) || diag > limit)
                        {
                            c.i++;
                        }
                        else
                        {
                            int a, b;
                            if (c.i % 2 != 0)
                            {
                                a = (_c.i + 1) / 2.0;
                                b = (2 * _c.j - _c.i - 1) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                            else
                            {
                                a = (_c.i + 2) / 2.0;
                                b = (2 * _c.j - _c.i - 2) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                        }
                        bd = get_block_data_by_coord(c);
                        return bd.core_data.d.row + get_offset_first_bound(bd);
                    }
                }
            }
            return 0;
        }

        int get_offset_second_bound(block_data &root)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int row, frag_level = root.frag_level, diag = _c.get_diag(), limit, exit;
            row = 1;

            if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    row = 0;
                else
                {
                    limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                    if (diag != (limit + 1) && _c.i == 1)
                        row = 0;
                }
            }
            else
            {
                limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
            }

            if (row)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag <= max_block_diag)
                    {
                        c.i--;
                        c.j -= 2;
                    }
                    else
                        c.j--;
                    bd = get_block_data_by_coord(c);
                    return bd.core_data.d.column + get_offset_second_bound(bd);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                            if (diag == (limit + 2))
                            {
                                c = _c;
                                c.i = _c.i - 1;
                                c.j = _c.j - 2;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.column + get_offset_second_bound(bd);
                            }
                            else
                            {
                                c = _c;
                                if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.i--;
                                    c.j -= 2;
                                }
                                else
                                {
                                    c.j--;
                                }
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.column + get_offset_second_bound(bd);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                c.i = 2 * _c.i - 1;
                                c.j = (2 * _c.j + c.i - 3) / 2.0;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.column + get_offset_second_bound(bd);
                            }
                            if (diag <= limit && max_frag != frag_level && exit == 0)
                            {
                                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                                if (diag == (limit + 2))
                                {
                                    exit = 1;
                                    c = _c;
                                    c.i = _c.i - 1;
                                    c.j = _c.j - 2;
                                    bd = get_block_data_by_coord(c);
                                    return bd.core_data.d.column + get_offset_second_bound(bd);
                                }
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.i--;
                                c.j -= 2;
                                bd = get_block_data_by_coord(c);
                                return bd.core_data.d.column + get_offset_second_bound(bd);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        if (diag == (limit - 1) || diag > limit)
                        {
                            c.j--;
                        }
                        else
                        {
                            int a, b;
                            if (_c.i % 2 != 0)
                            {
                                a = (_c.i + 1) / 2.0;
                                b = (2 * _c.j - _c.i - 3) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                            else
                            {
                                a = _c.i / 2.0;
                                b = 1 + ((2 * _c.j - _c.i - 4) / 2.0);
                                c.i = a;
                                c.j = b;
                            }
                        }
                        bd = get_block_data_by_coord(c);
                        return bd.core_data.d.column + get_offset_second_bound(bd);
                    }
                }
            }
            return 0;
        }

        void construct_all_block(const int &rank)
        {
            int id = rank;
            int i = 0;
            vector<block_data> all_need_block_data;
            while (id < max_block)
            {
                tab_block[i].data = tab_block_data[id];
                tab_block[i].need = get_need_block_data_list(tab_block_data[id], &all_need_block_data);

                if (!tab_block[i].need.empty())
                {
                    vector<block_data> copie((tab_block[i].need));
                    all_need_block_data.insert(
                        all_need_block_data.end(),
                        make_move_iterator(copie.begin()),
                        make_move_iterator(copie.end()));
                }
                tab_block[i].depend = get_depend_block_data_list(tab_block_data[id]);
                i++;
                id += max_process;
            }
        }

        vector<block_data> get_need_block_data_list(block_data &root, vector<block_data> *all_need_block_data)
        {
            vector<block_data> need_block_data_list;
            get_need_block_data_list_row(root, need_block_data_list, *all_need_block_data);
            get_need_block_data_list_column(root, need_block_data_list, *all_need_block_data);
            get_need_block_data_list_oblique(root, need_block_data_list, *all_need_block_data);
            bubble_sort_block_data_list(need_block_data_list);
            return need_block_data_list;
        }

        void bubble_sort_block_data_list(vector<block_data> &list)
        {
            block_data tmp;
            if (list.empty())
                return;

            for (unsigned int i = 0; i < list.size(); i++)
            {
                for (unsigned int j = i; j < list.size(); j++)
                {
                    if (list[i].core_data.id > list[j].core_data.id)
                    {
                        tmp = list[i];
                        list[i] = list[j];
                        list[j] = tmp;
                    }
                }
            }
        }

        void add_block_data_to_block_data_list_filter(vector<block_data> &block_data_list, block_data bd, const int &option, const int &level, vector<block_data> &all_need_block_data)
        {
            if (!find_block_data_to_block_data_list_and_update_level(all_need_block_data, bd, level))
            {
                add_block_data_to_block_data_list(block_data_list, bd, option, level);
            }
        }

        void add_block_data_to_block_data_list(vector<block_data> &block_data_list, block_data &bd, const int &option, const int &level)
        {
            if (block_data_list.empty())
            {
                if (level == 0)
                {
                    bd.is_row = 1;
                }
                else if (level == 1)
                {
                    bd.is_column = 1;
                }
                else
                {
                    bd.is_oblique = 1;
                }
                block_data_list.push_back(bd);
            }
            else
            {
                int is_present = 0;
                if (option == 0)
                    for (auto &path : block_data_list)
                    {
                        if (path.rank == bd.rank)
                        {
                            is_present = 1;
                            if (path.core_data.diag > bd.core_data.diag)
                            {
                                path = bd;
                            }
                            if (level == 0)
                            {
                                path.is_row = 1;
                            }
                            else if (level == 1)
                            {
                                path.is_column = 1;
                            }
                            else
                            {
                                path.is_oblique = 1;
                            }
                            break;
                        }
                    }
                if (is_present == 0)
                {
                    if (level == 0)
                    {
                        bd.is_row = 1;
                    }
                    else if (level == 1)
                    {
                        bd.is_column = 1;
                    }
                    else
                    {
                        bd.is_oblique = 1;
                    }
                    block_data_list.push_back(bd);
                }
            }
            return;
        }

        bool find_block_data_to_block_data_list(const vector<block_data> &block_data_list, const block_data &bd)
        {
            if (block_data_list.empty())
                return false;
            for (auto path : block_data_list)
            {
                if (path.core_data.id == bd.core_data.id)
                    return true;
            }
            return false;
        }

        bool find_block_data_to_block_data_list_and_update_level(const vector<block_data> &block_data_list, const block_data &bd, const int &level)
        {
            if (block_data_list.empty())
                return false;
            for (auto path : block_data_list)
            {
                if (path.core_data.id == bd.core_data.id)
                {
                    int is_update = 0;
                    for (int i = 0; (i < max_eval_block && is_update == 0); i++)
                    {
                        if (find_block_data_to_block_data_list(tab_block[i].need, bd))
                        {
                            is_update = 1;
                            for (auto &second_path : tab_block[i].need)
                            {
                                if (second_path.core_data.id == bd.core_data.id)
                                {
                                    if (level == 0)
                                    {
                                        second_path.is_row = 1;
                                    }
                                    else if (level == 1)
                                    {
                                        second_path.is_column = 1;
                                    }
                                    else
                                    {
                                        second_path.is_oblique = 1;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    return true;
                }
            }

            return false;
        }

        void get_need_block_data_list_row(block_data &root, vector<block_data> &block_data_list, vector<block_data> &all_need_block_data)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int row, frag_level = root.frag_level, diag = _c.get_diag(), limit, exit;
            row = 1;

            if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    row = 0;
                else
                {
                    limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                    if (diag != (limit + 1) && _c.i == 1)
                        row = 0;
                }
            }
            else
            {
                limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
            }

            if (row)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag <= max_block_diag)
                    {
                        c.i--;
                        c.j -= 2;
                    }
                    else
                        c.j--;
                    bd = get_block_data_by_coord(c);
                    if (root.rank != bd.rank)
                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                            if (diag == (limit + 2))
                            {
                                c = _c;
                                c.i = _c.i - 1;
                                c.j = _c.j - 2;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                                c.i = 2 * _c.i - 2;
                                c.j = _c.i + _c.j - 4;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                            }
                            else
                            {
                                c = _c;
                                if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.i--;
                                    c.j -= 2;
                                }
                                else
                                {
                                    c.j--;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                c.i = 2 * _c.i - 1;
                                c.j = (2 * _c.j + c.i - 3) / 2.0;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                            }
                            if (diag <= limit && max_frag != frag_level && exit == 0)
                            {
                                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                                if (diag == (limit + 2))
                                {
                                    exit = 1;
                                    c = _c;
                                    c.i = _c.i - 1;
                                    c.j = _c.j - 2;
                                    bd = get_block_data_by_coord(c);
                                    if (root.rank != bd.rank)
                                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                                    c.i = 2 * _c.i - 2;
                                    c.j = _c.i + _c.j - 4;
                                    bd = get_block_data_by_coord(c);
                                    if (root.rank != bd.rank)
                                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                                }
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.i--;
                                c.j -= 2;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        if (diag == (limit - 1) || diag > limit)
                        {
                            c.j--;
                        }
                        else
                        {
                            int a, b;
                            if (_c.i % 2 != 0)
                            {
                                a = (_c.i + 1) / 2.0;
                                b = (2 * _c.j - _c.i - 3) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                            else
                            {
                                a = _c.i / 2.0;
                                b = 1 + ((2 * _c.j - _c.i - 4) / 2.0);
                                c.i = a;
                                c.j = b;
                            }
                        }
                        bd = get_block_data_by_coord(c);
                        if (root.rank != bd.rank)
                            add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 0, all_need_block_data);
                    }
                }
            }
        }

        void get_need_block_data_list_column(block_data &root, vector<block_data> &block_data_list, vector<block_data> &all_need_block_data)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int column, frag_level = root.frag_level, diag = _c.get_diag(), limit = 0, exit;
            column = 1;

            c.i = _c.i;
            c.j = _c.j - (max_frag - frag_level) * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
            if (c.i == c.get_diag())
                column = 0;

            if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    ;
                else
                {
                    limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                }
            }
            else
            {
                limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
            }

            if (column)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag <= max_block_diag)
                        c.j--;
                    else
                        c.i++;
                    bd = get_block_data_by_coord(c);
                    if (root.rank != bd.rank)
                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                            if (diag == (limit + 2))
                            {
                                c = _c;
                                c.j--;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                                c.i = 2 * _c.i - 1;
                                c.j = _c.i + _c.j - 3;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                            }
                            else
                            {
                                c = _c;
                                if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                    c.j--;
                                else
                                    c.i++;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                c.i = 2 * _c.i;
                                c.j = _c.j + _c.i - 1;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                            }
                            if (diag <= limit && max_frag != frag_level && exit == 0)
                            {
                                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                                if (diag == (limit + 2))
                                {
                                    exit = 1;
                                    c = _c;
                                    c.j--;
                                    bd = get_block_data_by_coord(c);
                                    if (root.rank != bd.rank)
                                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                                    c.i = 2 * _c.i - 1;
                                    c.j = _c.i + _c.j - 3;
                                    bd = get_block_data_by_coord(c);
                                    if (root.rank != bd.rank)
                                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                                }
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.j--;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        if (diag == (limit - 1) || diag > limit)
                        {
                            c.i++;
                        }
                        else
                        {
                            int a, b;
                            if (c.i % 2 != 0)
                            {
                                a = (_c.i + 1) / 2.0;
                                b = (2 * _c.j - _c.i - 1) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                            else
                            {
                                a = (_c.i + 2) / 2.0;
                                b = (2 * _c.j - _c.i - 2) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                        }
                        bd = get_block_data_by_coord(c);
                        if (root.rank != bd.rank)
                            add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 1, all_need_block_data);
                    }
                }
            }
        }

        void get_need_block_data_list_oblique(block_data &root, vector<block_data> &block_data_list, vector<block_data> &all_need_block_data)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int oblique, frag_level = root.frag_level, diag = _c.get_diag(), limit, exit;
            oblique = 1;

            c.i = _c.i;
            c.j = _c.j - (max_frag - frag_level) * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
            if (c.i == c.get_diag())
                oblique = 0;

            if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    oblique = 0;
                else
                {
                    limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                    if (diag != (limit + 1) && _c.i == 1)
                        oblique = 0;
                }
            }
            else
            {
                limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                if (diag == limit)
                    oblique = 0;
            }

            if (oblique)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag == (max_block_diag + 1))
                    {
                        c.j -= 2;
                    }
                    else if (diag > max_block_diag)
                    {
                        c.i++;
                        c.j--;
                    }
                    else
                    {
                        c.i--;
                        c.j -= 3;
                    }
                    bd = get_block_data_by_coord(c);
                    if (root.rank != bd.rank)
                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                            if (diag == (limit + 2))
                            {
                                c.i = 2 * _c.i - 2;
                                c.j = _c.i + _c.j - 5;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                            }
                            else if (diag == (limit + 3))
                            {
                                c.i--;
                                c.j -= 3;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                            }
                            else
                            {
                                c = _c;
                                if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.i--;
                                    c.j -= 3;
                                }
                                else if ((diag - 1) == (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.j -= 2;
                                }
                                else
                                {
                                    c.i++;
                                    c.j--;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                            }
                        }
                        else
                        {
                            exit = 0;
                            c = _c;
                            if (diag == (limit + 1))
                            {
                                c.i = 2 * _c.i - 1;
                                c.j = _c.j + _c.i - 3;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                            }
                            if (diag <= limit && max_frag != frag_level && exit == 0)
                            {
                                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                                if (diag == (limit + 2))
                                {
                                    exit = 1;
                                    c.i = 2 * _c.i - 2;
                                    c.j = _c.i + _c.j - 5;
                                    bd = get_block_data_by_coord(c);
                                    if (root.rank != bd.rank)
                                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                                }
                                else if (diag == (limit + 3))
                                {
                                    exit = 1;
                                    c.i--;
                                    c.j -= 3;
                                    bd = get_block_data_by_coord(c);
                                    if (root.rank != bd.rank)
                                        add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                                }
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.i--;
                                c.j -= 3;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        if (diag == (limit - 1))
                        {
                            if ((max_block_diag == 2) || (max_block_diag == 3 && frag_level == 1))
                                c.j -= 2;
                            else
                            {
                                c.i++;
                                c.j--;
                            }
                        }
                        else if (diag == (limit + 1))
                        {
                            int a, b;
                            if (c.i % 2 != 0)
                            {
                                a = (_c.i + 1) / 2.0;
                                b = 1 + ((2 * _c.j - _c.i - 5) / 2.0);
                                c.i = a;
                                c.j = b;
                            }
                            else
                            {
                                a = (_c.i + 2) / 2.0;
                                b = (2 * _c.j - _c.i - 4) / 2.0;
                                c.i = a;
                                c.j = b;
                            }
                        }
                        else
                        {
                            c.i++;
                            c.j--;
                        }
                        bd = get_block_data_by_coord(c);
                        if (root.rank != bd.rank)
                            add_block_data_to_block_data_list_filter(block_data_list, bd, 1, 2, all_need_block_data);
                    }
                }
            }
        }

        vector<block_data> get_depend_block_data_list(block_data &root)
        {
            vector<block_data> dependBlock_data_list;
            get_depend_block_data_list_row(root, dependBlock_data_list);
            get_depend_block_data_list_column(root, dependBlock_data_list);
            get_depend_block_data_list_oblique(root, dependBlock_data_list);

            return dependBlock_data_list;
        }

        void get_depend_block_data_list_row(block_data &root, vector<block_data> &block_data_list)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int row, frag_level = root.frag_level, diag = _c.get_diag(), limit = 0, exit;
            row = 1;

            c.i = _c.i;
            c.j = _c.j + (max_frag - frag_level) * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
            if (c.j == (max_diag + 1))
                row = 0;

            if (diag >= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    ;
                else
                {
                    limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                    limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                }
            }
            else
            {
                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
            }

            if (row)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag < max_block_diag)
                    {
                        c.i++;
                        c.j += 2;
                    }
                    else
                    {
                        c.j++;
                    }
                    bd = get_block_data_by_coord(c);
                    if (root.rank != bd.rank)
                        add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (int)floor(max_block_diag / 2.0) + 2;
                            limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                            if (diag == (limit - 2))
                            {
                                c = _c;
                                c.j++;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                                c.i = 2 * _c.i - 1;
                                c.j = _c.j + _c.i + 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                            }
                            else
                            {
                                c = _c;
                                if (diag < (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.i++;
                                    c.j += 2;
                                }
                                else
                                {
                                    c.j++;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                c.i++;
                                c.j += 2;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                            }
                            if (diag == limit && exit == 0)
                            {
                                int a, b;
                                if (_c.i % 2 != 0)
                                {
                                    a = (_c.i + 1) / 2.0;
                                    b = (2 * _c.j - _c.i + 3) / 2.0;
                                    c.i = a;
                                    c.j = b;
                                }
                                else
                                {
                                    a = (_c.i + 2) / 2.0;
                                    b = (2 * _c.j - _c.i + 6) / 2.0;
                                    c.i = a;
                                    c.j = b;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                                exit = 1;
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.i++;
                                c.j += 2;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        exit = 0;
                        if (diag == (limit - 1))
                        {
                            c.i *= 2;
                            c.j += _c.i + 1;
                            exit = 1;
                            bd = get_block_data_by_coord(c);
                            if (root.rank != bd.rank)
                                add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                        }
                        if (diag >= limit && max_frag != frag_level && exit == 0)
                        {
                            limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                            limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                            if (diag == (limit - 2))
                            {
                                exit = 1;
                                c = _c;
                                c.j++;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                                c.i = 2 * _c.i - 1;
                                c.j = _c.j + _c.i + 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                            }
                        }
                        if (exit == 0)
                        {
                            c = _c;
                            c.j++;
                            bd = get_block_data_by_coord(c);
                            if (root.rank != bd.rank)
                                add_block_data_to_block_data_list(block_data_list, bd, 0, 0);
                        }
                    }
                }
            }
        }

        void get_depend_block_data_list_column(block_data &root, vector<block_data> &block_data_list)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int column, frag_level = root.frag_level, diag = _c.get_diag(), limit, exit;
            column = 1;

            if (diag >= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    column = 0;
                else
                {
                    limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                    limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                    if (diag != (limit - 1) && _c.i == 1)
                        column = 0;
                }
            }
            else
            {
                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
            }

            if (column)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag < max_block_diag)
                    {
                        c.j++;
                    }
                    else
                    {
                        c.i--;
                    }
                    bd = get_block_data_by_coord(c);
                    if (root.rank != bd.rank)
                        add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (int)floor(max_block_diag / 2.0) + 2;
                            limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                            if (diag == (limit - 2))
                            {
                                c = _c;
                                c.i--;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                                c.i = 2 * _c.i - 2;
                                c.j = _c.j + _c.i;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                            }
                            else
                            {
                                c = _c;
                                if (diag < (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.j++;
                                }
                                else
                                {
                                    c.i--;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                c.j++;
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                            }
                            if (diag == limit && exit == 0)
                            {
                                int a, b;
                                if (_c.i % 2 != 0)
                                {
                                    a = (_c.i + 1) / 2.0;
                                    b = (2 * _c.j - _c.i + 5) / 2.0;
                                    c.i = a;
                                    c.j = b;
                                }
                                else
                                {
                                    a = _c.i / 2.0;
                                    b = ((2 * _c.j - _c.i + 4) / 2.0) - 1;
                                    c.i = a;
                                    c.j = b;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                                exit = 1;
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.j++;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        exit = 0;
                        if (diag == (limit - 1))
                        {
                            c.i = 2 * _c.i - 1;
                            c.j = _c.i + _c.j;
                            exit = 1;
                            bd = get_block_data_by_coord(c);
                            if (root.rank != bd.rank)
                                add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                        }
                        if (diag >= limit && max_frag != frag_level && exit == 0)
                        {
                            limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                            limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                            if (diag == (limit - 2))
                            {
                                exit = 1;
                                c = _c;
                                c.i--;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                                c.i = 2 * _c.i - 2;
                                c.j = _c.j + _c.i;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                            }
                        }
                        if (exit == 0)
                        {
                            c = _c;
                            c.i--;
                            bd = get_block_data_by_coord(c);
                            if (root.rank != bd.rank)
                                add_block_data_to_block_data_list(block_data_list, bd, 0, 1);
                        }
                    }
                }
            }
        }

        void get_depend_block_data_list_oblique(block_data &root, vector<block_data> &block_data_list)
        {
            block_data bd;
            coord _c = root.core_data.c, c;
            int oblique, frag_level = root.frag_level, diag = _c.get_diag(), limit, exit;
            oblique = 1;

            c.i = _c.i;
            c.j = _c.j + (max_frag - frag_level) * (max_block_diag - (int)floor(max_block_diag / 2.0) + 1);
            if (c.j == (max_diag + 1))
                oblique = 0;

            if (diag >= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
            {
                if (frag_level == 0 && _c.i == 1)
                    oblique = 0;
                else
                {
                    limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level - 1) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                    limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                    if (diag != (limit - 1) && _c.i == 1)
                        oblique = 0;
                }
            }
            else
            {
                limit = (max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0)) + (max_frag - frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1) - 1;
                if (diag == limit)
                    oblique = 0;
            }

            if (oblique)
            {
                c = _c;
                if (max_frag == 0)
                {
                    if (diag >= max_block_diag)
                    {
                        c.i--;
                        c.j++;
                    }
                    else if (diag == (max_block_diag - 1))
                    {
                        c.j += 2;
                    }
                    else
                    {
                        c.i++;
                        c.j += 3;
                    }
                    bd = get_block_data_by_coord(c);
                    if (root.rank != bd.rank)
                        add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                }
                else
                {
                    if (diag <= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1) + (int)floor(max_block_diag / 2.0) - 1))
                    {
                        if (frag_level == 0)
                        {
                            limit = (int)floor(max_block_diag / 2.0) + 2;
                            limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                            if (diag == (limit - 2))
                            {
                                c.i = 2 * _c.i - 2;
                                c.j = _c.j + _c.i + 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                            else if (diag == (limit - 3))
                            {
                                c.i--;
                                c.j++;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                            else
                            {
                                c = _c;
                                if (diag >= (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)))
                                {
                                    c.i--;
                                    c.j++;
                                }
                                else if (diag == ((max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1))
                                {
                                    c.j += 2;
                                }
                                else
                                {
                                    c.i++;
                                    c.j += 3;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                        }
                        else
                        {
                            exit = 0;
                            if (diag == (limit + 1))
                            {
                                c = _c;
                                if (max_block_diag == 3 || max_block_diag == 2)
                                {
                                    c.j += 2;
                                }
                                else
                                {
                                    c.i++;
                                    c.j += 3;
                                }
                                exit = 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                            if (diag == (limit - 1) && exit == 0)
                            {
                                int a, b;
                                if (_c.i % 2 != 0)
                                {
                                    a = (_c.i + 1) / 2.0;
                                    b = (2 * _c.j - _c.i + 5) / 2.0;
                                    c.i = a;
                                    c.j = b;
                                }
                                else
                                {
                                    a = (_c.i + 2) / 2.0;
                                    b = (2 * _c.j - _c.i + 8) / 2.0;
                                    c.i = a;
                                    c.j = b;
                                }
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                                exit = 1;
                            }
                            if (exit == 0)
                            {
                                c = _c;
                                c.i++;
                                c.j += 3;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                        }
                    }
                    else
                    {
                        c = _c;
                        exit = 0;
                        if (diag == (limit - 1))
                        {
                            c.i = 2 * _c.i - 1;
                            c.j = _c.i + _c.j + 1;
                            exit = 1;
                            bd = get_block_data_by_coord(c);
                            if (root.rank != bd.rank)
                                add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                        }
                        if (diag >= limit && max_frag != frag_level && exit == 0)
                        {
                            limit = (int)floor(max_block_diag / 2.0) + 2 + (frag_level) * ((max_block_diag + (max_block_diag % 2 != 0 ? 1 : 0)) / 2.0 + 1);
                            limit += (max_block_diag + max_frag * ((int)ceil(max_block_diag / 2.0) + 1)) - 1;
                            if (diag == (limit - 2))
                            {
                                exit = 1;
                                c = _c;
                                c.i = 2 * _c.i - 2;
                                c.j = _c.j + _c.i + 1;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                            else if (diag == (limit - 3))
                            {
                                exit = 1;
                                c = _c;
                                c.i--;
                                c.j++;
                                bd = get_block_data_by_coord(c);
                                if (root.rank != bd.rank)
                                    add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                            }
                        }
                        if (exit == 0)
                        {
                            c = _c;
                            c.i--;
                            c.j++;
                            bd = get_block_data_by_coord(c);
                            if (root.rank != bd.rank)
                                add_block_data_to_block_data_list(block_data_list, bd, 0, 2);
                        }
                    }
                }
            }
        }

        inline int get_count(const int &id, const int &max_block, const int &max_process)
        {
            int count = (int)floor(max_block / (float)max_process);
            if (max_block % max_process == 0 || ((count * max_process + id + 1) > max_block))
                return count;
            else
                return count + 1;
        }

        inline int get_max_size_block(const int &n, const int &p, const int &type)
        {
            return (int)(type == 1 ? ceil((n + 1) / (float)p) : floor((n + 1) / (float)p));
        }

        inline int get_id_block_data_by_coord(const coord &c)
        {
            int i = c.i;
            int j = c.j;
            return (max_diag - i + 1) * (i - 1) + i * (i - 1) / 2.0 + (j - i - 1);
        }

    public:
        partitioning() {}
        vector<block> build(const string &a, const string &b, const int &_rank, const int &_max_process, int &_max_diag)
        {
            rank = _rank;
            max_process = _max_process;
            size_first_string = a.size();
            size_second_string = b.size();
            max_frag = 0;
            max_block_diag = max_process;
            _max_diag = 2 * max_block_diag - 1 + 2 * max_frag * (ceil(max_block_diag / 2.0) + 1);
            max_diag = _max_diag;
            max_block = pow(max_block_diag, 2) + max_frag * (max_block_diag + 1) * (max_block_diag + 2 * (max_block_diag % 2)) - max_frag * (ceil(max_block_diag / 2.0) * (ceil(max_block_diag / 2.0) - 1));
            max_eval_block = get_count(rank, max_block, max_process);

            tab_block_data.resize(max_block);
            coord c;
            if (max_frag > 0 && (max_block_diag % 2 != 0))
            {
                c.i = max_block_diag + 1;
            }
            else
            {
                c.i = max_block_diag;
            }
            c.j = max_diag + 1;

            tab_block_data_coord.resize(get_id_block_data_by_coord(c) + 1);
            tab_block.resize(max_eval_block);

            construct_all_block_data(rank);

            construct_all_block(rank);

            return tab_block;
        }
    };

    vector<block> build(const string &a, const string &b, const int &_rank, const int &_max_process, int &_max_diag)
    {
        partitioning p;
        return p.build(a, b, _rank, _max_process, _max_diag);
    }

} // namespace partitioning_classic

#endif // PARTITIONING_H
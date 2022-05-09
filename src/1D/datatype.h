#ifndef DATATYPE_H
#define DATATYPE_H
#include <vector>

class coord
{
public:
    int i;
    int j;

    coord(int i, int j) : i(i), j(j) {}
    coord() : i(0), j(0) {}
    int get_diag() { return j - i; }
};

class dim
{
public:
    int row;
    int column;

    dim() : row(0), column(0){};
    dim(int row, int column) : row(row), column(column){};
    int get_size() { return row * column; }
};

class bound
{
public:
    int begin;
    int end;

    bound() : begin(0), end(0){};
    bound(int begin, int end) : begin(begin), end(end){};
    int length() { return end - begin + 1; }
};

class block_core_data
{
public:
    int id;
    coord c;
    int address;
    dim d;
    int diag;
    bound first_bound;
    bound second_bound;

    block_core_data() {}
};

class block_data
{
public:
    block_core_data core_data;
    int frag_level;
    int rank;
    int is_row;
    int is_column;
    int is_oblique;

    block_data() {}
};

class block
{
public:
    block_data data;
    std::vector<block_data> need;
    std::vector<block_data> depend;

    block() {}
};

class const_block_data_attr
{
public:
    int step;
    int diag;
    int nb_block; // number of blocks before the diagonal block id (it is_ not easy to define)
    int tmp;
    int frag_level;
    int iter; // iterator
    int init_tmp;
    int eval_lower;

    const_block_data_attr() {}
    void init(int max_frag, int max_block_diag)
    {
        step = (max_frag != 0 ? max_block_diag + 1 + (max_block_diag % 2 != 0 ? 1 : 0) : 2 * max_block_diag - 1);
        diag = 0;
        nb_block = 0;
        tmp = 1;
        frag_level = max_frag;
        iter = 1;
        init_tmp = 0;
        eval_lower = 1; // default value is 1 because the evaluation begins to the lower triangular matrix
    }
    void save(int step, int diag, int nb_block, int tmp, int frag_level, int iter, int init_tmp, int eval_lower)
    {
        this->step = step;
        this->diag = diag;
        this->nb_block = nb_block;
        this->tmp = tmp;
        this->frag_level = frag_level;
        this->iter = iter;
        this->init_tmp = init_tmp;
        this->eval_lower = eval_lower;
    }
};
#endif // DATATYPE_H
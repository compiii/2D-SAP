#ifndef INDICE_H
#define INDICE_H

class indice
{
public:
    indice() { indice(0, 0, 0, 0); }

    indice(int _m1, int _n1, int _m2, int _n2) : m1(_m1), n1(_n1), m2(_m2), n2(_n2)
    {
    }

    int get_indice_mat(int k, int l)
    {
        return k * (n2 + 1) + l;
    }

    int get_indice_cube(int j, int k, int l)
    {
        return j * (m2 + 1) * (n2 + 1) + get_indice_mat(k, l);
    }

private:
    int m1;
    int n1;
    int m2;
    int n2;
};

#endif // INDICE_H

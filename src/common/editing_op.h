#ifndef EDITING_OP_H
#define EDITING_OP_H

inline int ins()
{
    return -1;
}

inline int del()
{
    return -1;
}

inline int sub(const char &a, const char &b)
{
    return (a == b) ? 2 : -1;
}
#endif
#ifndef CLASSIC_SEQ_EDIT_DISTANCE_H
#define CLASSIC_SEQ_EDIT_DISTANCE_H

#include "datatype.h"

#include <string>
#include <vector>

void display(const std::vector<std::vector<int>> &m);
int classic_seq_edit_distance(const std::string &a, const std::string &b);
int classic_seq_edit_distance(const std::string &a, const std::string &b, int **&scores);
void classic_seq_edit_distance(const std::string &a, const std::string &b, const block &bk, int **&scores);

#endif // SIMPLE_EDIT_DISTANCE_H

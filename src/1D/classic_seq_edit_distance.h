#ifndef CLASSIC_SEQ_EDIT_DISTANCE_H
#define CLASSIC_SEQ_EDIT_DISTANCE_H

#include <string>
#include <vector>

void display(const std::vector<std::vector<int>> &m);
int classic_seq_edit_distance(const std::string &a, const std::string &b);
int classic_seq_edit_distance(const std::string &a, const std::string &b, int **&scores);

#endif // SIMPLE_EDIT_DISTANCE_H

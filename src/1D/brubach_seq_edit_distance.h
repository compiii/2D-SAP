#ifndef BRUBACH_SEQ_EDIT_DISTANCE_H
#define BRUBACH_SEQ_EDIT_DISTANCE_H

#include "duo.h"
#include "bitwise_shorcuts.h"

#include <vector>
#include <map>
#include <unordered_map>

using brubach_lookup_table_t = std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<vlong, vlong>>>;

struct row_column
{
    std::vector<int> row;
    std::vector<int> column;
};

row_column brubach_seq_edit_distance(
    const std::string &a, const std::string &b,
    const row_column &initial_values,
    brubach_lookup_table_t &lookup_table, int &t, double &table_build_duration, double &efficiency);

std::vector<int> brubach_seq_edit_distance(
    const std::string &a, const std::string &b,
    const std::vector<int> &initial_values,
    brubach_lookup_table_t &lookup_table,
    double &table_build_duration, double &efficiency);

int brubach_seq_edit_distance(
    const std::string &a, const std::string &b,
    double &table_build_duration, double &efficiency);

int solution(const std::vector<std::vector<duo>> &scores);

#endif // BRUBACH_SEQ_EDIT_DISTANCE_H

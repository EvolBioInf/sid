#ifndef PROFILES_H
#define PROFILES_H

#include <array>
#include <cstdint>
#include <iostream>

// [#A, #C, #G, #T, #(A+C+G+T)]
using Profile = std::array<uint16_t, 5>;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;
const int COV = 4;

namespace std {
    std::ostream& operator<<(std::ostream& os, const Profile& p);
}

Profile parseRead(const char* read, char reference = 'n');

#endif

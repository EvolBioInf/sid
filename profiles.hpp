#ifndef PROFILES_H
#define PROFILES_H

#include <array>
#include <string>
#include <iostream>

// [#A, #C, #G, #T, #(A+C+G+T)]
using Profile = std::array<unsigned long int, 5>;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;
const int COV = 4;

Profile operator+(const Profile& lhs, const Profile& rhs);
Profile operator*(const Profile& lhs, const int mult);

namespace std {
    std::ostream& operator<<(std::ostream& os, const Profile& p);
}

Profile parseRead(const std::string read, char reference = 'n');

#endif

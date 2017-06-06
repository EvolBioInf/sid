#ifndef PROFILES_H
#define PROFILES_H

#include <array>
#include <string>
#include <iostream>

// [#A, #C, #G, #T, #(A+C+G+T)]
using Profile = std::array<int, 5>;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;
const int COV = 4;

Profile parseRead(const std::string read);

std::ostream& operator<<(std::ostream& os, const Profile& p);

namespace std {
    template<> struct hash<Profile> {
        size_t operator()(const Profile& p) const;
    };
}

//class Profile {
//public:
//    int a, c, g, t;
//    Profile(int aa, int cc, int gg, int tt);
//    Profile(const std::string read);
//    bool operator==(const Profile& other) const;
//    int operator[](const int index);
//};



#endif

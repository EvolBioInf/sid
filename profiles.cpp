#include <cstring>
#include <string>
#include <iostream>

#include "profiles.hpp"

using namespace std;

Profile parseRead(const string read) {
    Profile p {0, 0, 0, 0, 0};
    for(auto i = read.begin(); i != read.end(); ++i) {
        switch (*i) {
            case 'a': 
            case 'A': 
                ++p[A]; break;
            case 'c': 
            case 'C': 
                ++p[C]; break;
            case 'g': 
            case 'G': 
                ++p[G]; break;
            case 't': 
            case 'T': 
                ++p[T]; break;
            case '^':
                // skip next char
                ++i;
                break;
            case '+':
            case '-': {
                // parse following number, which indicates a range of insert/del bases
                size_t first_after_number = 0;
                int length = stoi(string(i+1, read.end()), &first_after_number);

                // skip parsed number + that number of bases after the number
                i += first_after_number + length;
                break;
            } default:
                continue;
        }
    }
    p[COV] = p[A] + p[C] + p[G] + p[T];
    return p;
}

ostream& operator<<(ostream& os, const Profile& p) {
    return os << "<" << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << ">";
}

size_t hash<Profile>::operator()(const Profile& p) const {
    return hash<int>{}(p[0]) ^ (hash<int>{}(p[1]) << 1) ^ (hash<int>{}(p[2]) << 2) ^ (hash<int>{}(p[3]) << 3);
}

/*Profile::Profile(int aa, int cc, int gg, int tt) : a{aa}, c{cc}, g{gg}, t{tt} {}

Profile::Profile(const string read) : a{0}, c{0}, g{0}, t{0} {
    for(auto i = read.begin(); i != read.end(); ++i) {
        switch (*i) {
            case 'a': 
            case 'A': 
                ++a; break;
            case 'c': 
            case 'C': 
                ++c; break;
            case 'g': 
            case 'G': 
                ++g; break;
            case 't': 
            case 'T': 
                ++t; break;
            case '^':
                // skip next char
                ++i;
                break;
            case '+':
            case '-': {
                // parse following number, which indicates a range of insert/del bases
                size_t first_after_number = 0;
                int length = stoi(string(i+1, read.end()), &first_after_number);

                // skip parsed number + that number of bases after the number
                i += first_after_number + length;
                break;
            } default:
                continue;
        }
    }
}

bool Profile::operator==(const Profile& other) const {
    return a == other.a && c == other.c && g == other.g && t == other.t;
}

int Profile::operator[](const int index) {
    switch (index) {
        case 0: return a;
        case 1: return c;
        case 2: return g;
        case 3: return t;
        default: return 0;
    }
}

*/

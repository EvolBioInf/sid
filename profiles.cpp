#include <cstdlib>
#include <cstring>
#include <iostream>

#include "profiles.hpp"

using namespace std;

bool inline tryIncrementBaseCount(char base, Profile& p) {
    switch (base) {
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
        default:
            return false;
    }
    return true;
}

Profile parseRead(const char* read, char reference) {
    Profile p {0, 0, 0, 0, 0};
    for(const char* base = read; *base != '\0'; ++base) {
        if (!tryIncrementBaseCount(*base, p)) {
            switch (*base) {
                case '^':
                    // skip next char
                    ++base;
                    break;
                case '+':
                case '-': {
                    // parse following number, which indicates a range of insert/del bases
                    char* first_after_number;
                    int length = strtol(base + 1, &first_after_number, 10);

                    // skip parsed number + that number of bases after the number
                    base = first_after_number + length - 1;
                    break;
                }
                case '.':
                case ',': {
                    tryIncrementBaseCount(reference, p);
                    break;
                }
                default:
                    break;
            }
        }
    }
    p[COV] = p[A] + p[C] + p[G] + p[T];
    return p;
}

ostream& std::operator<<(ostream& os, const Profile& p) {
    return os << p[0] << "." << p[1] << "." << p[2] << "." << p[3];
}

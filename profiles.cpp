#include <cstring>
#include <string>
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

Profile parseRead(const string& read, char reference) {
    Profile p {0, 0, 0, 0, 0};
    for(auto i = read.begin(); i != read.end(); ++i) {
        if (!tryIncrementBaseCount(*i, p)) {
            switch (*i) {
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

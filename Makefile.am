
SUBDIRS = . test

bin_PROGRAMS = sid sid2

sid_SOURCES = sid.cpp likelihoods.cpp likelihoods.hpp profiles.cpp profiles.hpp stats.hpp stats.cpp optimization.hpp
sid_CPPFLAGS = -std=c++11 -pedantic -Wall -Wextra -Wconversion

sid2_SOURCES = sid2.cpp call.hpp pileup.hpp lynch.hpp call.cpp pileup.cpp lynch.cpp stats.hpp stats.cpp
sid2_CPPFLAGS = -std=c++14 -pedantic -Wall -Wextra -Wconversion -fopenmp
sid2_LDFLAGS = -fopenmp

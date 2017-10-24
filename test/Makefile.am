AUTOMAKE_OPTIONS = serial-tests
TESTS = test test2
check_PROGRAMS = test test2

test_CPPFLAGS = -std=c++14 -I$(top_srcdir)
test_SOURCES = test.cpp test-profiles.cpp ../profiles.cpp

test2_CPPFLAGS = -std=c++14 -I$(top_srcdir)
test2_SOURCES = test.cpp
test2_SOURCES += test-pileup_parser.cpp ../pileup_parser.cpp
test2_SOURCES += test-call.cpp ../call.cpp

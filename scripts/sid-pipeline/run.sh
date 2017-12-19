#!/bin/bash
./parallel-run-sid.sh &&
./filter-exon-sites.sh &&
./identify-nonsyn.sh

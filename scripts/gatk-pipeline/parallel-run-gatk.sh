#!/bin/bash
parallel -u -t "$@" "./run-gatk.sh" ::: {1..19} X Y MT

#!/bin/bash
parallel -u -t "$@" "./run-sam.sh" ::: {1..19} X Y MT

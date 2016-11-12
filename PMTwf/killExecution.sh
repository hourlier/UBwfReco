#!/bin/bash

myProcess=($(ps aux | grep -i "python mac/runPMTwf_ana.py" | grep -v grep))
#echo ${myProcess[1]}
kill -9 ${myProcess[1]}
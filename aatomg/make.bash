#!/bin/bash
# gcc version 4.4.6 20110731 (Red Hat 4.4.6-3) (GCC) 
gcc -O3 -o ./bin/gor       ./src/gor.c      -lm
gcc -O3 -o ./bin/aatobase  ./src/aatobase.c
gcc -O3 -o ./bin/edmg0     ./src/edmg0.c
gcc -O3 -o ./bin/edmg1     ./src/edmg1.c
gcc -O3 -o ./bin/edmg2     ./src/edmg2.c
gcc -O3 -o ./bin/microgene ./src/microgene.c
gcc -O3 -o ./bin/ssscore   ./src/ssscore.c
gcc -O3 -o ./bin/uniqaa    ./src/uniqaa.c
gcc -O3 -o ./bin/scoreab   ./src/scoreab.c

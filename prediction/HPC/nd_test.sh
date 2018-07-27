#!/bin/sh
# options
#$ -m abe
#$ -M mtrachsel@wisc.edu
#$ -V
#$ -R y
#$ -q long


./cal_g num_sample = 50 num_warmup = 50 input file = elicitation_neus_certainty_median.dump output file = stepps_test.csv





#!/bin/sh
# options
#$ -N Matt
#$ -m be
#$ -M mtrachsel@wisc.edu
#$ -pe smp 15
#$ -V
#$ -R y
#$ -q long

./veg_od_mpp_nb_262  sample num_samples=2000 num_warmup=250 data file=veg_data_15_taxa_6796_cells_260_knots.dump output file=stepps_veg1.csv



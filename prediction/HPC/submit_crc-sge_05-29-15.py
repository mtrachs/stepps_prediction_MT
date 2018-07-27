
import subprocess

runs = [ ('pred_kw_kgamma_PL',
          './pred_kw_kgamma_262.exe \
          sample num_warmup=75 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_0to2000ypb_PL_umw_3by_v0.3.dump \
          output file=../output/12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v3.dump \
          random seed=42'),
         ('pred_kw_kgamma_G',
          './pred_kw_kgamma_262.exe \
          sample num_warmup=75 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v0.3.dump \
          output file=../output/12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v3.dump \
          random seed=42')
]

         # data file=../r/dump/12taxa_6341cells_120knots_0to2000ypb_PL_umw_1by_v0.3.dump \
         #  output file=../output/12taxa_6341cells_120knots_0to2000ypb_PL_umw_1by_v0.3.csv \

qsub = """
#!/bin/sh
# options
#$ -N {name}
#$ -m be
#$ -M andria.dawson@gmail.com
#$ -pe smp {threads}
#$ -V
#$ -R y
#$ -q *@@bio

cd $HOME/Documents/projects/stepps-prediction/cpp

fsync $SGE_STDOUT_PATH &
fsync ../output/12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v3.dump &
fsync ../output/12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v3.dump &

export OMP_NUM_THREADS={threads}
{command}
"""

dry_run = False

for name, command in runs:
#    sub = qsub.format(queue="low.q", walltime="672:00:00", command=run, threads=1)
    sub = qsub.format(command=command, threads=12, name=name)
    with open(name + ".qsub", 'w') as f:
        f.write(sub)
    print "submitting:", name
    if not dry_run:
        subprocess.check_call(['qsub', name + '.qsub'])
    else:
        print sub

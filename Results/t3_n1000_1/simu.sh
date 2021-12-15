#!/bin/bash
#PBS -l nodes=1:ppn=32,walltime=99:00:00
#PBS -N t3_n1000_1
#PBS -l mem=120gb
#PBS -o jobout/out.out
#PBS -e jobout/out.err
#PBS -q sph_sbasu

module load compilers/r-4.1.0
cd /home/jtu22/VBJM3/t3_n1000_1/

N=50
for ii in {1..100}; do
        # .. do your stuff here
	((i=ii%N)); ((i++==0)) && wait
        echo "starting task $ii.."
        Rscript runSimu2.R $ii &
done

# wait for pending jobs
wait
echo "all done"

##Rscript test.R $PBS_ARRAYID Rscript test.R $PBS_ARRAYID


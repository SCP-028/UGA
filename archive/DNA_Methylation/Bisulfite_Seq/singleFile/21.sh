#PBS -S /bin/bash
#PBS -q batch
#PBS -N 21
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=240:00:00
#PBS -l mem=4gb
#PBS -M jzhou0317@gmail.com
#PBS -m bae

cd $PBS_O_WORKDIR

#module load R/3.2.3
#R CMD BATCH data_download.R

/home/yz73026/packages/bin/python3.5 21.py

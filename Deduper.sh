#!/bin/bash


#SBATCH --account=bgmp                    
#SBATCH --partition=bgmp               
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32GB                        

conda activate bgmp_star

/usr/bin/time ./VanGordon_deduper.py -f "./C1_SE_uniqAlign.sorted.sam" -o "./results/DeduplicatedC1_SE_uniqAlign.sam" -u "./STL96.txt" -s "./results/SummaryDeduplicated.txt" -d "./results/Duplicates.sam"
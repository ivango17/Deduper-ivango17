#!/bin/bash


#SBATCH --account=bgmp                    
#SBATCH --partition=compute               
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32GB                        

conda activate base

/usr/bin/time ./VanGordon_deduper.py -f "./C1_SE_uniqAlign.sorted.sam" -o "./DeduplicatedC1_SE_uniqAlign.sam" -u "./STL96.txt" -s "./SummaryDeduplicated.txt"
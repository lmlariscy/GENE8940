#!/bin/bash
#SBATCH --job-name=simulate_reads_rsv	               
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=4		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=2:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL 

#set working directory
cd /work/gene8940/lml38336/project/references/rsv/primers

#primer pair 1
echo ">RSVA_F1" >> rsv_1.fas
echo "ACGCGAAAAAATGCGTACTACAAAC" >> rsv_1.fas
echo ">RSVA_R1" >> rsv_1.fas
echo "CTGMACCATAGGCATTCATAAACA" >> rsv_1.fas

#primer pair 2
echo ">RSVA_F2" >> rsv_2.fas
echo "ATGGGAGARGTRGCTCCAGAATA" >> rsv_2.fas
echo ">RSVA_R2" >> rsv_2.fas
echo "CGTGTAGCTGTRTGYTTCCAA" >> rsv_2.fas

#primer pair 3
echo ">RSVA_F3" >> rsv_3.fas
echo "GCTATGGCAAGACTYAGGAATG" >> rsv_3.fas
echo ">RSVA_R3" >> rsv_3.fas
echo "TTGAGRTCTAACACTTTGCTGGT" >> rsv_3.fas

#primer pair 4
echo ">RSVA_F4" >> rsv_4.fas
echo "AGCAAATTYTGGCCYTAYTTTAC" >> rsv_4.fas
echo ">RSVA_R4" >> rsv_4.fas
echo "CTCATAGCAACACATGCTGATTG" >> rsv_4.fas

#primer pair 5
echo ">RSVA_F5" >> rsv_5.fas
echo "TGATGCATCAATATCTCAAGTCA" >> rsv_5.fas
echo ">RSVA_R5" >> rsv_5.fas
echo "GRCCTATDCCTGCATACTC" >> rsv_5.fas

#primer pair 6
echo ">RSVA_F6" >> rsv_6.fas
echo "TGGACCATWGAAGCYATATCA" >> rsv_6.fas
echo ">RSVA_R6" >> rsv_6.fas
echo "AGTGTCAAAAACTAATRTCTCGT" >> rsv_6.fas
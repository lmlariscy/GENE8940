#!/bin/bash
#SBATCH --job-name=create_rsv_primer_files              
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=4		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=2:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL 

#set working directory
cd /work/gene8940/lml38336/project/references/covid/primers

#1
echo ">SARS-COV-2_F1" >> covid_1.fas
echo "AACAAACCAACCAACTTTCGATCTC" >> covid_1.fas
echo ">SARS-COV-2_R#" >> covid_#.fas
echo "CTTCTACTAAGCCACAAGTGCCA" >> covid_1.fas

#2
echo ">SARS-COV-2_F2" >> covid_2.fas
echo "TTTACAGGTTCGCGACGTGC" >> covid_2.fas
echo ">SARS-COV-2_R2" >> covid_2.fas
echo "ATAAGGATCAGTGCCAAGCTCG" >> covid_2.fas

#3
echo ">SARS-COV-2_F3" >> covid_3.fas
echo "GTAATAAAGGAGCTGGTGGCCA" >> covid_3.fas
echo ">SARS-COV-2_R3" >> covid_3.fas
echo "GCCAATTTAATTTCAAAAGGTGTCTGC" >> covid_3.fas

#4
echo ">SARS-COV-2_F4" >> covid_4.fas
echo "GTGTATACTGCTGCCGTGAACA" >> covid_4.fas
echo ">SARS-COV-2_R4" >> covid_4.fas
echo "ACAACAGCATTTTGGGGTAAGTAAC" >> covid_4.fas

#5
echo ">SARS-COV-2_F5" >> covid_5.fas
echo "TGAAACTTCATGGCAGACGGG" >> covid_5.fas
echo ">SARS-COV-2_R5" >> covid_5.fas
echo "TTGATGTTGACTTTCTCTTTTTGGAGT" >> covid_5.fas

#6
echo ">SARS-COV-2_F6" >> covid_6.fas
echo "CGTGCTAGCGCTAACATAGGTT" >> covid_6.fas
echo ">SARS-COV-2_R6" >> covid_6.fas
echo "AACACGCACAGAATTTTGAGCAG" >> covid_6.fas

#7
echo ">SARS-COV-2_F7" >> covid_7.fas
echo "ACTGAGTCCTCTTTATGCATTTGC" >> covid_7.fas
echo ">SARS-COV-2_R7" >> covid_7.fas
echo "CCACCGACAATTTCACAAGCAC" >> covid_7.fas

#8
echo ">SARS-COV-2_F8" >> covid_8.fas
echo "GCTTGAAGAGAAGTTTAAGGAAGGTG" >> covid_8.fas
echo ">SARS-COV-2_R8" >> covid_8.fas
echo "GGTTGTTCTAATGGTTGTAAATCACCA" >> covid_8.fas

#9
echo ">SARS-COV-2_F9" >> covid_9.fas
echo "TCTTCTTAGAGGGAGAAACACTTCC" >> covid_9.fas
echo ">SARS-COV-2_R9" >> covid_9.fas
echo "CACAGGCGAACTCATTTACTTCTG" >> covid_9.fas

#10
echo ">SARS-COV-2_F10" >> covid_10.fas
echo "TGAATATCACTTTTGAACTTGATGAAAGGATTG" >> covid_10.fas
echo ">SARS-COV-2_R10" >> covid_10.fas
echo "GGTTGAAGAGCAGCAGAAGTG" >> covid_10.fas

#11
echo ">SARS-COV-2_F11" >> covid_11.fas
echo "AGAAGAGTTTGAGCCATCAACTCA" >> covid_11.fas
echo ">SARS-COV-2_R11" >> covid_11.fas
echo "TTTAAGGCTCCTGCAACACCTC" >> covid_11.fas

#12
echo ">SARS-COV-2_F12" >> covid_12.fas
echo "TGCAGACATTGTGGAAGAAGCT" >> covid_12.fas
echo ">SARS-COV-2_R12" >> covid_12.fas
echo "CAGCTAAGTAGACATTTGTGCGAAC" >> covid_12.fas

#13
echo ">SARS-COV-2_F13" >> covid_13.fas
echo "AGCACGAAGTTCTACTTGCACC" >> covid_13.fas
echo ">SARS-COV-2_R13" >> covid_13.fas
echo "GATGTCAATGTCACTAACAAGAGTGG" >> covid_13.fas

#14
echo ">SARS-COV-2_F14" >> covid_14.fas
echo "TGGAAGAAACTAAGTTCCTCACAGAA" >> covid_14.fas
echo ">SARS-COV-2_R14" >> covid_14.fas
echo "CATGTGCAAGCATTTCTCGCAA" >> covid_14.fas

#15
echo ">SARS-COV-2_F15" >> covid_15.fas
echo "AAAAGTGCCTTTTACATTCTACCATCT" >> covid_15.fas
echo ">SARS-COV-2_R15" >> covid_15.fas
echo "GCATCAGGTGAAGAAACAGAAACTG" >> covid_15.fas

#16
echo ">SARS-COV-2_F16" >> covid_16.fas
echo "TGTAACACATGGCTTAAATTTGGAAGAA" >> covid_16.fas
echo ">SARS-COV-2_R16" >> covid_16.fas
echo "CACAACTTGCGTGTGGAGGTTA" >> covid_16.fas

#17
echo ">SARS-COV-2_F17" >> covid_17.fas
echo "TGACAATCTTAAGACACTTCTTTCTTTGAG" >> covid_17.fas
echo ">SARS-COV-2_R17" >> covid_17.fas
echo "TTCAACTCTATTTGTTGGAGTGTTAACAA" >> covid_17.fas

#18
echo ">SARS-COV-2_F18" >> covid_18.fas
echo "TGGAAATACCCACAAGTTAATGGTTTAAC" >> covid_18.fas
echo ">SARS-COV-2_R18" >> covid_18.fas
echo "GCTTGTTTACCACACGTACAAGG" >> covid_18.fas

#19
echo ">SARS-COV-2_F19" >> covid_19.fas
echo "AAGCTGTTATGTACATGGGCACA" >> covid_19.fas
echo ">SARS-COV-2_R19" >> covid_19.fas
echo "TGTCCAACTTAGGGTCAATTTCTGT" >> covid_19.fas

#20
echo ">SARS-COV-2_F20" >> covid_20.fas
echo "ACAAAGAAAACAGTTACACAACAACCA" >> covid_20.fas
echo ">SARS-COV-2_R20" >> covid_20.fas
echo "ACGTGGCTTTATTAGTTGCATTGTT" >> covid_20.fas

#21
echo ">SARS-COV-2_F21" >> covid_21.fas
echo "CACTACACACCCTCTTTTAAGAAAGG" >> covid_21.fas
echo ">SARS-COV-2_R21" >> covid_21.fas
echo "GTAAGACTAGAATTGTCTACATAAGCAGC" >> covid_21.fas

#22
echo ">SARS-COV-2_F22" >> covid_22.fas
echo "GTAGGAGACATTATACTTAAACCAGCAAA" >> covid_22.fas
echo ">SARS-COV-2_R22" >> covid_22.fas
echo "CCGACACTCTTAACAGTATTCTTTGC" >> covid_22.fas

#23
echo ">SARS-COV-2_F23" >> covid_23.fas
echo "AAACCGTGTTTGTACTAATTATATGCCTT" >> covid_23.fas
echo ">SARS-COV-2_R23" >> covid_23.fas
echo "AGAATCTAAACCACTAAGACAAACACTAC" >> covid_23.fas

#24
echo ">SARS-COV-2_F24" >> covid_24.fas
echo "GGTTACAGAGAAGGCTATTTGAACTCT" >> covid_24.fas
echo ">SARS-COV-2_R24" >> covid_24.fas
echo "ACAACATGCACATAACTTTTCCATACA" >> covid_24.fas

#25
echo ">SARS-COV-2_F25" >> covid_25.fas
echo "CAAATGGCCCCGATTTCAGCTA" >> covid_25.fas
echo ">SARS-COV-2_R25" >> covid_25.fas
echo "TGGATGGAACCATTCTTCACTGT" >> covid_25.fas

#26
echo ">SARS-COV-2_F26" >> covid_26.fas
echo "GCGAGAGACTTGTCACTACAGTT" >> covid_26.fas
echo ">SARS-COV-2_R26" >> covid_26.fas
echo "GAGTTTTTCCATTGGTACGTTAAAAGTTG" >> covid_26.fas

#27
echo ">SARS-COV-2_F27" >> covid_27.fas
echo "CTGATGTTGGTGATAGTGCGGA" >> covid_27.fas
echo ">SARS-COV-2_R27" >> covid_27.fas
echo "AATGTTGTGACTTTTTGCTACCTGC" >> covid_27.fas

#28
echo ">SARS-COV-2_F28" >> covid_28.fas
echo "TGAAAACATGACACCCCGTGAC" >> covid_28.fas
echo ">SARS-COV-2_R28" >> covid_28.fas
echo "TGACACCACCATCAATAGCCTTG" >> covid_28.fas

#29
echo ">SARS-COV-2_F29" >> covid_29.fas
echo "CTTGTGTTCCTTTTTGTTGCTGC" >> covid_29.fas
echo ">SARS-COV-2_R29" >> covid_29.fas
echo "AGCCAAAACACAAGCTGATGTTG" >> covid_29.fas

#30
echo ">SARS-COV-2_F30" >> covid_30.fas
echo "ACCTAGAGTTTTTAGTGCAGTTGGT" >> covid_30.fas
echo ">SARS-COV-2_R30" >> covid_30.fas
echo "CTACACCACAGAAAACTCCTGGT" >> covid_30.fas

#31
echo ">SARS-COV-2_F31" >> covid_31.fas
echo "CCTTGAAGGTTCTGTTAGAGTGGT" >> covid_31.fas
echo ">SARS-COV-2_R31" >> covid_31.fas
echo "AATGAGTAAACTGGTGTTAAACAGAGTAC" >> covid_31.fas

#32
echo ">SARS-COV-2_F32" >> covid_32.fas
echo "GAGCTTTTGGTGAATACAGTCATGTAG" >> covid_32.fas
echo ">SARS-COV-2_R32" >> covid_32.fas
echo "GAGGTAATAGCACATCACTACGCA" >> covid_32.fas

#33
echo ">SARS-COV-2_F33" >> covid_33.fas
echo "GTACTTTTGAAGAAGCTGCGCTG" >> covid_33.fas
echo ">SARS-COV-2_R33" >> covid_33.fas
echo "TGTCTTGGACAGTAAACTACGTCATC" >> covid_33.fas

#34
echo ">SARS-COV-2_F34" >> covid_34.fas
echo "TCCCATCTGGTAAAGTTGAGGGT" >> covid_34.fas
echo ">SARS-COV-2_R34" >> covid_34.fas
echo "CCACATGAACCATTAAGGAATGAACC" >> covid_34.fas

#35
echo ">SARS-COV-2_F35" >> covid_35.fas
echo "GTGTTAGCTTGTTACAATGGTTCACC" >> covid_35.fas
echo ">SARS-COV-2_R35" >> covid_35.fas
echo "AGGTCCTAGTATGTCAACATGGTCT" >> covid_35.fas

#36
echo ">SARS-COV-2_F36" >> covid_36.fas
echo "CAATCGATTTACCACAACTCTTAATGACT" >> covid_36.fas
echo ">SARS-COV-2_R36" >> covid_36.fas
echo "ACCCATAGCAAAAGGTAAAAAGGC" >> covid_36.fas

#37
echo ">SARS-COV-2_F37" >> covid_37.fas
echo "CACACCACTGGTTGTTACTCACA" >> covid_37.fas
echo ">SARS-COV-2_R37" >> covid_37.fas
echo "GTGTCAAGACATTCATAAGTGTCCAC" >> covid_37.fas

#38
echo ">SARS-COV-2_F38" >> covid_38.fas
echo "GACTGTGTTATGTATGCATCAGCTG" >> covid_38.fas
echo ">SARS-COV-2_R38" >> covid_38.fas
echo "CCTGTGTAGAAACTAAGTAATCATAAACACC" >> covid_38.fas

#39
echo ">SARS-COV-2_F39" >> covid_39.fas
echo "GCTATTTTTGTACTTGTTACTTTGGCC" >> covid_39.fas
echo ">SARS-COV-2_R39" >> covid_39.fas
echo "CCCTGCATGGAAAGCAAAACAG" >> covid_39.fas

#40
echo ">SARS-COV-2_F40" >> covid_40.fas
echo "TGTCCAGTTACACAATGACATTCTCT" >> covid_40.fas
echo ">SARS-COV-2_R40" >> covid_40.fas
echo "ACTTTTGCCCTCTTGTCCTCAG" >> covid_40.fas

#41
echo ">SARS-COV-2_F41" >> covid_41.fas
echo "ATTTGACCGTGATGCAGCCAT" >> covid_41.fas
echo ">SARS-COV-2_R41" >> covid_41.fas
echo "AAGAGGCCATGCTAAATTAGGTGAA" >> covid_41.fas

#42
echo ">SARS-COV-2_F42" >> covid_42.fas
echo "TGGTACAACATTTACTTATGCATCAGC" >> covid_42.fas
echo ">SARS-COV-2_R42" >> covid_42.fas
echo "TGTCTGTAACAAACCTACAAGGTGG" >> covid_42.fas

#43
echo ">SARS-COV-2_F43" >> covid_43.fas
echo "GGATTTGAAATGGGCTAGATTCCCT" >> covid_43.fas
echo ">SARS-COV-2_R43" >> covid_43.fas
echo "CGATGCACCACCAAAGGATTCT" >> covid_43.fas

#44
echo ">SARS-COV-2_F44" >> covid_44.fas
echo "GGGGACAACCAATCACTAATTGTG" >> covid_44.fas
echo ">SARS-COV-2_R44" >> covid_44.fas
echo "CATCAGTACTAGTGCCTGTGCC" >> covid_44.fas

#45
echo ">SARS-COV-2_F45" >> covid_45.fas
echo "TAAACGGGTTTGCGGTGTAAGT" >> covid_45.fas
echo ">SARS-COV-2_R45" >> covid_45.fas
echo "TCACAATTACCTTCATCAAAATGCCT" >> covid_45.fas

#46
echo ">SARS-COV-2_F46" >> covid_46.fas
echo "AGAATAGACGGTGACATGGTACC" >> covid_46.fas
echo ">SARS-COV-2_R46" >> covid_46.fas
echo "TCTACAACAGGAACTCCACTACCT" >> covid_46.fas

#47
echo ">SARS-COV-2_F47" >> covid_47.fas
echo "TGGTGTACTGACATTAGATAATCAAGATCT" >> covid_47.fas
echo ">SARS-COV-2_R47" >> covid_47.fas
echo "TGGAACACCATCAACAAATATTTTTCTCA" >> covid_47.fas

#48
echo ">SARS-COV-2_F48" >> covid_48.fas
echo "ACTGTTTGGATGACAGATGCATTC" >> covid_48.fas
echo ">SARS-COV-2_R48" >> covid_48.fas
echo "CAGAACTTCCTTCCTTAAAGAAACCC" >> covid_48.fas

#49
echo ">SARS-COV-2_F49" >> covid_49.fas
echo "ACAATGTTGCTTTTCAAACTGTCAAAC" >> covid_49.fas
echo ">SARS-COV-2_R49" >> covid_49.fas
echo "GGGATGACATTACGTTTTGTATATGCG" >> covid_49.fas

#50
echo ">SARS-COV-2_F50" >> covid_50.fas
echo "CATTTAATAAATGGGGTAAGGCTAGACTTT" >> covid_50.fas
echo ">SARS-COV-2_R50" >> covid_50.fas
echo "GAGCAAGAACAAGTGAGGCCAT" >> covid_50.fas

#51
echo ">SARS-COV-2_F51" >> covid_51.fas
echo "GCAAATTCTATGGTGGTTGGCAC" >> covid_51.fas
echo ">SARS-COV-2_R51" >> covid_51.fas
echo "GTCTGTGTTGTAAATTGCGGACA" >> covid_51.fas

#52
echo ">SARS-COV-2_F52" >> covid_52.fas
echo "CTGTCACGGCCAATGTTAATGC" >> covid_52.fas
echo ">SARS-COV-2_R52" >> covid_52.fas
echo "GGATCTGGGTAAGGAAGGTACACA" >> covid_52.fas

#53
echo ">SARS-COV-2_F53" >> covid_53.fas
echo "ACTAAAGGACCTCATGAATTTTGCTC" >> covid_53.fas
echo ">SARS-COV-2_R53" >> covid_53.fas
echo "GCAAAGAACACAAGCCCCAAC" >> covid_53.fas

#54
echo ">SARS-COV-2_F54" >> covid_54.fas
echo "ACATGATGAGTTAACAGGACACATG" >> covid_54.fas
echo ">SARS-COV-2_R54" >> covid_54.fas
echo "CCAAAAACTTGTCCATTAGCACACA" >> covid_54.fas

#55
echo ">SARS-COV-2_F55" >> covid_55.fas
echo "AATGCTCCAGGTTGTGATGTCA" >> covid_55.fas
echo ">SARS-COV-2_R55" >> covid_55.fas
echo "ACACGATAACCAGTAAAGACATAATTTCG" >> covid_55.fas

#56
echo ">SARS-COV-2_F56" >> covid_56.fas
echo "ACTGTACGTGAAGTGCTGTCTG" >> covid_56.fas
echo ">SARS-COV-2_R56" >> covid_56.fas
echo "TGACTCTTACCAGTACCAGGTGG" >> covid_56.fas

#57
echo ">SARS-COV-2_F57" >> covid_57.fas
echo "GGCTTATACCCAACACTCAATATCTCA" >> covid_57.fas
echo ">SARS-COV-2_R57" >> covid_57.fas
echo "CTGGCATTGACAACACTCAAATCA" >> covid_57.fas

#58
echo ">SARS-COV-2_F58" >> covid_58.fas
echo "TGCCTGAGACGACAGCAGATAT" >> covid_58.fas
echo ">SARS-COV-2_R58" >> covid_58.fas
echo "TGTGGCCTGTTAATTGCAGATGA" >> covid_58.fas

#59
echo ">SARS-COV-2_F59" >> covid_59.fas
echo "GCTTAAAGCACATAAAGACAAATCAGC" >> covid_59.fas
echo ">SARS-COV-2_R59" >> covid_59.fas
echo "TCCTACGTGGAATTTCAAGACTTGT" >> covid_59.fas

#60
echo ">SARS-COV-2_F60" >> covid_60.fas
echo "ACAGATTTAATGTTGCTATTACCAGAGC" >> covid_60.fas
echo ">SARS-COV-2_R60" >> covid_60.fas
echo "TAGCATGACACCCCTCGACAT" >> covid_60.fas

#61
echo ">SARS-COV-2_F61" >> covid_61.fas
echo "ACCCTAACATGTTTATCACCCGC" >> covid_61.fas
echo ">SARS-COV-2_R61" >> covid_61.fas
echo "GCTCAGGTCCTATTTTCACAAAATACTT" >> covid_61.fas

#62
echo ">SARS-COV-2_F62" >> covid_62.fas
echo "GTGACACACTTAAAAATCTCTCTGACAG" >> covid_62.fas
echo ">SARS-COV-2_R62" >> covid_62.fas
echo "CCGCATTAATCTTCAGTTCATCACC" >> covid_62.fas

#63
echo ">SARS-COV-2_F63" >> covid_63.fas
echo "TAGGTGTCTAGCTGTCCACGAG" >> covid_63.fas
echo ">SARS-COV-2_R63" >> covid_63.fas
echo "CCAGGCAAGTTAAGGTTAGATAGCA" >> covid_63.fas

#64
echo ">SARS-COV-2_F64" >> covid_64.fas
echo "GCCTATTTTGGAATTGCAATGTCGA" >> covid_64.fas
echo ">SARS-COV-2_R64" >> covid_64.fas
echo "GTATCAAATTGTTTGTAAACCCACAAGC" >> covid_64.fas

#65
echo ">SARS-COV-2_F65" >> covid_65.fas
echo "GTCTGTAGACATCATGCTAATGAGTACA" >> covid_65.fas
echo ">SARS-COV-2_R65" >> covid_65.fas
echo "GCTGGAGCATCTCTTTTGTAGTCC" >> covid_65.fas

#66
echo ">SARS-COV-2_F66" >> covid_66.fas
echo "AACCAGTACCAGAGGTGAAAATACTC" >> covid_66.fas
echo ">SARS-COV-2_R66" >> covid_66.fas
echo "TTTCTACTCTGAGTAAAGTAAGTTTCAGGT" >> covid_66.fas

#67
echo ">SARS-COV-2_F67" >> covid_67.fas
echo "CAAACAAGCTAGTCTTAATGGAGTCAC" >> covid_67.fas
echo ">SARS-COV-2_R67" >> covid_67.fas
echo "AACACACACACTTAGATGAACCTGT" >> covid_67.fas

#68
echo ">SARS-COV-2_F68" >> covid_68.fas
echo "GACTAGCTAAACGTTTTAAGGAATCACC" >> covid_68.fas
echo ">SARS-COV-2_R68" >> covid_68.fas
echo "GCGACATTCATCATTATGCCTTTAGG" >> covid_68.fas

#69
echo ">SARS-COV-2_F69" >> covid_69.fas
echo "CGGGTGTTGCTATGCCTAATCT" >> covid_69.fas
echo ">SARS-COV-2_R69" >> covid_69.fas
echo "TTTGTAACATTTTTAGTCTTAGGGTCGTAC" >> covid_69.fas

#70
echo ">SARS-COV-2_F70" >> covid_70.fas
echo "TTGATTGGTGATTGTGCAACTGTAC" >> covid_70.fas
echo ">SARS-COV-2_R70" >> covid_70.fas
echo "AGAATAGGAAGACAACTGAATTGGATTTG" >> covid_70.fas

#71
echo ">SARS-COV-2_F71" >> covid_71.fas
echo "GGCAAACCACGCGAACAAATAG" >> covid_71.fas
echo ">SARS-COV-2_R71" >> covid_71.fas
echo "TGAGGATCTGAAAACTTTGTCAGGG" >> covid_71.fas

#72
echo ">SARS-COV-2_F72" >> covid_72.fas
echo "GTGATGTTCTTGTTAACAACTAAACGAAC" >> covid_72.fas
echo ">SARS-COV-2_R72" >> covid_72.fas
echo "GTAGCGTTATTAACAATAAGTAGGGACTG" >> covid_72.fas

#73
echo ">SARS-COV-2_F73" >> covid_73.fas
echo "AGAGGCTGGATTTTTGGTACTACT" >> covid_73.fas
echo ">SARS-COV-2_R73" >> covid_73.fas
echo "ACCTAGTGATGTTAATACCTATTGGCA" >> covid_73.fas

#74
echo ">SARS-COV-2_F74" >> covid_74.fas
echo "TGGACCTTGAAGGAAAACAGGG" >> covid_74.fas
echo ">SARS-COV-2_R74" >> covid_74.fas
echo "TGATAGATTCCTTTTTCTACAGTGAAGGA" >> covid_74.fas

#75
echo ">SARS-COV-2_F75" >> covid_75.fas
echo "GAAAATGGAACCATTACAGATGCTGT" >> covid_75.fas
echo ">SARS-COV-2_R75" >> covid_75.fas
echo "TTTGCCCTGGAGCGATTTGT" >> covid_75.fas

#76
echo ">SARS-COV-2_F76" >> covid_76.fas
echo "ATGTCTATGCAGATTCATTTGTAATTAGAGGT" >> covid_76.fas
echo ">SARS-COV-2_R76" >> covid_76.fas
echo "GTCCACAAACAGTTGCTGGTG" >> covid_76.fas

#77
echo ">SARS-COV-2_F77" >> covid_77.fas
echo "CAAACCTTTTGAGAGAGATATTTCAACTGA" >> covid_77.fas
echo ">SARS-COV-2_R77" >> covid_77.fas
echo "CACTGACACCACCAAAAGAACATG" >> covid_77.fas

#78
echo ">SARS-COV-2_F78" >> covid_78.fas
echo "CTGAGTCTAACAAAAAGTTTCTGCCTT" >> covid_78.fas
echo ">SARS-COV-2_R78" >> covid_78.fas
echo "GGATTGACTAGCTACACTACGTGC" >> covid_78.fas

#79
echo ">SARS-COV-2_F79" >> covid_79.fas
echo "ACCCATTGGTGCAGGTATATGC" >> covid_79.fas
echo ">SARS-COV-2_R79" >> covid_79.fas
echo "AATTGGTGGTGTTTTGTAAATTTGTTTGAC" >> covid_79.fas

#80
echo ">SARS-COV-2_F80" >> covid_80.fas
echo "CCGTGCTTTAACTGGAATAGCTG" >> covid_80.fas
echo ">SARS-COV-2_R80" >> covid_80.fas
echo "GCAAATGGTATTTGTAATGCAGCAC" >> covid_80.fas

#81
echo ">SARS-COV-2_F81" >> covid_81.fas
echo "TGCTCAATACACTTCTGCACTGT" >> covid_81.fas
echo ">SARS-COV-2_R81" >> covid_81.fas
echo "TGAAGTCTGCCTGTGATCAACC" >> covid_81.fas

#82
echo ">SARS-COV-2_F82" >> covid_82.fas
echo "TGCACAAGCTTTAAACACGCTT" >> covid_82.fas
echo ">SARS-COV-2_R82" >> covid_82.fas
echo "CACGAGGAAAGTGTGCTTTTCC" >> covid_82.fas

#83
echo ">SARS-COV-2_F83" >> covid_83.fas
echo "GCATGTGACTTATGTCCCTGCA" >> covid_83.fas
echo ">SARS-COV-2_R83" >> covid_83.fas
echo "AGATTCATTTAAATTCTTGGCAACCTCA" >> covid_83.fas

#84
echo ">SARS-COV-2_F84" >> covid_84.fas
echo "GTTGATTTAGGTGACATCTCTGGCA" >> covid_84.fas
echo ">SARS-COV-2_R84" >> covid_84.fas
echo "AGCATCCTTGATTTCACCTTGCT" >> covid_84.fas

#85
echo ">SARS-COV-2_F85" >> covid_85.fas
echo "ATGAAGACGACTCTGAGCCAGT" >> covid_85.fas
echo ">SARS-COV-2_R85" >> covid_85.fas
echo "CTGCAAGAAGTAGACTAAAGCATAAAGAT" >> covid_85.fas

#86
echo ">SARS-COV-2_F86" >> covid_86.fas
echo "TGTTGTTTGTAACAGTTTACTCACACC" >> covid_86.fas
echo ">SARS-COV-2_R86" >> covid_86.fas
echo "TCAATTGAGTTGAGTACAGCTGGT" >> covid_86.fas

#87
echo ">SARS-COV-2_F87" >> covid_87.fas
echo "GTGGTTATACTGAAAAATGGGAATCTGG" >> covid_87.fas
echo ">SARS-COV-2_R87" >> covid_87.fas
echo "AATCGAAGCGCAGTAAGGATGG" >> covid_87.fas

#88
echo ">SARS-COV-2_F88" >> covid_88.fas
echo "TTATGTACTCATTCGTTTCGGAAGAG" >> covid_88.fas
echo ">SARS-COV-2_R88" >> covid_88.fas
echo "ACAAAAACCTATTCCTGTTGGCATAG" >> covid_88.fas

#89
echo ">SARS-COV-2_F89" >> covid_89.fas
echo "TAGGTTTCCTATTCCTTACATGGATTTGT" >> covid_89.fas
echo ">SARS-COV-2_R89" >> covid_89.fas
echo "CTAGATGGTGTCCAGCAATACGAAG" >> covid_89.fas

#90
echo ">SARS-COV-2_F90" >> covid_90.fas
echo "ATTCTTCTCAACGTGCCACTCC" >> covid_90.fas
echo ">SARS-COV-2_R90" >> covid_90.fas
echo "ATTAGTAATATCTCTGCTATAGTAACCTGAAAG" >> covid_90.fas

#91
echo ">SARS-COV-2_F91" >> covid_91.fas
echo "TCCAGTAGCAGTGACAATATTGCTT" >> covid_91.fas
echo ">SARS-COV-2_R91" >> covid_91.fas
echo "AGTGCAAATTTGTTATCAGCTAGAGG" >> covid_91.fas

#92
echo ">SARS-COV-2_F92" >> covid_92.fas
echo "CACTACCAAGAGTGTGTTAGAGGTAC" >> covid_92.fas
echo ">SARS-COV-2_R92" >> covid_92.fas
echo "GTTCAAGTGAGAACCAAAAGATAATAAGC" >> covid_92.fas

#93
echo ">SARS-COV-2_F93" >> covid_93.fas
echo "TTGTTGCGGCAATAGTGTTTATAACA" >> covid_93.fas
echo ">SARS-COV-2_R93" >> covid_93.fas
echo "TGGGTGATTTAGAACCAGCCTC" >> covid_93.fas

#94
echo ">SARS-COV-2_F94" >> covid_94.fas
echo "ACCCGTGTCCTATTCACTTCTATTC" >> covid_94.fas
echo ">SARS-COV-2_R94" >> covid_94.fas
echo "TTATTGGGTAAACCTTGGGGCC" >> covid_94.fas

#95
echo ">SARS-COV-2_F95" >> covid_95.fas
echo "GTGCGTTGTTCGTTCTATGAAGAC" >> covid_95.fas
echo ">SARS-COV-2_R95" >> covid_95.fas
echo "ACCATCTTGGACTGAGATCTTTCATT" >> covid_95.fas

#96
echo ">SARS-COV-2_F96" >> covid_96.fas
echo "AGATGACCAAATTGGCTACTACCG" >> covid_96.fas
echo ">SARS-COV-2_R96" >> covid_96.fas
echo "CCATTGCCAGCCATTCTAGCA" >> covid_96.fas

#97
echo ">SARS-COV-2_F97" >> covid_97.fas
echo "TTCCTCATCACGTAGTCGCAAC" >> covid_97.fas
echo ">SARS-COV-2_R97" >> covid_97.fas
echo "CGACATTCCGAAGAACGCTGA" >> covid_97.fas

#98
echo ">SARS-COV-2_F98" >> covid_98.fas
echo "CCAGGAACTAATCAGACAAGGAACT" >> covid_98.fas
echo ">SARS-COV-2_R98" >> covid_98.fas
echo "TTTAGGCCTGAGTTGAGTCAGC" >> covid_98.fas

#99
echo ">SARS-COV-2_F99" >> covid_99.fas
echo "CTTCTTCCTGCTGCAGATTTGGA" >> covid_99.fas
echo ">SARS-COV-2_R99" >> covid_99.fas
echo "GCTATTAAAATCACATGGGGATAGCAC" >> covid_99.fas



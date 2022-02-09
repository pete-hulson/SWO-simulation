# SWO-simulation
Some code to perform length sampling

2 Scripts included here:
1. 'Sim Length Smplng.R' is the base code to perform resampling
2. 'est_SzPop.R' is the function that expands the length samples to estimates at pop'n at age

Structure of subfolders under sampling script:
- Data: holds files 'agepop.csv': RACE pop'n estimates at age, 'CPUE.csv': survey CPUE observations by haul, 'lfreq.csv': length frequency data, 'sizepop.csv': RACE pop'n estimates at size, 'specimen.csv': specimen dataset (with ages), 'strata.csv': file that defines strata of surveys
- Functions: where 'est_SzPop.R' would be placed
- Results: the script will create regional, then species-specific folders to write results

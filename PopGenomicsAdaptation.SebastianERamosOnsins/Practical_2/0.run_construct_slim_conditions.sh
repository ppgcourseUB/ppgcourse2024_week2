#### RUN DIFFERENT CONDITIONS. OBTAIN SFS IN SYN AND NONSYN. CALCULATE MKTasymptotic. ####
#### GENERAK CONDITIONS:
#### Two populations: ONE target plus ONE outgroup.
#### 1.The initial population run for 5Ne generations to achieve some equilibrium mutation-selection-drift
#### 2.Include NEGATIVE AND/OR POSITIVE SELECTION
#### 3.Split target and outgroup for 10Ne generations

rm ./run_slim_conditions.sh

#fixed paraneters
Ne=750; L=500000;
ind_sample_size=25;

# CONDITION 0:
#NEUTRAL.
FILEOUT="'./00_slim_SFS_SNM.txt'"
rate_ben=0; s_backg_ben=0;
rate_del=0; s_backg_del=0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim \& >> ./run_slim_conditions.sh

# CONDITION 1:
#NEGATIVE SELECTION.
FILEOUT="'./01_slim_SFS_NEG.txt'"
rate_ben=0;    s_backg_ben=0;
rate_del=0.85; s_backg_del=-0.005;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim \& >> ./run_slim_conditions.sh

# CONDITION 2:
#BENEFICIAL SELECTION.
FILEOUT="'./02_slim_SFS_POS.txt'"
rate_ben=0.1; s_backg_ben=0.002;
rate_del=0.0; s_backg_del=0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim \& >> ./run_slim_conditions.sh

# CONDITION 3:
#NEGATIVE SELECTION PLUS BENEFICIAL SELECTION.
FILEOUT="'./03_slim_SFS_NEG_POS.txt'"
rate_ben=0.10; s_backg_ben=+0.002;
rate_del=0.85; s_backg_del=-0.005;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim \& >> ./run_slim_conditions.sh

echo "sh ./run_slim_conditions.sh"
sh ./run_slim_conditions.sh

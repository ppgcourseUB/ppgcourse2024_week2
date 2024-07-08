# Neutral:
slim -t -m -d "freq_sel_init=0.000125" -d "freq_sel_end=1.0" -d "s_beneficial=0" -d "ind_sample_size=50" -d "file_output1='./00.slim.snm.practical1.output'" ./sweeps_ind_vcf.slim
# Complete Sweep:
slim -t -m -d "freq_sel_init=0.000125" -d "freq_sel_end=1.0" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./01.slim.selsweep.practical1.output'" ./sweeps_ind_vcf.slim
# Incomplete Sweep:
slim -t -m -d "freq_sel_init=0.000125" -d "freq_sel_end=0.75" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./02.slim.incomplete_selsweep.practical1.output'" ./sweeps_ind_vcf.slim
# Sweep from Standing variation:
slim -t -m -d "freq_sel_init=0.25" -d "freq_sel_end=1.0" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./03.slim.standing_selsweep.practical1.output'" ./sweeps_ind_vcf.slim

#plot selective trajectory
R --vanilla < ./run_curve_selection.R

#plot sliding windows heterozigosity
R --vanilla < ./run_sliding_windows_pi.R

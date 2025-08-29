./FES_from_Reweighting_multiT_funnel.py  --sigma 0.01 --colvar CVS_chn_unbiased_run1 --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_ref_rmsd_1.dat --max 0.85 --min 0.015;
./FES_from_Reweighting_multiT_funnel.py  --sigma 0.01 --colvar CVS_chn_unbiased_run2 --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_ref_rmsd_2.dat --max 0.85 --min 0.015;
./FES_from_Reweighting_multiT_funnel.py  --sigma 0.01 --colvar CVS_chn_unbiased_run3 --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_ref_rmsd_3.dat --max 0.85 --min 0.015;

grep 'DeltaF ' fes_ref_rmsd_?.dat | awk '{print -$4}' > deltaF.dat
awk '{values[NR] = $1; sum += $1} END {avg = sum / NR; for (i = 1; i <= NR; i++) {sumsq += (values[i] - avg) * (values[i] - avg);} stddev = sqrt(sumsq / (NR-1)); print avg, stddev}' deltaF.dat > deltaFall.dat;


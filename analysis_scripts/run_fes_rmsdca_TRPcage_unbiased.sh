rm fes_blocks*
./FES_from_Reweighting_multiT_funnel.py  --sigma 0.02 --colvar COLVAR_200us --cv rmsd_ca --bin 150 --temp 320 --deltaFat 0.5 --blocks 3 --out fes_ref_rmsd.dat --max 2.0 --min 0.03;
grep 'DeltaF ' fes_ref_rmsd_?.dat | awk '{print -$4}' > deltaF.dat
awk '{values[NR] = $1; sum += $1} END {avg = sum / NR; for (i = 1; i <= NR; i++) {sumsq += (values[i] - avg) * (values[i] - avg);} stddev = sqrt(sumsq / (NR-1)); print avg, stddev}' deltaF.dat > deltaFall.dat;

./FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR_200us --cv armsd19 --bin 120 --temp 320 --deltaFat 3.5 --blocks 3 --out fes_ref_armsd19.dat;
./FES_from_Reweighting_multiT_funnel.py --sigma 0.5 --colvar COLVAR_200us --cv diffHB_compact --bin 200 --temp 320 --deltaFat -10 --blocks 3 --out fes_ref_HB.dat --min -63 --max 7;
./FES_from_Reweighting_multiT_funnel.py --sigma 0.25 --colvar COLVAR_200us --cv cmap_compact --bin 200 --temp 320 --deltaFat 8 --blocks 3 --out fes_ref_SC.dat --min -15 --max 20;

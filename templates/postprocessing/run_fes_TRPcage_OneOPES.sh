

cd ../OneOPES_b30_a/0
rm fes_stride*
./../../FES_from_Reweighting.py  --sigma 0.02 --colvar COLVAR.* --cv rmsd_ca --bin 150 --temp 320 --deltaFat 0.5 --out fes_stride.dat --stride 20000 --skip 14362 --max 2.0 --min 0.03;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_b30_b/0
rm fes_stride*
./../../FES_from_Reweighting.py  --sigma 0.02 --colvar COLVAR.* --cv rmsd_ca --bin 150 --temp 320 --deltaFat 0.5 --out fes_stride.dat --stride 20000 --skip 20000 --max 2.0 --min 0.03;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_b30_c/0
rm fes_stride*
./../../FES_from_Reweighting.py  --sigma 0.02 --colvar COLVAR.* --cv rmsd_ca --bin 150 --temp 320 --deltaFat 0.5 --out fes_stride.dat --stride 20000 --skip 20000 --max 2.0 --min 0.03;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_b30_d/0
rm fes_stride*
./../../FES_from_Reweighting.py  --sigma 0.02 --colvar COLVAR.* --cv rmsd_ca --bin 150 --temp 320 --deltaFat 0.5 --out fes_stride.dat --stride 20000 --skip 20000 --max 2.0 --min 0.03;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_b30_e/0
rm fes_stride*
./../../FES_from_Reweighting.py  --sigma 0.02 --colvar COLVAR.* --cv rmsd_ca --bin 150 --temp 320 --deltaFat 0.5 --out fes_stride.dat --stride 20000 --skip 20000 --max 2.0 --min 0.03;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_b30

paste ../OneOPES_b30_a/0/deltaF.dat ../OneOPES_b30_b/0/deltaF.dat ../OneOPES_b30_c/0/deltaF.dat ../OneOPES_b30_d/0/deltaF.dat ../OneOPES_b30_e/0/deltaF.dat > deltaF_list.dat;

awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*20000+20000, sum}' deltaF_list.dat > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF-1; sum3 = sqrt(sum2); print sum3}' deltaF_list.dat > temp2;
paste temp1 temp2 > deltaFall.dat;
rm temp*

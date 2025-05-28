
cd a
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaF.dat
cd ../b
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaF.dat
cd ../c
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaF.dat
cd ../d
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaF.dat
cd ../e
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaF.dat
cd ..

paste a/delta* b/delta* c/delta* c/delta* e/delta* > deltaF_list.dat;

awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*20000+20000, sum}' deltaF_list.dat > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF-1; sum3 = sqrt(sum2); print sum3}' deltaF_list.dat > temp2;
paste temp1 temp2 > deltaFall.dat;
rm temp*

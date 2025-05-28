cd ../OneOPES_a/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv diffHB_compact --bin 150 --temp 340 --skip 20000 --deltaFat -10 --out fes_HB.dat --min -50 --max 5;
./../../FES_from_Reweighting.py --sigma 0.2 --colvar COLVAR.* --cv cmap_compact --bin 100 --temp 340 --skip 20000 --deltaFat 5 --blocks 3 --out fes_SC.dat --min -5 --max 8;
./../../FES_from_Reweighting.py --sigma 0.5,0.2 --colvar COLVAR.* --cv 5,6 --bin 150,100 --skip 20000 --temp 340 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-50, -5" --max 5,8;
tail -n 15402 fes_blocks_2d_HB_SC.dat > temp_fes_2d_a.dat
tail -n 151 fes_HB.dat > temp_fes_HB_a.dat
tail -n 101 fes_SC.dat > temp_fes_SC_a.dat
cp temp* ../../OneOPES/.
#
cd ../../OneOPES_b/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv diffHB_compact --bin 150 --temp 340 --skip 20000 --deltaFat -10 --out fes_HB.dat --min -50 --max 5;
./../../FES_from_Reweighting.py --sigma 0.2 --colvar COLVAR.* --cv cmap_compact --bin 100 --temp 340 --skip 20000 --deltaFat 5 --blocks 3 --out fes_SC.dat --min -5 --max 8;
./../../FES_from_Reweighting.py --sigma 0.5,0.2 --colvar COLVAR.* --cv 5,6 --bin 150,100 --skip 20000 --temp 340 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-50, -5" --max 5,8;
tail -n 15402 fes_blocks_2d_HB_SC.dat > temp_fes_2d_b.dat
tail -n 151 fes_HB.dat > temp_fes_HB_b.dat
tail -n 101 fes_SC.dat > temp_fes_SC_b.dat
cp temp* ../../OneOPES/.
#
cd ../../OneOPES_c/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv diffHB_compact --bin 150 --temp 340 --skip 20000 --deltaFat -10 --out fes_HB.dat --min -50 --max 5;
./../../FES_from_Reweighting.py --sigma 0.2 --colvar COLVAR.* --cv cmap_compact --bin 100 --temp 340 --skip 20000 --deltaFat 5 --blocks 3 --out fes_SC.dat --min -5 --max 8;
./../../FES_from_Reweighting.py --sigma 0.5,0.2 --colvar COLVAR.* --cv 5,6 --bin 150,100 --skip 20000 --temp 340 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-50, -5" --max 5,8;
tail -n 15402 fes_blocks_2d_HB_SC.dat > temp_fes_2d_c.dat
tail -n 151 fes_HB.dat > temp_fes_HB_c.dat
tail -n 101 fes_SC.dat > temp_fes_SC_c.dat
cp temp* ../../OneOPES/.
#
cd ../../OneOPES_d/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv diffHB_compact --bin 150 --temp 340 --skip 20000 --deltaFat -10 --out fes_HB.dat --min -50 --max 5;
./../../FES_from_Reweighting.py --sigma 0.2 --colvar COLVAR.* --cv cmap_compact --bin 100 --temp 340 --skip 20000 --deltaFat 5 --blocks 3 --out fes_SC.dat --min -5 --max 8;
./../../FES_from_Reweighting.py --sigma 0.5,0.2 --colvar COLVAR.* --cv 5,6 --bin 150,100 --skip 20000 --temp 340 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-50, -5" --max 5,8;
tail -n 15402 fes_blocks_2d_HB_SC.dat > temp_fes_2d_d.dat
tail -n 151 fes_HB.dat > temp_fes_HB_d.dat
tail -n 101 fes_SC.dat > temp_fes_SC_d.dat
cp temp* ../../OneOPES/.
#
cd ../../OneOPES_e/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv diffHB_compact --bin 150 --temp 340 --skip 20000 --deltaFat -10 --out fes_HB.dat --min -50 --max 5;
./../../FES_from_Reweighting.py --sigma 0.2 --colvar COLVAR.* --cv cmap_compact --bin 100 --temp 340 --skip 20000 --deltaFat 5 --blocks 3 --out fes_SC.dat --min -5 --max 8;
./../../FES_from_Reweighting.py --sigma 0.5,0.2 --colvar COLVAR.* --cv 5,6 --bin 150,100 --skip 20000 --temp 340 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-50, -5" --max 5,8;
tail -n 15402 fes_blocks_2d_HB_SC.dat > temp_fes_2d_e.dat
tail -n 151 fes_HB.dat > temp_fes_HB_e.dat
tail -n 201 fes_SC.dat > temp_fes_SC_e.dat
cp temp* ../../OneOPES/.

paste temp_fes_2d_a.dat temp_fes_2d_b.dat temp_fes_2d_c.dat temp_fes_2d_d.dat temp_fes_2d_e.dat | awk '{ if ($1 == "" || $2 == "" || $3 == "" || $4 == "") print ""; else print $1, $2, ($3+$7+$11+$15+$19)/5 }' > fes_2d_HB_SC_average.dat
paste temp_fes_HB_a.dat temp_fes_HB_b.dat temp_fes_HB_c.dat temp_fes_HB_d.dat temp_fes_HB_e.dat | awk '{ if ($1 == "" || $2 == "" || $3 == "") print ""; else print $1, ($2+$4+$6+$8+$10)/5 }' > fes_HB_average.dat
paste temp_fes_SC_a.dat temp_fes_SC_b.dat temp_fes_SC_c.dat temp_fes_SC_d.dat temp_fes_SC_e.dat | awk '{ if ($1 == "" || $2 == "" || $3 == "") print ""; else print $1, ($2+$4+$6+$8+$10)/5 }' > fes_SC_average.dat


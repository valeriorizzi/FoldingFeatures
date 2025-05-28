cd ../../OneOPES_b30_a/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv 5 --bin 200 --temp 320 --skip 20000 --out fes_HB.dat --min -63 --max 7 &
./../../FES_from_Reweighting.py --sigma 0.25 --colvar COLVAR.* --cv 6 --bin 200 --temp 320 --skip 20000 --out fes_SC.dat --min -15 --max 20 &
./../../FES_from_Reweighting.py --sigma 0.5,0.25 --colvar COLVAR.* --cv 5,6 --bin 200,200 --skip 20000 --temp 320 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-63, -15" --max 7,20;
tail -n 40602 fes_blocks_2d_HB_SC.dat > temp_fes_2d_g.dat
tail -n 201 fes_HB.dat > temp_fes_HB_g.dat
tail -n 201 fes_SC.dat > temp_fes_SC_g.dat
cp temp* ../../OneOPES_b30/.
cd ../../OneOPES_b30
#
cd ../OneOPES_b30_b/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv 5 --bin 200 --temp 320 --skip 20000 --out fes_HB.dat --min -63 --max 7 &
./../../FES_from_Reweighting.py --sigma 0.25 --colvar COLVAR.* --cv 6 --bin 200 --temp 320 --skip 20000 --out fes_SC.dat --min -15 --max 20 &
./../../FES_from_Reweighting.py --sigma 0.5,0.25 --colvar COLVAR.* --cv 5,6 --bin 200,200 --skip 20000 --temp 320 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-63, -15" --max 7,20;
tail -n 40602 fes_blocks_2d_HB_SC.dat > temp_fes_2d_b.dat
tail -n 201 fes_HB.dat > temp_fes_HB_b.dat
tail -n 201 fes_SC.dat > temp_fes_SC_b.dat
cp temp* ../../OneOPES_b30/.
#
cd ../../OneOPES_b30_c/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv 5 --bin 200 --temp 320 --skip 20000 --out fes_HB.dat --min -63 --max 7 &
./../../FES_from_Reweighting.py --sigma 0.25 --colvar COLVAR.* --cv 6 --bin 200 --temp 320 --skip 20000 --out fes_SC.dat --min -15 --max 20 &
./../../FES_from_Reweighting.py --sigma 0.5,0.25 --colvar COLVAR.* --cv 5,6 --bin 200,200 --skip 20000 --temp 320 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-63, -15" --max 7,20;
tail -n 40602 fes_blocks_2d_HB_SC.dat > temp_fes_2d_c.dat
tail -n 201 fes_HB.dat > temp_fes_HB_c.dat
tail -n 201 fes_SC.dat > temp_fes_SC_c.dat
cp temp* ../../OneOPES_b30/.
#
cd ../../OneOPES_b30_d/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv 5 --bin 200 --temp 320 --skip 20000 --out fes_HB.dat --min -63 --max 7 &
./../../FES_from_Reweighting.py --sigma 0.25 --colvar COLVAR.* --cv 6 --bin 200 --temp 320 --skip 20000 --out fes_SC.dat --min -15 --max 20 &
./../../FES_from_Reweighting.py --sigma 0.5,0.25 --colvar COLVAR.* --cv 5,6 --bin 200,200 --skip 20000 --temp 320 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-63, -15" --max 7,20;
tail -n 40602 fes_blocks_2d_HB_SC.dat > temp_fes_2d_d.dat
tail -n 201 fes_HB.dat > temp_fes_HB_d.dat
tail -n 201 fes_SC.dat > temp_fes_SC_d.dat
cp temp* ../../OneOPES_b30/.
#
cd ../../OneOPES_b30_e/0
rm temp*
./../../FES_from_Reweighting.py --sigma 0.5 --colvar COLVAR.* --cv 5 --bin 200 --temp 320 --skip 20000 --out fes_HB.dat --min -63 --max 7 &
./../../FES_from_Reweighting.py --sigma 0.25 --colvar COLVAR.* --cv 6 --bin 200 --temp 320 --skip 20000 --out fes_SC.dat --min -15 --max 20 &
./../../FES_from_Reweighting.py --sigma 0.5,0.25 --colvar COLVAR.* --cv 5,6 --bin 200,200 --skip 20000 --temp 320 --blocks 3 --out fes_blocks_2d_HB_SC.dat --min "-63, -15" --max 7,20;
tail -n 40602 fes_blocks_2d_HB_SC.dat > temp_fes_2d_e.dat
tail -n 201 fes_HB.dat > temp_fes_HB_e.dat
tail -n 201 fes_SC.dat > temp_fes_SC_e.dat
cp temp* ../../OneOPES_b30/.
#


paste temp_fes_2d_a.dat temp_fes_2d_b.dat temp_fes_2d_c.dat temp_fes_2d_d.dat temp_fes_2d_e.dat | awk '{ if ($1 == "" || $2 == "" || $3 == "" || $4 == "") print ""; else print $1, $2, ($3+$7+$11+$15+$19)/5 }' > fes_2d_HB_SC_average.dat
paste temp_fes_HB_a.dat temp_fes_HB_b.dat temp_fes_HB_c.dat temp_fes_HB_d.dat temp_fes_HB_e.dat | awk '{ if ($1 == "" || $2 == "" || $3 == "") print ""; else print $1, ($2+$4+$6+$8+$10)/5 }' > fes_HB_average.dat
paste temp_fes_SC_a.dat temp_fes_SC_b.dat temp_fes_SC_c.dat temp_fes_SC_d.dat temp_fes_SC_e.dat | awk '{ if ($1 == "" || $2 == "" || $3 == "") print ""; else print $1, ($2+$4+$6+$8+$10)/5 }' > fes_SC_average.dat

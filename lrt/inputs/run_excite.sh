cp excite*.gro a.gro
bash simulation_lrt.sh
bash correction_lrt.sh
bash make_frames.sh
echo 15  47 48 | gmx_mpi energy -f md.edr -o md_T.xvg 
rm \#*
pdb_tidy mdframe.pdb | pdb_selchain -A,B | pdb_rplchain -A:H | pdb_rplchain -B:L > fv.pdb
python cdr_metrics.py fv.pdb 2FJG_fv_chothia.pdb 2FJG_fv_chothia.pdb

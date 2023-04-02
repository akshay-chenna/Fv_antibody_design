mkdir temp
for i in {01..91} 
do
	mpirun -np 1 $ROSETTA_BIN/identify_cdr_clusters.mpi.linuxgccrelease -input_ab_scheme Chothia -in:file:l x$i -in:path renumbered -scorefile clusters.sc -allow_omega_mismatches_for_north_clusters -out:path:pdb temp/ &
done
wait
rm -rf temp

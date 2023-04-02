mkdir mdframes/
echo 1 | gmx_mpi trjconv -f md_corrected.xtc -o mdframes/mdframe.pdb -s frame_cluster.tpr -n index.ndx -split 1 -b 101

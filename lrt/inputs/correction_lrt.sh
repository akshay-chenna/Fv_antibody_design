echo 1 0 | gmx_mpi trjconv -f a.gro -s md.tpr -o frame_cluster.gro -pbc cluster
gmx_mpi grompp -f md_lrt.mdp -c frame_cluster.gro -o frame_cluster.tpr -maxwarn 1
echo 0 | gmx_mpi trjconv -f md.xtc -o md_corrected.xtc -s frame_cluster.tpr -pbc nojump
echo 0 | gmx_mpi trjconv -f md_corrected.xtc -o mdframe.pdb -b 200  -s frame_cluster.tpr
echo 0 | gmx_mpi trjconv -f md.gro -o md_corrected.gro -s frame_cluster.tpr
echo 1 | gmx_mpi trjconv -f md_corrected.xtc -o mdframe_solute.pdb -b 200 -s frame_cluster.tpr

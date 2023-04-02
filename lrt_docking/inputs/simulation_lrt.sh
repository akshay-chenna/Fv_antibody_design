export OMP_NUM_THREADS=4
export CUDA_MPS_PIPE_DIRECTORY=$PWD
mkdir -p $CUDA_MPS_PIPE_DIRECTORY
nvidia-cuda-mps-control -d

gmx_mpi grompp -f md_lrt.mdp -c a.gro -r a.gro -p topol.top -o md.tpr -po mdout.mdp -maxwarn 1 
CUDA_VISIBLE_DEVICES=1 mpirun -np 1 gmx_mpi mdrun -v -s md.tpr -o md.trr -x md.xtc -cpo md.cpt -e md.edr -g md.log -c md.gro -ntomp 4 -nstlist 150 -nb gpu -bonded gpu -pme gpu -update gpu

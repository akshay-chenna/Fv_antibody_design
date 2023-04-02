rm inputs/*npy
mv *npy inputs/.
for i in {5..30..5}  
do 
	mkdir $i
	cp -r inputs/* $i/.  
	cp e*_${i}K.gro $i/.
	cd $i 
	nohup bash run_excite.sh &
	cd ..
done
wait
python rmsd_vlrt.py

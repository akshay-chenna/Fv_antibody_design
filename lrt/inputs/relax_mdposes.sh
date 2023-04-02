for i in {49..98}
do
	python relax_mdposes.py mdframes/mdframe${i}.pdb 1 &
done
for i in {99..99}
do
	python relax_mdposes.py mdframes/mdframe${i}.pdb 3 &
done
wait

i=1
while read -r l 
do 
	python active-passive.py uniq_fv_9k_renumbered/$l active-passive_${i}.txt 
	i=$((i+1)) 
done < list.txt

mkdir active-passive_ab
mv active-passive*txt active-passive_ab/ 

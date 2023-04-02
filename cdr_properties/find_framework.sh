i=1
while read -r l 
do 
	python find_framework.py uniq_fv_9k_renumbered/$l framework_${i}.txt 
	i=$((i+1)) 
done < list.txt

mkdir frameworks
mv framework*txt frameworks/ 

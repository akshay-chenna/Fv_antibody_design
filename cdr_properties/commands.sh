ln -s ~/datasets/oas_paired95_flat-igfold .
cp oas_paired95_flat-igfold/oas_igfold.fasta . 
grep ">" oas_igfold.fasta  | sed 's/^>//g' | sed 's/$/.pdb/g' >> list.txt

mkdir renumbered
python anarci_.py # Obtained from igfold/utils. Modified to parallelize on 90 cores.
ls renumbered/ >> list.txt

#split -n 90 --numeric-suffixes=1 list.txt
split -l 1166 --numeric-suffixes=1 list.txt

nohup bash run_cluster_identify.sh & # Identifies clusters
mkdir lists
mv x* list.txt lists/.
rm *log

cat 1000_names.txt | cut -d . -f 1 | while read -r l ; do grep $l clusters.sc >> 1000_clusters.sc ; done
# 1000_names are 1000 sequentially unique. Based on kmeans clustering on the sequence similarity score.



### Analysis
awk '{ print $7 }' clusters.sc | grep H3 | sort -u >> H3_clusters.txt
awk '{ print $7 }' 1000_clusters.sc | grep H3 | sort -u >> 1000_H3_clusters.txt
cat 1000_H3_clusters.txt H3_clusters.txt | sort | uniq -u >> uniq_not-in_1000.txt
shuf -n 10000 lists/list.txt > 10k_random_names.txt

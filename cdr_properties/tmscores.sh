ln -s ../datasets/oas_paired95_flat-igfold .
cp ../tests/seq_identity_calc/oas_igfold_names.txt .
sed -i 's/^>//' oas_igfold_names.txt
sed -i 's/$/.pdb/g' oas_igfold_names.txt
cp ../tests/seq_identity_calc/1000_names.txt .
sed -i 's/$/.pdb/g' 1000_names.txt
~/apps/tmalign/TMalign -ter 1 -outfmt 2 -dir oas_paired95_flat-igfold/ 1000_names.txt >> 1000_tmscores.txt & .

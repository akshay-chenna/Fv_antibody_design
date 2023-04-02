################## CLEAN PDB #####################
pdb_selchain -V 2FJG.pdb | pdb_reres -1 | pdb_keepcoord > chain_v.pdb
pdb_selchain -W 2FJG.pdb | pdb_reres -1 | pdb_keepcoord > chain_w.pdb
pdb_selchain -L,H 2FJG.pdb | pdb_keepcoord > chain_LH.pdb
pdb_merge chain_v.pdb chain_w.pdb chain_LH.pdb | pdb_delhetatm | pdb_tidy >> 2FJG_modified.pdb
################## FIND residues in contact with antibody ###############
get_static_contacts.py --structure 2FJG_modified.pdb --output contacts.tsv --itypes all --sele "chain H L"  --sele2 "chain V W"
get_contact_frequencies.py --input_files contacts.tsv --output_file freq.tsv
awk '{print $2 }' freq.tsv | sort -u  >> antigen_contacts.txt   # Might have to edit antigen_contacts.txt file
awk '{print $1 }' freq.tsv | sort -u  >> antibody_contacts.txt
################# Convert HL  ###################################
pdb_selchain -H chain_LH.pdb | pdb_tofasta | tail -n +2 | tr -d "\n" > H_fasta.txt
pdb_selchain -L chain_LH.pdb | pdb_tofasta | tail -n +2 | tr -d "\n" > L_fasta.txt
python get_fv_lengths.py >> lengths.txt
################ Extract Fvs of H and L ########################
pdb_selchain -L chain_LH.pdb | pdb_selres -:`head -1 lengths.txt` > L_fv.pdb
pdb_selchain -H chain_LH.pdb | pdb_selres -:`tail -1 lengths.txt` > H_fv.pdb
pdb_merge H_fv.pdb L_fv.pdb | pdb_delhetatm | pdb_tidy > 2FJG_fv.pdb
######## Renumber pdb to Chothia ################
python anarci_.py 2FJG_fv.pdb 2FJG_fv_chothia.pdb

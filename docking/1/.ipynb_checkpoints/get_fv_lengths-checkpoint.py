from abnumber import Chain
with open('L_fasta.txt','r') as f:
        x = f.readlines()
print(len(Chain(x[0], scheme='chothia')))

with open('H_fasta.txt','r') as f:
        x = f.readlines()
print(len(Chain(x[0], scheme='chothia')))

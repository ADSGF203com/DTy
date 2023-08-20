# Get allelic and HierCC callings
python DTy.py -q examples/ERR1129648.fasta -o examples/ERR1129648.alleles

head -4 examples/ERR1129648.alleles

# Generate tab-delimited allelic profiles
python alleles2profile.py -p examples/examples examples/*.alleles

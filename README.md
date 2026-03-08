# Malaria Case Study

---

## Genome Cleaning and Gene Prediction

First, remove small scaffolds and predict genes. 

```bash
# Remove small scaffolds from the raw genome
python removeScaffold.py Haemoproteus_tartakovskyi.raw.genome 40 clean.genome 3000

# Predict genes using GeneMark-ES
gmes_petap.pl --ES --sequence Haemoproteus_tartakovskyi.raw.genome --min_contig 5000 --cores 30 &

# Clean headers in the genome file
sed "s/ GC=.*//" clean.genome > clean2.genome

# Clean headers in the GTF file
cat genemark.gtf | sed "s/  length=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht.gtf

# Parse GFF to generate nucleotide (.fna) and protein (.faa) fasta files
perl gffParse.pl -i clean2.genome -g Ht.gtf -c -p -b Ht

# Organize outputs
mkdir Ht_result
mv Ht.faa Ht.fna Ht.gtf Ht_result
```

## Contamination Filtering via BLAST and Taxonomy Parsing


```bash
# Create a local SwissProt BLAST database
makeblastdb -in SwissProt.fasta -dbtype prot -out SwissProt.blastq

# Run BLASTP against the SwissProt database
nohup blastp -query Ht_result/Ht.faa -db SwissProt -out Ht.blastp -num_threads 20

# Link and parse taxonomy data to identify target/contaminant scaffolds
ln -s /resources/binp29/Data/malaria/taxonomy.dat taxonomy.dat
cp ../../../resources/binp29/Data/malaria/datParser.py .
python datParser.py Ht.blastp Ht_result/Ht.fna taxonomy.dat uniprot_sprot.dat > scaffolds.txt

# Filter the genome and proteome based on the parsed scaffolds
# (Requires custom script: clean_genome.py)
python clean_genome.py Ht_result/Ht.fna scaffolds.txt Ht_result/Ht_clean.fna
python clean_genome.py Ht_result/Ht.faa scaffolds.txt Ht_result/Ht_clean.faa

# Re-run GeneMark on the cleaned genome
nohup gmes_petap.pl --ES --sequence Ht_clean.fna --cores 30 
```

## Orthology Inference

```bash
# Parse GFF files for all comparison species
perl gffParse.pl -i Plasmodium_berghei.genome -g PutGenomeHere/P_berghei.gtf -c -p -b Pb
perl gffParse.pl -i Plasmodium_yoelii.genome -g PutGenomeHere/Plasmodium_yoelii.gtf -c -p -b Py
perl gffParse.pl -i Plasmodium_cynomolgi.genome -g PutGenomeHere/cynomolgi.gtf -c -p -b Pc
perl gffParse.pl -i Plasmodium_knowlesi.genome -g PutGenomeHere/knowlesi.gtf -c -p -b Pk
perl gffParse.pl -i Plasmodium_faciparum.genome -g PutGenomeHere/Pfalciparum.gtf -c -p -b Pf
perl gffParse.pl -i Toxoplasma_gondii.genome -g PutGenomeHere/Tg.gff -c -p -b Tg
perl gffParse.pl -i Plasmodium_vivax.genome -g PutGenomeHere/vivax.gtf -c -p -b Pv

# Clean FASTA files to remove invalid amino acid characters
sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' *.faa

# Run ProteinOrtho for orthology detection
nohup proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &

# Organize FASTA files
mkdir faa_all
mv *.faa faa_all/
cp *.faa faa_all/
```

## Phylogenomics and Tree Construction

```bash
# Extract single-copy orthologs present in all 8 species
mkdir ortho_results
mkdir faa
mv *.faa faa/
awk 'NR==1 || ($1 == 8 && $2 == 8)' myproject.proteinortho.tsv > core_8_8.tsv

# Align sequences using MAFFT
cd ortho_results
mkdir aligned_files

for file in *.fasta; do
    mafft --auto "$file" > "aligned_files/${file%.fasta}.aln"
done

# Install RAxML via Conda
cd aligned_files
conda install -c bioconda raxml

# Build Maximum-Likelihood trees using RAxML
for file in *.aln; do
    basename=${file%.aln}
    
    raxmlHPC -s "$file" -n "${basename}.tre" -m PROTGAMMABLOSUM62 -p 12345
done

# Generate a consensus tree using PHYLIP
mkdir tree_file
mv *.tre tree_file/
cd tree_file/

# Concatenate all best trees
cat RAxML_bestTree.*.tre > intree

# Run PHYLIP consense to build the final tree
phylip consense
```

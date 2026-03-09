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
```
## Write a python script to clean the genome


Usage: python clean_genome.py Ht.fna scaffolds.txt Ht_clean.fna

```python
import sys

def filter_fasta(fasta_file, black_list_file, output_file):
    # 1. Read IDs to be removed into a set
    with open(black_list_file, 'r') as f:
        # One ID per line, stripping whitespace and newlines
        remove_ids = {line.strip() for line in f if line.strip()}

    print(f"Processing {fasta_file}, preparing to remove {len(remove_ids)} host 
scaffolds...")

    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        keep = False
        removed_count = 0
        kept_count = 0

        for line in f_in:
            if line.startswith('>'):
                # Check if the entire header line contains any ID from the blacklist
                is_blacklisted = any(r_id in line for r_id in remove_ids)

                if is_blacklisted:
                    keep = False
                    removed_count += 1
                else:
                    keep = True
                    kept_count += 1
                    f_out.write(line)
            else:
                # If the current scaffold is not blacklisted, write the sequence lines
                if keep:
                    f_out.write(line)

    print(f"Sequences kept: {kept_count}")
    print(f"Sequences removed: {removed_count}")
    print(f"Result saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python clean_genome.py <input_fasta> <blacklist_file> <output_fasta>")
    else:
        filter_fasta(sys.argv[1], sys.argv[2], sys.argv[3])

```
Than run the python script to clean genome.

```bash
# Filter the genome and proteome based on the parsed scaffolds
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

# Remove invalid amino acid characters
sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' *.faa

# Run Proteinortho for orthology detection
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

# Align sequences
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

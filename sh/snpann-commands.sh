# requires correctly formatted GFF. Use script `gbk2gff-bcft-csq.py` to convert a standard genbank file into this format

#bcftools csq -f ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.fna -g ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.csq.gff Pf_ANC.ASM922.bwa-mem.octopus.vcf -Ov -o sample.csq.vcf

## snpEff has a very peculiar requirement for building its reference database

# fasta needs suffix .fa not .fna or .fasta

# ref can be whatever name you choose. Also must change this in the config.file and need to change the chromosome name ref.MYCHROMNAME.codonTable : 

ref=ASM922v1
data/${ref}/genes.gff
data/genomes/${ref}.fa

# snpeff_data/
# ├── ASM922v1
# │   ├── genes.gff
# │   └── snpEffectPredictor.bin
# ├── genomes
# │   └── ASM922v1.fa
# └── snpeff.config

cd snpeff_data

snpEff build -c snpeff.config -dataDir . -gff3 ASM922v1

snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c snpeff.config -dataDir . ASM922v1 ../Pf_E3.ASM922.bwa-mem.octopus.vcf > ../Pf_E3.ASM922.bwa-mem.octopus.snpeff.vcf

cd ..

# make snippy formatted tsv
snippy-vcf_to_tab --ref ../PsFluSBW25/PsFluSBW25.fna --gff ../PsFluSBW25/PsFluSBW25.gff --vcf Pf_E2.ASM922.bwa-mem.octopus.snpeff.vcf > Pf_E2.ASM922.final.tsv

# get other stats from vcf file
bcftools query -f '%CHROM %FILTER %POS %AN %REF %ALT "[%PS]" "[%GT]" "[%ADP]" "[%AF]" "[%MAP_HF]"\n' Pf_E1.ASM922.bwa-mem.octopus.vcf > Pf_E1.ASM922.hapfreq.tsv

# snippy can work with bcftools csq but it doesn't parse/report information about the intragenic regions. Snippy works better with snpEff generated annotations

snippy-vcf_to_tab --ref ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.fna --gff ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.gff --vcf sample.csq.vcf > bcftools_csq.tsv

# works better with snpEff generated annotations
snippy-vcf_to_tab --ref ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.fna --gff ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.gff --vcf snpeff_ann.vcf > snpeff_ann.tsv
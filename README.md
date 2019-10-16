kgp2anc
=======

A protocol to estimate ancestry from the main five ethnicities (European, African, East Asian, Native American, South Asian) using principal components analysis together with individuals from the 1000 Genomes project phase 3 starting from raw Illumina genotype array data. For any feedback, send an email to giulio.genovese@gmail.com

[![](https://img.youtube.com/vi/XkCCE8jC1F4/0.jpg)](https://www.youtube.com/watch?v=XkCCE8jC1F4)

Installation
============

Install basic tools (Debian/Ubuntu specific):
```
sudo apt install wget gzip unzip samtools python3-pandas r-base-core
```

Optionally, you can install these libraries to activate further HTSlib features:
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res/kgp && cd /tmp
```

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Compile latest version of HTSlib (optionally disable bz2, gcs, and lzma) and BCFtools (make sure you are using gcc version 5 or newer)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/bcftools $HOME/bin/
```

Make sure the directory with the plugins is available to bcftools
```
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

Install latest development version of plink 1.9
```
wget http://s3.amazonaws.com/plink1-assets/dev/plink_linux_x86_64.zip
unzip -od $HOME/bin/ plink_linux_x86_64_dev.zip plink
```

Install latest version of GCTA
```
wget http://cnsgenomics.com/software/gcta/gcta_1.91.6beta.zip
unzip -ojd $HOME/bin/ gcta_1.91.6beta.zip gcta_1.91.6beta/gcta64
```

Install scripts
```
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/kgp2anc/master/{vcf2plink.py,markerqc.sh,kgpmerge.sh,kgp2pc.py,pc2anc.R}
chmod a+x $HOME/bin/{vcf2plink.py,markerqc.sh,kgpmerge.sh,kgp2pc.py,pc2anc.R}
```

Lift over list of <a href="http://doi.org/10.1016/j.ajhg.2008.06.005">long-range LD</a> regions to be excluded from principal component analysis from reference build NCBI 36 to reference build GRCh37 and GRCh38
```
cd $HOME/res
wget https://raw.githubusercontent.com/freeseek/kgp2anc/master/ld.hg18.bed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg{19,38}.over.chain.gz
ln -s hg18ToHg19.over.chain.gz hg18ToHg37.over.chain.gz
for build in 37 38; do
  ./liftOver \
    -minMatch=.4 \
    ld.hg18.bed \
    hg18ToHg$build.over.chain.gz \
    /dev/stdout \
    /dev/stderr | \
    sed 's/^chr//' > ld.grch$build.bed
done
/bin/rm hg18ToHg37.over.chain.gz
```

Download 1000 Genomes project phase 3 population assignments
```
cd $HOME/res/kgp
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
(echo -e "IID\tPOP"; cut -f2,7 20130606_g1k.ped | tail -n+2 | sort) > kgp.pop
```

Download global ancestry estimates from 1000 Genomes project phase 3
```
cd $HOME/res/kgp
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140818_ancestry_deconvolution/{ACB,ASW,CLM,MXL,PEL,PUR}_phase3_ancestry_deconvolution.zip
for pop in ACB ASW CLM MXL PEL PUR; do
  unzip -jo ${pop}_phase3_ancestry_deconvolution.zip $pop/PopPhased/lai_global_$pop.txt
done
(echo -e "IID\tAFR\tEUR\tNAT\tUNK"; tail -qn+2 lai_global_{ACB,ASW,CLM,MXL,PEL,PUR}.txt) > kgp.anc
```

Download list of Illumina Omni 2.5 and Affymetrix 6.0 markers
```
cd $HOME/res/kgp
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/Omni25_genotypes_2141_samples.b37.v2.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/coriell_affy6_intensities/ALL.wgs.nhgri_coriell_affy_6.20131213.snps_indels.chip.genotypes.vcf.gz
bcftools query -f "%CHROM\t%POS\n" Omni25_genotypes_2141_samples.b37.v2.vcf.gz | sort | uniq > kgp.omni25.grch37.chr.pos
bcftools query -f "%CHROM\t%POS\n" ALL.wgs.nhgri_coriell_affy_6.20131213.snps_indels.chip.genotypes.vcf.gz | sort | uniq > kgp.affy60.grch37.chr.pos
join -t@ kgp.omni25.grch37.chr.pos kgp.affy60.grch37.chr.pos > kgp.array.grch37.chr.pos
```

Download GRCh37 human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/res/human_g1k_v37.fasta
samtools faidx $HOME/res/human_g1k_v37.fasta
```

Download 1000 Genomes project phase 3 genotypes for GRCh37 (notice this loop takes several hours to run)
```
cd $HOME/res/kgp
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{{1..22}.phase3_shapeit2_mvncall_integrated_v5a,X.phase3_shapeit2_mvncall_integrated_v1b,Y.phase3_integrated_v2a}.20130502.genotypes.vcf.gz{,.tbi}
for chr in {1..22} X Y; do
  bcftools view --no-version -Ou -c 2 ALL.chr${chr}.phase3*integrated_v[125][ab].20130502.genotypes.vcf.gz | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -o ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf -d none -f $HOME/res/human_g1k_v37.fasta && \
  bcftools index -f ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf
done
```

Subset phase 3 GRCh37 genotypes to array markers (takes a few hours)
```
cd $HOME/res/kgp
bcftools concat --no-version -Ov ALL.chr{{1..22},X}.phase3_integrated.20130502.genotypes.bcf | \
  awk 'NR==FNR {x[$1":"$2]++} NR>FNR && ($0~"^#" || $1":"$2 in x)' kgp.array.grch37.chr.pos - | \
  bcftools view --no-version -Ob -o kgp.array.grch37.bcf && \
  bcftools index -f kgp.array.grch37.bcf
```

Convert GRCh37 genotypes to plink format
```
cd $HOME/res/kgp
$HOME/bin/vcf2plink.py \
  --vcf kgp.array.grch37.bcf \
  --ref $HOME/res/human_g1k_v37.fasta \
  --build b37 \
  --impute-sex .6 .6 \
  --out kgp.array.grch37
```

Download GRCh38 human genome reference
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Download 1000 Genomes project phase 3 genotypes for GRCh38 (notice this loop takes several hours to run)
```
cd $HOME/res/kgp
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X,Y}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
for chr in {1..22} X Y; do
  (bcftools view --no-version -h ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
    grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
  bcftools view --no-version -H -c 2 ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
  grep -v "[0-9]|\.\|\.|[0-9]" | sed 's/^/chr/') | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -o ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -d none -f $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && \
  bcftools index -f ALL.chr${chr}_GRCh38.genotypes.20170504.bcf
done
```

Subset phase 3 GRCh38 genotypes to array markers (takes a few hours)
```
cd $HOME/res/kgp
bcftools concat --no-version -Ov ALL.chr{{1..22},X}_GRCh38.genotypes.20170504.bcf | \
  awk 'NR==FNR {x["chr"$1":"$2]++} NR>FNR && $0!~"^#" {y=$8; gsub(".*GRCH37_POS=","",y); gsub(";.*","",y)}
  NR>FNR && ($0~"^#" || $1":"y in x)' kgp.array.grch37.chr.pos - | \
  bcftools view --no-version -Ob -o kgp.array.grch38.bcf && \
  bcftools index -f kgp.array.grch38.bcf
```

Convert GRCh38 genotypes to plink format
```
cd $HOME/res/kgp
$HOME/bin/vcf2plink.py \
  --vcf kgp.array.grch38.bcf \
  --ref $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --build b38 \
  --impute-sex .6 .6 \
  --out kgp.array.grch38
```

Compute ancestry estimates
==========================

Setup parameters
```
build="38" # or build="37"
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
vcf="..." # the VCF file containing your sample's genotypes
pfx="..." # a prefix string for all the temporary files that will be generated
```

To convert Illumina or Affymetrix data to VCF format, use the <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> tools.

Convert VCF file to plink format
```
$HOME/bin/vcf2plink.py --build b$build --ref $ref --out $pfx --impute-sex .6 .6 --vcf $vcf
```

Perform quality control
```
$HOME/bin/markerqc.sh $pfx out/$pfx $HOME/res/ld.grch$build.bed
```

Merge cohort with 1000 Genomes project
```
$HOME/bin/kgpmerge.sh $pfx kgp/$pfx $HOME/res/kgp/kgp.array.grch$build out/$pfx.prune.in
```

Compute principal components
```
$HOME/bin/kgp2pc.py --grm-bin kgp/$pfx --fam kgp/$pfx.fam --out kgp/$pfx --pop $HOME/res/kgp/kgp.pop
```

Run a linear regression model to estimate main five ethnicities
```
$HOME/bin/pc2anc.R kgp/$pfx.all.pca $HOME/res/kgp/kgp.anc kgp/$pfx.anc.tsv
```

Adjust ancestry estimates to force values within closed interval [0,1]
```
awk -F"\t" -v OFS="\t" 'NR>1 {for (i=4; i<=NF; i++) {if ($i<0) $i="0.0000"; if ($i>1) $i="1.0000"}} {print}' kgp/$pfx.anc.tsv | tr -d " " > kgp/$pfx.adjanc.tsv
```

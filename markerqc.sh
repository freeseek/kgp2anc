#!/bin/bash
###
#  The MIT License
#
#  Copyright (c) 2016-2018 Giulio Genovese
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

# performs basic quality control of markers
# requires plink prefix of the input dataset
# requires plink prefix of the output dataset
# optional file with list of samples to be removed
# optional minimum frequency of missingness allowed (I use .02 for Illumina, and .05 for Affymetrix)

set -e -o pipefail

if [ $# -lt 3 ]; then
  echo "About:   Perform marker quality control using plink (Jul 26th 2018)"
  echo "Usage:   markerqc.sh <input prefix> <output prefix> <long range regions> [remove range] [lmiss threshold]"
  exit 1
fi

in="$1"
out="$2"
ld="$3"

if [ ! -f $in.bed ] || [ ! -f $in.bim ] || [ ! -f $in.fam ]; then
  printf "Input file %s.bed or %s.bim or %s.fam not present\n" $in $in $in >&2
  exit 1
fi

if [ ! -f $ld ] || [ ! -s $ld ]; then
  printf "File range regions file %s not present or empty\n" $ld >&2
  exit 1
fi

# makes sure the ouput directory exists
dir=$(dirname $out)
mkdir -p $dir

if [ $# -ge 4 ]; then
  if [ ! -f $4 ] || [ ! -s $4 ]; then
    printf "File %s not present or empty\n" $4 >&2
    exit 1
  fi
  opt="--remove $4"
fi

if [ $# -ge 5 ]; then
  lmiss="$5"
else
  lmiss=".02"
fi

# identifies markers not on autosomes
awk '$1=="0" || $4=="0" {print $2}' $in.bim | sort > $out.noref

# identifies markers with excess missingness
plink --bfile $in --keep-allele-order $opt --freq counts gz --out $out
gzip -cd $out.frq.counts.gz | awk 'NR>1 && ($5=="NA" || $5==0 || $6==0) {print $2}' | sort > $out.nofrq
gzip -cd $out.frq.counts.gz | awk 'NR>1 && $5!="NA" && $5!=0 && $5!=1 && ($5/($5+$6)<.005 || $5/($5+$6)>.995) {print $2}' | sort > $out.frq.ex

# identifies markers with too elevated levels of missingness
plink --bfile $in --keep-allele-order $opt --exclude $out.nofrq --missing gz --out $out
gzip -cd $out.lmiss.gz | awk -v f=$lmiss 'NR>1 && $5>=f {print $2}' | sort > $out.miss.ex

# generates historgrams of missingness for individuals and for markers
gzip -cd $out.lmiss.gz | awk 'NR>1 {x[$3]++} END {for (i in x) print i"\t"x[i]}' | sort -k1,1n | awk '{SUM+=$2; print $0"\t"SUM}' > $out.lhist
gzip -cd $out.imiss.gz | awk 'NR>1 {x[$4]++} END {for (i in x) print i"\t"x[i]}' | sort -k1,1n | awk '{SUM+=$2; print $0"\t"SUM}' > $out.ihist

# identifies markers failing Hardy Weinberg equilibrium test due to excess heterozygousness
plink --bfile $in --keep-allele-order $opt --hardy gz --out $out
gzip -cd $out.hwe.gz | awk 'NR>1 && $3~"ALL" && ($9<1e-6 && $7>$8) {print $2}' | sort > $out.hwe.ex

# identifies autosomal markers that associate with sex
plink --bfile $in --autosome-xy --pheno $in.fam --mpheno 3 --assoc fisher --out $out.sex
awk 'NR>1 && $8<1e-6 {print $2}' $out.sex.assoc.fisher | sort > $out.sex.ex
gzip -f $out.sex.assoc.fisher

# generate list of markers falling within long range regions
plink --bfile $in --make-set $ld --write-set --out $out
grep -v "^LD$\|^END$\|^$" $out.set | sort > $out.ld.ex

# generates file with list of markers to be removed
if [ -f $out.nopass ]; then
  cat $out.{nopass,noref,nofrq,{frq,miss,hwe,sex}.ex} | sort | uniq > $out.ex
else
  cat $out.{noref,nofrq,{frq,miss,hwe,sex}.ex} | sort | uniq> $out.ex
fi
join -a1 -a2 $out.ex $out.ld.ex > $out.pca.ex

# compute set of independent markers
plink --bfile $in --keep-allele-order --exclude $out.pca.ex --maf .01 --indep 50 5 2 --out $out

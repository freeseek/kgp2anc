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

# merges a plink dataset with KGP samples
# requires plink prefix of the input dataset
# requires plink prefix of the output dataset
# optional file with list of markers to be extracted (usually a prune.in file)
# optional file with list of samples to be removed

set -e -o pipefail

if [ $# -lt 3 ]; then
  echo "About:   Merge plink dataset with 1000 Genomes project phase 3 dataset (Jul 26th 2018)"
  echo "Usage:   kgpmerge.sh <input prefix> <output prefix> <kgp prefix> [extract range] [remove range]"
  exit 1
fi

in="$1"
out="$2"
kgp="$3"

if [ ! -f $in.bed ] || [ ! -f $in.bim ] || [ ! -f $in.fam ]; then
  printf "Input file %s.bed or %s.bim or %s.fam not present\n" $in $in $in >&2
  exit 1
fi

if [ ! -f $kgp.bed ] || [ ! -f $kgp.bim ] || [ ! -f $kgp.fam ]; then
  printf "Reference panel files %s.bed or %s.bim or %s.fam not present\n" $kgp $kgp $kgp >&2
  exit 1
fi

# makes sure the ouput directory exists
dir=$(dirname $out)
mkdir -p $dir

# if list of markers with at least 1% frequency in KGP set does not exist, create it
if [ ! -f $kgp.maf.001 ] || [ ! -s $kgp.maf.001 ]; then
  if [ ! -f $kgp.frq.counts.gz ] || [ ! -s $kgp.frq.counts.gz ]; then
    plink --bfile $kgp --keep-allele-order --freq counts gz --out $kgp
  fi
  gzip -cd $kgp.frq.counts.gz | tail -n+2 | awk '$5/($5+$6)>.01 || $6/($5+$6)>.01 {print $2}' | sort > $kgp.maf.001
fi

# generate a list of markers to extract from each dataset
if [ $# -ge 4 ]; then
  if [ ! -f $4 ] || [ ! -s $4 ]; then
    printf "File %s not present or empty\n" $4 >&2
    exit 1
  fi
  cat $4
else
  awk '{print $2}' $in.bim
fi | sort | join - $kgp.maf.001 > $out.extract

if [ $# -ge 5 ]; then
  if [ ! -f $5 ] || [ ! -s $5 ]; then
    printf "File %s not present or empty\n" $5 >&2
    exit 1
  fi
  opt="--remove $5"
fi

# extract markers from the provided dataset and the KGP dataset and merge the two datasets
plink --bfile $in --keep-allele-order --extract $out.extract $opt --make-bed --out $out.set
plink --bfile $kgp --keep-allele-order --extract $out.extract --make-bed --out $out.kgp
plink --bfile $out.set --keep-allele-order --bmerge $out.kgp --make-bed --out $out

# cleanup
rm $out.extract $out.set.bed $out.set.bim $out.set.fam $out.set.log $out.kgp.bed $out.kgp.bim $out.kgp.fam $out.kgp.log

# compute similarity matrix
plink --bfile $out --keep-allele-order --autosome --make-grm-bin --out $out --thread-num 10

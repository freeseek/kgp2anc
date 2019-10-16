#!/usr/bin/env python3
"""
   The MIT License

   Copyright (c) 2016-2018 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
"""

import argparse, os, sys
from subprocess import call, Popen, PIPE

parser = argparse.ArgumentParser(description = 'vcf2plink.py: converts a VCF file to plink format (Jul 26th 2018)', add_help = False, usage = 'vcf2plink.py --ref <fasta> --build <b37/b38> [options]')
parser.add_argument('--vcf', metavar = '<in.vcf.gz>', type = str, default = '-', help = 'Specify a VCF file to be converted [stdin]')
parser.add_argument('--ref', metavar = '<file>', type = str, required = True, help = 'reference sequence (e.g. human_g1k_v37.fasta)')
parser.add_argument('--build', metavar = '<str>', type = str, required = True, nargs = '*', help = 'reference build code (see "plink --help --split-x")')
parser.add_argument('--filter', metavar = '<expr>', type = str, help = 'exclude sites for which the expression is true (e.g. "FORMAT/DP<10 || FORMAT/GQ<20") (see bcftools manual for details)')
parser.add_argument('--set-GTs', metavar = '<char>', type = str, default = '.', help = 'set genotypes of failed samples to missing (.) or ref (0) (see "bcftools filter -h") [.]')
parser.add_argument('--out', metavar = '[prefix]', type = str, default = 'plink', help = 'Specify prefix for output files (see "plink --help --out") [plink]')
parser.add_argument('--mem', metavar = '<int>', type = int, help = 'main workspace size, in GB')
parser.add_argument('--impute-sex', metavar = '{female max F} {male min F}', type = float, nargs = '*', help = 'Impute sex from chromosome X (see "plink --help --impute-sex")')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

# make sure the ouput directory exists
outdir = os.path.dirname(args.out)
if outdir != '' and not os.path.isdir(outdir):
  os.makedirs(outdir)

# invoke bcftools to preprocess VCF file
if args.vcf == '-':
  cmd = ['bcftools', 'norm', '-Ou', '-m', '-any']
  sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + ' |\n')
  p1 = Popen(cmd, stdin = sys.stdin.buffer, stdout = PIPE, stderr = sys.stderr.buffer)
else:
  cmd = ['bcftools', 'norm', '-Ou', '-m', '-any', args.vcf]
  sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + ' |\n')
  p1 = Popen(cmd, stdout = PIPE, stderr = sys.stderr.buffer)

cmd = ['bcftools', 'norm', '-Ou', '-f', args.ref]
sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + ' |\n')
p2 = Popen(cmd, stdin = p1.stdout, stdout = PIPE, stderr = sys.stderr.buffer)

if args.filter:
  cmd = ['bcftools', 'filter', '-Ou', '-e', args.filter, '--set-GTs', args.set_GTs]
  sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + ' |\n')
  p2a = Popen(cmd, stdin = p2.stdout, stdout = PIPE, stderr = sys.stderr.buffer)

cmd = ['bcftools', 'annotate', '-Ob', '-x', 'ID', '-I', '+%CHROM:%POS:%REF:%ALT']
sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + ' |\n')
p3 = Popen(cmd, stdin = p2a.stdout if args.filter else p2.stdout, stdout = PIPE, stderr = sys.stderr.buffer)

# invoke plink to convert VCF file
cmd = ['plink', '--make-bed',
  '--keep-allele-order',
  '--bcf', '/dev/stdin',
  '--vcf-idspace-to', '_',
  '--const-fid',
  '--allow-extra-chr', '0',
  '--split-x'] + args.build + ['no-fail',
  '--out', args.out] + (['--memory', str(1024 * args.mem)] if args.mem else [])
sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + '\n')
if call(cmd, stdin = p3.stdout, stdout = sys.stdout.buffer, stderr = sys.stderr.buffer):
    raise Exception('vcf2plink.py: Error while converting VCF file to plink')

# impute sex using chromosome X heterozygosity
if args.impute_sex:
  cmd = ['plink', '--make-bed',
  '--bfile', args.out,
  '--keep-allele-order',
  '--impute-sex', 'ycount', str(args.impute_sex[0]), str(args.impute_sex[1]), '10000', '0',
  '--out', args.out]
  sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + '\n')
  if call(cmd, stdout = sys.stdout.buffer, stderr = sys.stderr.buffer):
    raise Exception('vcf2plink.py: Error while imputing sex')
  cmd = ['rm', args.out + '.bed~', args.out + '.bim~', args.out + '.fam~']
  sys.stderr.write('vcf2plink.py: ' + ' '.join(cmd) + '\n')
  if call(cmd, stdout = sys.stdout.buffer, stderr = sys.stderr.buffer):
    raise Exception('vcf2plink.py: Error while removing temporary files')

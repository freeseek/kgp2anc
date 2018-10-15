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

import argparse, os, sys, pandas as pd, numpy as np
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description = 'kgp2pc.py: compute PC components from merged KGP datasets (Jul 26th 2018)', add_help = False, usage = 'kgp2pc.py [options]')
parser.add_argument("--grm-bin", metavar = "[prefix]", type = str, help = "Prefix for binary GRM file")
parser.add_argument("--pop", metavar = "[filename]", type = str, required = True, help = "File with population information")
parser.add_argument('--fam', metavar = '[filename]', type = str, required = True, help = 'Specify full name of .fam file')
parser.add_argument("--pca", metavar = "[int]", type = int, default = 20, help = "number of PCs")
parser.add_argument("--out", metavar = "[prefix]", type = str, default = "plink", help = "Specify prefix for output files")
parser.add_argument("--groups", metavar = "[groups]", type = str, default = 'ALL,AFAM,EUR', help = "Groups to be used for the principal component computations")
parser.add_argument("--remove", metavar = "[filename]", type = str, help = "Exclude all samples named in the file")
parser.add_argument('--xlsx', action = 'store_true', default = False, help = 'Whether the output table is an xlsx file.')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

if not os.path.isfile(args.grm_bin + '.grm.id') or not os.path.isfile(args.grm_bin + '.grm.bin'):
  raise Exception('kgp2pc.py: Input file ' + args.grm_bin + '.grm.id or ' + args.grm_bin + '.grm.bin not present')

if args.remove and not os.path.isfile(args.remove):
  raise Exception('kgp2pc.py: Input file ' + args.remove + ' not present')

names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHE']
dtypes = {'FID': str, 'IID': str, 'PAT': str, 'MAT': str, 'SEX': int, 'PHE': int}
try:
  df_fam = pd.read_csv(args.fam, delim_whitespace = True, header = None, names = names, dtype = dtypes)
except:
  raise Exception('kgp2pc.py: Failed to open ' + args.fam)
dtypes = {'FID': str, 'IID': str, 'POP': str}
try:
  df_pop = pd.read_csv(args.pop, delim_whitespace = True, dtype = dtypes)
except:
  raise Exception('kgp2pc.py: Failed to open ' + args.pop)
df = df_fam[['FID', 'IID', 'SEX']].merge(df_pop, how = 'left')
df.loc[df['POP'].isnull(), 'POP'] = 'SET'

pop = dict()
pop['EAS'] = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV']
pop['EUR'] = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']
pop['AFR'] = ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB']
pop['AMR'] = ['MXL', 'PUR', 'CLM', 'PEL']
pop['SAS'] = ['GIH', 'PJL', 'BEB', 'STU', 'ITU']
for lbl in pop:
  pop[lbl] = pop[lbl] + [lbl + '-' + x for x in pop[lbl]]
pop['SET'] = ['SET']
pop['ALL'] = pop['EAS'] + pop['EUR'] + pop['AFR'] + pop['AMR'] + pop['SAS']
pop['AFAM'] = pop['EUR'] + pop['AFR']

idx = dict()
for lbl in pop:
  idx[lbl] = np.array([x in pop[lbl] for x in df['POP']])

# invoke gcta64 to compute the principal components
for group in args.groups.split(','):
  if not group in idx:
    sys.stderr.write('kgp2pc.py: Group ' + group + ' undefined\n')
    continue
  sys.stderr.write('kgp2pc.py: Computing PC for group ' + group + '\n')
  cmd = ["gcta64",
    '--grm-bin', args.grm_bin,
    '--keep', '/dev/stdin',
    '--pca', str(args.pca),
    '--out', args.out + '.' + group.lower()]
  if args.remove:
    cmd += ['--remove', args.remove]

  sys.stderr.write('kgp2pc.py: ' + ' '.join(cmd) + '\n')
  p = Popen(cmd, stdin = PIPE, stdout = sys.stdout.buffer, stderr = sys.stderr.buffer, universal_newlines = True)
  df[np.logical_or(idx[group], idx['SET'])][['FID', 'IID']].to_csv(p.stdin, sep = '\t', index = False)
  p.communicate()
  if p.returncode:
    raise Exception('kgp2pc.py: Problems running gcta64')

  # add principal components information
  try:
    eigenval = pd.read_csv(args.out + '.' + group.lower() + '.eigenval', header = None, squeeze = True)
  except:
    raise Exception('kgp2pc.py: Failed to open ' + args.out + '.' + group.lower() + '.eigenval')
  names = ['FID', 'IID'] + ['PC' + str(i + 1) for i in range(args.pca)]
  dtypes = {'FID': str, 'IID': str}
  dtypes.update({'PC' + str(i + 1): float for i in range(args.pca)})
  try:
    eigenvec = pd.read_csv(args.out + '.' + group.lower() + '.eigenvec', delim_whitespace = True, header = None, names = names, dtype = dtypes)
  except:
    raise Exception('kgp2pc.py: Failed to open ' + args.out + '.' + group.lower() + '.eigenvec')
  for i in range(args.pca):
    eigenvec['PC' + str(i + 1)] *= eigenval[i]
  df_pca = df[np.logical_or(idx[group], idx['SET'])].merge(eigenvec, how = 'left')

  # generate output table
  if (args.xlsx):
    writer = pd.ExcelWriter(args.out + '.' + group.lower() + '.xlsx', engine = 'xlsxwriter')
    df_pca.to_excel(writer, sheet_name = 'Sheet1', index = False)
    writer.save()
  else:
    df_pca.to_csv(args.out + '.' + group.lower() + ".pca", sep = "\t", index = False)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:19:20 2020
@author: mlangleib
"""
import threading
import subprocess
import argparse
import glob
import os
from Bio import SeqIO

# computing PAML alignments with PAL2NAL
def pal2nal(seq_index):
  # define output name 
  paml_output = open('{0}.paml'.format(seq_index), 'w')
  # define command
  command = ['pal2nal.pl',
             '{0}.msa.aln'.format(seq_index),
             '{0}.seqs_nuc'.format(seq_index),
             '-output', 'paml', '-blockonly']
  # run sunprocess
  subprocess.run(command, stdout = paml_output)

# listing already computed ete3 calculations, in order to avoid re-computing them, as dict {<monophiletic_group_code>: <directory>}
codeml_dirs_computed = [{x[0].replace('.codeml_dir', ''): x[0]} for x in os.walk('../data/') if x[0].endswith('.codeml_dir')]

# listing monophyletic groups and running
monophyletic_groups = [msa_file.rpartition('.')[0] for msa_file in glob.glob('../data/*/' + '*.msa')]
# checking length //// makes sense
len(monophyletic_groups)

# creating clustal alignments
for group in monophyletic_groups:
  if os.path.isfile('{0}.msa.gbloMask_filtrada'.format(group)):
    SeqIO.convert('{0}.msa.gbloMask_filtrada'.format(group), 'fasta', '{0}.msa.aln'.format(group), 'clustal')
    subprocess.run("sed -i s/\\*//g {0}.msa.aln".format(group).split(' '))

for i in range(0, len(monophyletic_groups), 30):
    # defining step
    threads = []
    for index in range(0, 30, 1):
      if (i+index) < len(monophyletic_groups):
         t = threading.Thread(target=pal2nal, args=(monophyletic_groups[i+index],))
         threads.append(t)
         t.start()
         
    for t in threads:
         t.join()


# ordering paml files by number of seqs; filtering really big paml files
paml_files = glob.glob('../data/*/*paml')

# checking length of the list of generated files //// SOMETHING WEIRD, more than expected
len(paml_files) 

# performing codeml calculation on a few to try
monophyletic_groups = sorted(monophyletic_groups, key = lambda x: len(list(SeqIO.parse('{0}.msa'.format(x), 'fasta'))))

lista_prueba = monophyletic_groups[0:100]

def ete_calc(seq_index):
  # define log
  log_out = open('{0}.ete3_out'.format(seq_index), 'w') 
  # define command
  command = ['ete3',
             'evol',
             '-t', '{0}.msa.tree'.format(seq_index),
             '--alg', '{0}.paml'.format(seq_index),
             '--models', 'M0', 'M1', 'M2',
             '-o', '{0}.codeml_dir'.format(seq_index),
             '--clear_all',
             '--cpu', '1']
  subprocess.run(command, stdout = log_out)

for i in range(0, len(monophyletic_groups), 30):
    threads = []
    for index in range(0, 30, 1):
     if (i+index) < len(monophyletic_groups):
        t = threading.Thread(target=ete_calc, args=(monophyletic_groups[i+index],))
        threads.append(t)
        t.start()
        
    for t in threads:
         t.join()


# nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database /home/mauricio/programas_lab/databases/nr/nr --dbtype seqdb --qtype seq --target_orthologs all --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.diamondCOGs

# Con todos filos, todos los ortologos
#nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database /home/mauricio/programas_lab/databases/nr/nr --dbtype seqdb --qtype seq --target_orthologs all --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.diamondCOGs
# Con todos filos, ortologos unoauno
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --dbtype seqdb --qtype seq --target_orthologs one2one --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.all.one2one.diamondCOGs
# Con filo euk, todos
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --dbtype seqdb --qtype seq --target_orthologs all --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.euk.all.diamondCOGs
# Con filo euk, uno a uno
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --dbtype seqdb --qtype seq --target_orthologs one2one --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.euk.one2one.diamondCOGs
# Con filo animals, todos
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --tax_scope meNOG --dbtype seqdb --qtype seq --target_orthologs all --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.animals.all.diamondCOGs
# Con filo animals, uno a uno
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --tax_scope meNOG --dbtype seqdb --qtype seq --target_orthologs one2one --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.animals.one2one.diamondCOGs

# con filo bilateria, todos
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --tax_scope biNOG --dbtype seqdb --qtype seq --target_orthologs all --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.bilateria.all.diamondCOGs
# con filo bilateria, uno a uno
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --tax_scope biNOG --dbtype seqdb --qtype seq --target_orthologs one2one --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.bilateria.one2one.diamondCOGs
# Con filo nematodos, todos
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --tax_scope nemNOG --dbtype seqdb --qtype seq --target_orthologs all --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.nematodos.all.diamondCOGs
# Con filo nematodos, uno a uno
nohup python2.7 /home/mauricio/programas_lab/eggnog-mapper-0.99.3/emapper.py --database euk --tax_scope nemNOG --dbtype seqdb --qtype seq --target_orthologs one2one --go_evidence non-electronic -m diamond -i PROTEOMA_TOTAL.fasta -o PROTEOMA_TOTAL.nematodos.one2one.diamondCOGs

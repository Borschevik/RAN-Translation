import os
from Bio import SeqIO
path_to_fasta = str(input())
reserve_folder = str(input())
path_to_write_file = str(input())


def translation(mRNA):
    codon_dict = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                  'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L',
                  'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q',
                  'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
                  'AUG': 'M', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K',
                  'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
                  'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E',
                  'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GAG': 'E'}
    stop_list = ['UAA', 'UAG', 'UGA']
    str_rna = input()
    protein = str()
    l = len(str_rna)
    for i in range(0, l, 3):
        key = str_rna[int(i):int(i + 3)]
        if key in codon_dict:
            am_ac = codon_dict.get(key)
            protein += str(am_ac)
        else:
            break
    return protein

def read_fast(file_in, dir_out):
    for record in SeqIO.parse(open(file_in), "fasta"):
        f_out = os.path.join(dir_out, record.id + '.fasta')
        SeqIO.write([record], open(f_out, 'w'), "fasta")


def file_list(y):
    proteinfiels_f = []
    for d, dirs, files in os.walk(y):
        for f in files:
                proteinfiels = ''
                proteinfiels = os.path.join(d,f)
                proteinfiels_f.append(proteinfiels)
    return(proteinfiels_f)


read_fast(path_to_fasta,reserve_folder)
protein = file_list(reserve_folder)
#Analys of multifasta to find sequenses with amino acids repeats
name_list = list()
seq_l = list()
for el in protein:
    seq = str('')
    y = 0
    with open(el, 'r') as file:
        for line in file:
            if line.startswith('>'):
                name = line
            else:
                seq += str(line)
    list_rep = ('QAGR','MEWNG','LPAC','AAA','QQQ','SSS',
                'GGG' ,'PPP' ,'LLL','GPGP','GAGA','GRGR',
                'PRPR','GPGP')
    #These repeats take by analysis of known RAN-translation sites
    for el in list_rep:
        while y == 0:
            if el in seq:
                name_list.append(name)
                seq_l.append(seq)
                y = 1
            else:
                break
#Writing check multifasta
with open(path_to_write_file, 'w') as f:
    for i in range(0,len(name_list)):
        bb = str(name_list[i])
        cc = str(seq_l[i])
        f.write(bb)
        f.write(cc)

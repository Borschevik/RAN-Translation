from Bio import SeqIO

#input from fasta file
handle = open("Disease_db.fasta")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

print(len(records))


with open ('names_filtred_short_last_VERYLAST.fasta') as file:
    with open ('names_filtred_short_last_VERYLAST_RNA.fasta', 'w') as ouf:
        for record in SeqIO.parse(file, "fasta"):
            rec = record.id
            seq = record.seq.transcribe()
            ouf.write('>')
            ouf.write(str(rec))
            ouf.write('\n')
            #ouf.write('NNNNN')
            ouf.write(str(seq))
            ouf.write('\n')
    #seq.translate()
i=1
with open ('names_filtred_short_last_VERYLAST_RNA.fasta') as file:
    with open ('names_filtred_short_last_VERYLAST_Prot.fasta', 'w') as ouf:
        for record in SeqIO.parse(file, "fasta"):
            rec = record.id
            seq = record.seq.translate()
            #print('RNA')
            #print(record.seq)
            #print('prot')
            #print(seq)
            ouf.write('>')
            #ouf.write(str(i))
            ouf.write(str(rec))
            ouf.write('\n')
            ouf.write(str(seq))
            #ouf.write(str(i))
            ouf.write('\n')

            #ouf.write('>')
            #ouf.write(str(i))
            #ouf.write(str(rec))
            #ouf.write('\n')


            i = i+1

#messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", generic_rna)
#messenger_rna.translate()


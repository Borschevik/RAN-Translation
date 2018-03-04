from Bio import SeqIO

#input from fasta file
handle = open("names_filtred_short_last_VERYLAST.fasta")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()
print(len(records))

with open ('names_filtred_short_last_VERYLAST.fasta') as file:
    with open ('names_filtred_short_last_VERYLAST_edited.fasta', 'w') as ouf:
        for record in SeqIO.parse(file, "fasta"):
            rec = record.id
            #seq = record.seq.transcribe()
            ouf.write('>')
            ouf.write(str(rec))
            ouf.write('\n')
            #ouf.write('NNNNN')
            #ouf.write(str(seq))
            #ouf.write('\n')
    #seq.translate()

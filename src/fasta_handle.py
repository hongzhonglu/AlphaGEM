from Bio import SeqIO


def rename_genes_in_fasta(input_fasta, output_fasta):
    generecords=SeqIO.parse(input_fasta, "fasta")
    records=[]
    for record in generecords:
        try:
           record.id=record.id.split('SerialNumber')[1]
        except:
            pass
        try:
            record.id='sp|'+record.id.split('|')[1]+'|'+record.id.split('|')[1]
        except:
            record.id='sp|'+record.id+'|'+record.id
        record.name=''
        record.description=''
        record.seq=str(record.seq).replace('*','')
        records.append(record)
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
        output_handle.close()
def handle(fasta,name):
    rename_genes_in_fasta(fasta,f'./working/{name}/{name}.fasta')
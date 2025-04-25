from Bio import SeqIO
from Bio import Seq
from Bio import Align
def predate(name,refname):
    global refgenes
    global tgenes
    refgenes={x.name[3:x.id.find('|',3)]:x for x in SeqIO.parse(f'data_available/{refname}.fasta','fasta')}
    tgenes={x.name[3:x.id.find('|',3)]:x for x in SeqIO.parse(f'./working/{name}/{name}.fasta','fasta')}
def align(nameref,nametar):
    seq1=refgenes[nameref]
    seq2=tgenes[nametar]
    alignment=Align.PairwiseAligner().align(seqA=seq1,seqB=seq2)
    return alignment.score/alignment[0].length
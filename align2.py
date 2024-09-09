from Bio import SeqIO
from Bio import Seq
from Bio import Align
def predate(name,refname):
    global refgenes
    global tgenes
    refgenes=[x.name[3:x.id.find('|',3)] for x in SeqIO.parse(f'ziyuan/{refname}.fasta','fasta')]
    tgenes=[x.name[3:x.id.find('|',3)] for x in SeqIO.parse(f'ziyuan/{name}.fasta','fasta')]
def align(nameref,nametar,name,refname):
    seq1=[x for x in SeqIO.parse(f'ziyuan/{refname}.fasta','fasta')][refgenes.index(nameref)]
    seq2=[x for x in SeqIO.parse(f'ziyuan/{name}.fasta','fasta')][tgenes.index(nametar)]
    alignment=Align.PairwiseAligner().align(seqA=seq1,seqB=seq2)
    return alignment.score/alignment[0].length
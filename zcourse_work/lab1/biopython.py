#imporing module
from Bio.Seq import Seq

#creating Seq object
dna = Seq("ATGAAATTT")
print(f"DNA: {dna}")

#Transcription
rna = dna.transcribe()
print(f"RNA: {rna}")

#Translation
prot = rna.translate()
print(f"Protein: {prot}")

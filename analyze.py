from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


class Analyze:
    def pairwise_alignment(self, seq_1, seq_2):
        """Perform pairwise alignment"""
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        alignments = aligner.align(seq_1, seq_2)
        return alignments

    def parse_record(self, file):
        """Parse fasta file"""
        stringio = StringIO(file.getvalue().decode("utf-8"))
        info = stringio.readline().strip()[1:]
        info_id = info.split("|")[0]
        sequence = []
        for line in stringio:
            if line.startswith(">"):
                info = line[1:].strip()
            else:
                sequence.append(line.strip())
        return {"info": info, "id": info_id, "sequence": "".join(sequence)}

    def seq2fasta(self, sequence, info):
        """Convert sequence to fasta format"""
        seq = Seq(sequence)
        seq_record = SeqRecord(seq, id=info, description="")
        return seq_record.format("fasta")

    def transcribe(self, sequence):
        """Transcribe DNA to mRNA"""
        dna = Seq(sequence)
        return dna.transcribe()

    def translate(self, sequence, to_stop=False):
        """Translate DNA to protein"""
        dna = Seq(sequence)
        return dna.translate(to_stop=to_stop)

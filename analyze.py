from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


class Analyze:
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
        seq_record = SeqRecord(seq, id=info, description=info.split("|")[1])
        return seq_record.format("fasta")

    def transcribe(self, sequence):
        """Transcribe DNA to mRNA"""
        dna = Seq(sequence)
        return dna.transcribe()

    def translate(self, sequence):
        """Translate DNA to protein"""
        dna = Seq(sequence)
        return dna.translate()

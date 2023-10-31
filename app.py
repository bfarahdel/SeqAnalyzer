from Bio.Seq import Seq
from io import StringIO
import streamlit as st


# Parse records
def parse_record(file):
    stringio = StringIO(file.getvalue().decode("utf-8"))
    sequence = []
    for line in stringio:
        if line.startswith(">"):
            info = line[1:].strip()
        else:
            sequence.append(line.strip())
    return {"info": info, "sequence": "".join(sequence)}


# Transcription
def transcribe(sequence):
    dna = Seq(sequence)
    return dna.transcribe()


# Streamlit app
def app():
    st.set_page_config(page_title="SeqAnalyzer", page_icon="🧬", layout="centered")
    st.title("SeqAnalyzer")
    files = st.file_uploader("Upload a fasta file", accept_multiple_files=True)

    if files:
        records = []
        for file in files:
            records.append(parse_record(file))
        sequence_names = [record["info"] for record in records]
        sequence_names.sort()
        selected_sequences = st.multiselect("Select sequences", sequence_names)

        # Transcription
        if st.button("Transcription"):
            if selected_sequences:
                st.write("Transcription Results")
                for record in records:
                    expander = st.expander(label=record["info"])
                    mrna = transcribe(record["sequence"])
                    with expander:
                        st.download_button("Download", str(mrna))
                        st.write(mrna)
            else:
                st.write("Select sequence(s)")


if __name__ == "__main__":
    app()

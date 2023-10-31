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

# Translation
def translate(sequence):
    dna = Seq(sequence)
    return dna.translate()


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

        functions = ["Transcription", "Translation"]
        sorted_functions = sorted(functions)
        function = st.selectbox("Select an option", sorted_functions)

        if st.button("Analyze"):
            # Transcription
            if function == "Transcription":
                if selected_sequences:
                    st.write("Transcription Results")
                    for record in records:
                        record_info = record["info"]
                        expander = st.expander(label=record["info"])
                        mrna = transcribe(record["sequence"])
                        with expander:
                            st.download_button("Download", str(mrna), file_name=f"transcription_{record_info}.txt")
                            st.write(mrna)
                else:
                    st.write("Select sequence(s)")

            # Translation
            if function == "Translation":
                if selected_sequences:
                    st.write("Translation Results")
                    for record in records:
                        record_info = record["info"]
                        expander = st.expander(label=record_info)
                        protein = translate(record["sequence"])
                        with expander:
                            st.download_button("Download", str(protein), file_name=f"translation_{record_info}.txt")
                            st.write(protein)
                else:
                    st.write("Select sequence(s)")



if __name__ == "__main__":
    app()

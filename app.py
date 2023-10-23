import streamlit as st
from Bio import SeqIO


# Function to parse fasta file and return sequence names
def get_sequence_names(file):
    sequence_names = []
    for record in SeqIO.parse(file, "fasta"):
        sequence_names.append(record.id)
    return sequence_names


# Streamlit app
def app():
    st.title("Sequence Analyzer")
    st.write("Upload a fasta file to analyze sequences.")
    files = st.file_uploader("Upload a fasta file", accept_multiple_files=True)

    # Display the sequence name of each file in a multi-select box
    if files:
        sequence_names = []
        for file in files:
            sequence_names += get_sequence_names(file.name)
        sequence_names = list(set(sequence_names))
        sequence_names.sort()
        selected_sequences = st.multiselect("Select sequences", sequence_names)

        # Display the sequence names that were selected
        if selected_sequences:
            st.write("Selected sequences:")
            for sequence in selected_sequences:
                st.write(sequence)


if __name__ == "__main__":
    app()

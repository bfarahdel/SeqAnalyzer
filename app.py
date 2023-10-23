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


# Streamlit app
def app():
    st.title("Sequence Analyzer")
    st.write("Upload a fasta file to analyze sequences.")
    files = st.file_uploader("Upload a fasta file", accept_multiple_files=True)

    # Display the sequence name of each file in a multi-select box
    if files:
        records = []
        for file in files:
            records.append(parse_record(file))
        sequence_names = [record["info"] for record in records]
        sequence_names.sort()
        selected_sequences = st.multiselect("Select sequences", sequence_names)

        # Display the sequence names that were selected
        if selected_sequences:
            st.write("Selected sequences:")
            for sequence in selected_sequences:
                st.write(sequence)


if __name__ == "__main__":
    app()

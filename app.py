from analyze import Analyze
from plot_analysis import PlotAnalysis
import streamlit as st


def app():
    st.set_page_config(page_title="SeqAnalyzer", page_icon="🧬", layout="centered")
    st.title("SeqAnalyzer")
    files = st.file_uploader("Upload a fasta file", accept_multiple_files=True)
    analyze = Analyze()
    plot_analysis = PlotAnalysis()

    if files:
        records = []
        for file in files:
            records.append(analyze.parse_record(file))
        sequence_names = [record["info"] for record in records]
        sequence_names.sort()
        selected_sequences = st.multiselect("Select sequences", sequence_names)

        functions = ["Transcription", "Translation"]
        sorted_functions = sorted(functions)
        function = st.selectbox("Select an option", sorted_functions)

        if st.button("Analyze", use_container_width=True):
            # Transcription
            if function == "Transcription":
                if selected_sequences:
                    st.write("Transcription Results")
                    for record in records:
                        record_info = record["info"]
                        expander = st.expander(label=record["info"])
                        mrna = analyze.transcribe(record["sequence"])
                        mrna_fasta = analyze.seq2fasta(mrna, record["info"])
                        with expander:
                            st.download_button(
                                "Download",
                                mrna_fasta,
                                file_name=f"transcription_{record_info}.fasta",
                            )
                            transcription_plot = plot_analysis.plot_sequences(
                                [{"id": record_info, "seq": mrna}],
                                method="transcription",
                            )
                            st.bokeh_chart(transcription_plot, use_container_width=True)
                else:
                    st.write("Select sequence(s)")

            # Translation
            if function == "Translation":
                if selected_sequences:
                    st.write("Translation Results")
                    for record in records:
                        record_info = record["info"]
                        expander = st.expander(label=record_info)
                        protein = analyze.translate(record["sequence"])
                        protein_fasta = analyze.seq2fasta(protein, record["info"])
                        with expander:
                            st.download_button(
                                "Download",
                                protein_fasta,
                                file_name=f"translation_{record_info}.fasta",
                            )
                            translation_plot = plot_analysis.plot_sequences(
                                [{"id": record_info, "seq": protein}],
                                method="translation",
                            )
                            st.bokeh_chart(translation_plot, use_container_width=True)
                else:
                    st.write("Select sequence(s)")


if __name__ == "__main__":
    app()

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
        records = {}
        for file in files:
            file_parsed = analyze.parse_record(file)
            records[file_parsed["info"]] = file_parsed
        sequence_names = list(records.keys())
        sequence_names.sort()
        selected_sequences = st.multiselect(
            "Select sequences",
            sequence_names,
            placeholder="Select sequences to analyze",
        )

        analyses = ["Sequence Viewer", "Transcription", "Translation"]
        sorted_analyses = sorted(analyses)
        analysis = st.selectbox("Select a method of analysis", sorted_analyses)

        if analysis == "Translation":
            to_stop = st.checkbox("Stop at first stop codon")

        # Sequence Viewer
        if st.button("Analyze", use_container_width=True):
            if analysis == "Sequence Viewer":
                if selected_sequences:
                    st.write("Sequence(s)")
                    for selected in selected_sequences:
                        record = records[selected]
                        record_info = record["info"]
                        expander = st.expander(label=record["info"])
                        with expander:
                            sequence_plot = plot_analysis.plot_sequences(
                                [{"id": record_info, "seq": record["sequence"]}],
                                method="standard",
                            )
                            st.bokeh_chart(sequence_plot, use_container_width=True)
                else:
                    st.write("No sequence(s) selected")

            # Transcription
            if analysis == "Transcription":
                if selected_sequences:
                    st.write("Transcription Results")
                    for selected in selected_sequences:
                        record = records[selected]
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
                    st.write("No sequence(s) selected")

            # Translation
            if analysis == "Translation":
                if selected_sequences:
                    st.write("Translation Results")
                    for selected in selected_sequences:
                        record = records[selected]
                        record_info = record["info"]
                        expander = st.expander(label=record_info)
                        protein = analyze.translate(record["sequence"], to_stop=to_stop)
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
                    st.write("No sequence(s) selected")


if __name__ == "__main__":
    app()

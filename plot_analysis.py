from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.palettes import Category20
import numpy as np


class PlotAnalysis:
    def assign_colors(self, seqs, method="standard"):
        """Assign colors to elements or amino acids"""
        seq_str = [i for s in list(seqs) for i in s]
        if method == "standard":
            color_palette = {
                "A": "red",
                "T": "green",
                "G": "orange",
                "C": "blue",
                "-": "white",
            }
        elif method == "transcription":
            color_palette = {
                "A": "red",
                "U": "green",
                "G": "orange",
                "C": "blue",
                "-": "white",
            }
        elif method == "translation":
            # Assign one letter amino acid codes to each color in Category20 palette
            amino_acids = [
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y",
            ]
            color_palette = {}
            color_palette["*"] = "tan"
            for i, aa in enumerate(amino_acids):
                color_palette[aa] = Category20[20][i]

        colors = [color_palette[i] for i in seq_str]
        return colors

    def mini_plot_sequences(self, input_seqs, method="standard", plot_width=800):
        """View Bokeh plot of the sequences"""
        seqs = [rec["seq"] for rec in (input_seqs)]
        ids = [rec["id"] for rec in input_seqs]
        seq_str = [i for s in list(seqs) for i in s]
        colors = PlotAnalysis().assign_colors(seqs, method=method)
        N = len(seqs[0])
        S = len(seqs)

        x = np.arange(1, N + 1)
        y = np.arange(0, S, 1)
        # Creates a 2D grid from 1D arrays
        xx, yy = np.meshgrid(x, y)
        # Flattens the arrays
        gx = xx.ravel()
        gy = yy.flatten()
        recty = gy + 0.5

        source = ColumnDataSource(
            dict(x=gx, y=gy, recty=recty, text=seq_str, colors=colors)
        )
        x_range = Range1d(0, N + 1, bounds="auto")
        tools = "reset, save"

        plt_seq = figure(
            title=None,
            width=plot_width,
            height=50,
            x_range=x_range,
            y_range=ids,
            tools=tools,
            min_border=0,
            toolbar_location="below",
        )
        element_box = Rect(
            x="x",
            y="recty",
            width=1,
            height=1,
            fill_color="colors",
            line_color=None,
            fill_alpha=0.6,
        )
        plt_seq.add_glyph(source, element_box)
        plt_seq.yaxis.visible = False
        plt_seq.grid.visible = False

        plt_grid = gridplot([[plt_seq]], toolbar_location="below")
        return plt_grid

    def plot_sequences(self, input_seqs, method="standard", plot_width=800):
        """View Bokeh plot of the sequences"""
        seqs = [rec["seq"] for rec in (input_seqs)]
        ids = [rec["id"] for rec in input_seqs]
        text = [i for s in list(seqs) for i in s]
        colors = PlotAnalysis().assign_colors(seqs, method=method)
        N = len(seqs[0])
        S = len(seqs)

        x = np.arange(1, N + 1)
        y = np.arange(0, S, 1)
        # Creates a 2D grid from 1D arrays
        xx, yy = np.meshgrid(x, y)
        # Flattens the arrays
        gx = xx.ravel()
        gy = yy.flatten()
        recty = gy + 0.5
        source = ColumnDataSource(
            dict(x=gx, y=gy, recty=recty, text=text, colors=colors)
        )
        plot_height = len(seqs) * 50 + 100
        view_range = (0, 50)
        tools = "xwheel_pan, xpan, xwheel_zoom, reset, save"

        plt_seq = figure(
            width=plot_width,
            height=plot_height,
            x_range=view_range,
            y_range=ids,
            tools=tools,
            min_border=0,
            toolbar_location="below",
        )
        element = Text(
            x="x",
            y="y",
            text="text",
            text_align="center",
            text_color="black",
            text_font_size="14pt",
        )
        element_box = Rect(
            x="x",
            y="recty",
            width=1,
            height=1,
            fill_color="colors",
            line_color=None,
            fill_alpha=0.4,
        )
        plt_seq.add_glyph(source, element)
        plt_seq.add_glyph(source, element_box)

        plt_seq.grid.visible = False
        plt_seq.yaxis.minor_tick_line_width = 0
        plt_seq.yaxis.major_tick_line_width = 0
        plt_seq.xaxis.major_label_text_font_size = "10pt"
        plt_seq.yaxis.major_label_text_font_size = "0pt"

        plt_grid = gridplot([[plt_seq]], toolbar_location="below")
        return plt_grid

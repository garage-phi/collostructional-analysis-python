import os
import sys
import ipywidgets as widgets
from IPython.display import display
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    # Adjust the import path to include the core module
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
except NameError as e:
    # For Google Colab or other environments
    print(f"If you are using Google Colab or other environment detected,",
          f"please ensure that the core module is accessible.")

try:
    # 1. repository structure with package
    from core.collostructional_analysis import CollostructionalAnalysisMain
except ImportError:
    # 2. Loading in Colab or flat structure (when files are in the same location)
    from collostructional_analysis import CollostructionalAnalysisMain

"""
Variables:
* N: Total Corpus Size
* C: Construction Frequency
* W: Word Frequency
* r: Overlap Ratio (Percentage) (%) ... Simpson's Coefficient
Trajectory Modes:
* Sweep Overlap (fix W, vary a)
* Sweep Word Freq (fix r, vary W)
Memo:
* a: Co-occurrence Frequency
    a = min(C, W) * (r / 100.0)
* b: C - a
* c: W - a
* d: N - C - c
Note on rounding effects:
  Although a, b, c, and d are computed as floating-point values,
  they are rounded to integers for association metric calculations.
  Consequently, the resulting trajectories may deviate slightly
  from those computed using floating-point arithmetic.
"""

# ==========================================
# CollexemeUnitNCWr:
# Data management for a single data point
# ==========================================


class CollexemeUnitNCWr:
    """
    Class: CollexemeUnitNCWr
    It maintains the state and calculation logic for a single word/data point.
    Given N (total corpus size), C% (construction frequency percentage),
    W% (word frequency percentage), and r% (overlap ratio percentage),
    it calculates the actual frequencies and association metrics.
    It also generates trajectory data based on the selected mode.
    1. Sweep Overlap: Fix W, vary a (overlap)
    2. Sweep Word Freq: Fix r, vary W
    """

    def __init__(
        self, name, color, initial_C_pct=20.0, initial_W_pct=5.0, initial_r=50.0
    ):
        self.name = name
        self.color = color

        # State
        self.C_pct = initial_C_pct
        self.W_pct = initial_W_pct
        self.r_pct = initial_r
        self.traj_mode = "Sweep Overlap"

        # Cache
        self.current_metrics = {}
        self.trajectory_df = pd.DataFrame()
        self.a_current = 0
        self.C_val = 0
        self.W_val = 0

    def recalculate(self, N, signed_metrics=True, need_fisher=True):
        """Recalculate values based on the current N"""
        self.C_val = int(N * (self.C_pct / 100.0))
        self.W_val = int(N * (self.W_pct / 100.0))
        if self.W_val == 0:
            self.W_val = 1

        # Constraint: r_pct \in [0, 100]
        if self.r_pct < 0.0:
            self.r_pct = 0.0
        elif self.r_pct > 100.0:
            self.r_pct = 100.0

        # Current Point
        # a calculation: a = min(C, W) * (r / 100.0)
        # Ensure a does not exceed min(C, W)
        self.a_current = int(min(self.C_val, self.W_val) * (self.r_pct / 100.0))
        self.current_metrics = self._calc_single_point(
            N,
            self.C_val,
            self.W_val,
            self.a_current,
            signed_metrics=signed_metrics,
            need_fisher=need_fisher,
        )

        # Trajectory
        self.trajectory_df = self._get_trajectory_data(
            N, signed_metrics=signed_metrics, need_fisher=need_fisher
        )

    def _calc_single_point(
        self, N, C, W, a, signed_metrics: bool = True, need_fisher: bool = False
    ):
        b = C - a
        c = W - a
        d = N - C - c
        res = AssociationStatsKernel.calculate_all_metrics(
            a,
            b,
            c,
            d,
            total_corpus_size=N,
            signed_metrics=signed_metrics,
            include_fisher=need_fisher,
        )

        # Clearly label and store values here
        res["a"] = a
        res["W"] = W
        res["C"] = C
        # For convenience, also calculate and store the Overlap Ratio
        res["Overlap_Ratio_Calc"] = (a / min(W, C) * 100) if min(W, C) > 0 else 0
        return res

    def _get_trajectory_data(
        self, N, signed_metrics: bool = True, need_fisher: bool = False
    ):
        results = []
        mode = self.traj_mode

        if mode == "Sweep Overlap":
            # Word Freq (W) fixed, vary Overlap (a)
            # -> Taking W as the axis should result in "points"
            max_a = min(self.C_val, self.W_val)
            a_values = np.linspace(0, max_a, 200).astype(int)
            a_values = np.unique(a_values)
            for a in a_values:
                # W: Fixxed (self.W_val)
                results.append(
                    self._calc_single_point(
                        N,
                        self.C_val,
                        self.W_val,
                        a,
                        signed_metrics=signed_metrics,
                        need_fisher=need_fisher,
                    )
                )

        elif mode == "Sweep Word Freq":
            # Overlap Ratio (r) fixed, vary Word Freq (W)
            # -> Taking W as the axis should result in a curve
            w_values = np.unique(np.geomspace(1, N, 300).astype(int))
            for w in w_values:
                max_possible_a = min(self.C_val, w)
                a = int(max_possible_a * (self.r_pct / 100.0))
                if (self.C_val + w - a) > N:
                    continue

                # W varies
                results.append(
                    self._calc_single_point(
                        N,
                        self.C_val,
                        w,
                        a,
                        signed_metrics=signed_metrics,
                        need_fisher=need_fisher,
                    )
                )

        return pd.DataFrame(results)


# ==========================================
# CollostructionalComparisonSimulatorNCWr:
# Simulator UI for comparing two data points
# ==========================================


class CollostructionalComparisonSimulatorNCWr:
    """
    Class: CollostructionalComparisonSimulatorNCWr
    
    If you want to change the initial values, modify the __init__ method.
    If you want to change the initial modes, modify 
    the self.w_A_mode and self.w_B_mode default values in create_widgets().
    """

    def __init__(self):
        self.unit_A = CollexemeUnitNCWr(
            "Data A", "#4477AA", initial_C_pct=20.0, initial_W_pct=5.0
        )
        self.unit_B = CollexemeUnitNCWr(
            "Data B", "#EE6677", initial_C_pct=20.0, initial_W_pct=10.0
        )

        # Clearly define axis options (Remove Freq Raw, separate into a, W, C)
        self.metric_options = {
            "Fisher (-log10 p)": "FYE",
            "Log Likelihood (G2)": "LLR",
            "Log Odds Ratio": "LOGODDSRATIO",
            "PMI": "PMI",
            "Pearson Residuals": "PEARSONRESID",
            "Delta P (C->W)": "DELTAPC2W",
            "Co-occurrence (a)": "a",  # Co-occurrence count
            "Word Freq (W)": "W",  # Word frequency
            "Const Freq (C)": "C",  # Construction frequency
        }
        self.traj_options = ["Sweep Overlap", "Sweep Word Freq"]

        self.create_widgets()

    def create_widgets(self):
        style = {"description_width": "100px"}
        slider_layout = widgets.Layout(width="98%")

        # --- Global ---
        self.w_N = widgets.IntText(value=10000, description="Corpus N:", style=style)

        # --- Unit A Widgets ---
        self.w_A_C = widgets.FloatSlider(
            value=20.0,
            min=0.1,
            max=100.0,
            description="Const %:",
            style=style,
            layout=slider_layout,
        )
        self.w_A_W = widgets.FloatSlider(
            value=5.0,
            min=0.1,
            max=100.0,
            description="Word %:",
            style=style,
            layout=slider_layout,
        )
        self.w_A_r = widgets.FloatSlider(
            value=50.0,
            min=0.0,
            max=100.0,
            description="Overlap r%:",
            style=style,
            layout=slider_layout,
        )
        self.w_A_mode = widgets.Dropdown(
            options=self.traj_options,
            value="Sweep Overlap",
            description="Mode:",
            style=style,
            layout=widgets.Layout(width="90%"),
        )

        # --- Unit B Widgets ---
        self.w_B_C = widgets.FloatSlider(
            value=20.0,
            min=0.1,
            max=100.0,
            description="Const %:",
            style=style,
            layout=slider_layout,
        )
        self.w_B_W = widgets.FloatSlider(
            value=10.0,
            min=0.1,
            max=100.0,
            description="Word %:",
            style=style,
            layout=slider_layout,
        )
        self.w_B_r = widgets.FloatSlider(
            value=50.0,
            min=0.0,
            max=100.0,
            description="Overlap r%:",
            style=style,
            layout=slider_layout,
        )
        self.w_B_mode = widgets.Dropdown(
            options=self.traj_options,
            value="Sweep Overlap",
            description="Mode:",
            style=style,
            layout=widgets.Layout(width="90%"),
        )

        # --- Axis Config ---
        self.w_xaxis = widgets.Dropdown(
            options=self.metric_options, value="PMI", description="X Axis:", style=style
        )
        self.w_yaxis = widgets.Dropdown(
            options=self.metric_options, value="FYE", description="Y Axis:", style=style
        )
        self.w_detail = widgets.Checkbox(value=True, description="Show Details")
        self.w_signed = widgets.Checkbox(value=True, description="Signed Metrics")

        self.out_view = widgets.Output()

        # Observers
        widget_list = [
            self.w_N,
            self.w_xaxis,
            self.w_yaxis,
            self.w_detail,
            self.w_signed,
            self.w_A_C,
            self.w_A_W,
            self.w_A_r,
            self.w_A_mode,
            self.w_B_C,
            self.w_B_W,
            self.w_B_r,
            self.w_B_mode,
        ]
        for w in widget_list:
            w.observe(self.update_plot, names="value")

    def update_state(self, signed_metrics=True, need_fisher=True):
        N = self.w_N.value
        # Unit A
        self.unit_A.C_pct = self.w_A_C.value
        self.unit_A.W_pct = self.w_A_W.value
        self.unit_A.r_pct = self.w_A_r.value
        self.unit_A.traj_mode = self.w_A_mode.value
        self.unit_A.recalculate(
            N, signed_metrics=signed_metrics, need_fisher=need_fisher
        )
        # Unit B
        self.unit_B.C_pct = self.w_B_C.value
        self.unit_B.W_pct = self.w_B_W.value
        self.unit_B.r_pct = self.w_B_r.value
        self.unit_B.traj_mode = self.w_B_mode.value
        self.unit_B.recalculate(
            N, signed_metrics=signed_metrics, need_fisher=need_fisher
        )

    def find_peaks(self, df, x_col, y_col):
        if len(df) < 5:
            return []
        y = df[y_col].values
        x = df[x_col].values
        
        if np.var(y) < 1e-9 and np.var(x) < 1e-9:
            return []

        # Since a, b, c, and d are rounded to integers, 
        # small fluctuations may occur.
        # we assess the range using only the central 90% of the data.
        # If the main body is effectively flat, 
        # we skip peak detection to ignore these artifacts.

        cut_len = int(len(y) * 0.05)  # side:5%
        if cut_len > 0:
            # slice off 5% from both ends
            y_middle = y[cut_len : -cut_len]
        else:
            y_middle = y

        if len(y_middle) == 0:
            return []

        y_range = np.max(y_middle) - np.min(y_middle)
        if y_range < 1e-2: 
            return []

        peaks = []
        grads = np.gradient(y)
        for i in range(1, len(grads) - 1):
            if grads[i - 1] > 0 and grads[i] < 0:
                peaks.append((x[i], y[i], "Max"))
            elif grads[i - 1] < 0 and grads[i] > 0:
                peaks.append((x[i], y[i], "Min"))
        return peaks

    def update_plot(self, change=None):
        self.out_view.clear_output(wait=True)

        x_key = self.w_xaxis.value
        y_key = self.w_yaxis.value
        need_fisher = True  # for detailed stats
        # need_fisher = x_key == "FYE" or y_key == "FYE" # You can use this logic if needed

        signed_metrics = self.w_signed.value
        self.update_state(signed_metrics=signed_metrics, need_fisher=need_fisher)

        with self.out_view:
            fig, ax = plt.subplots(figsize=(12, 7))

            def plot_unit(unit, marker_curr):
                df = unit.trajectory_df
                if not df.empty:
                    # Clean data: remove NaN and Inf
                    df_clean = df.replace([np.inf, -np.inf], np.nan).dropna(
                        subset=[x_key, y_key]
                    )

                    if len(df_clean) > 0:
                        # Check if all points have the same value (switch to point plot)
                        is_single_point = (df_clean[x_key].nunique() <= 1) and (
                            df_clean[y_key].nunique() <= 1
                        )

                        label_txt = f"{unit.name}: {unit.traj_mode}"

                        if is_single_point:
                            """
                            When the trajectory is a "point" 
                            (e.g., when plotting fixed variables),
                            we use scatter to plot a single point.
                            """
                            ax.scatter(
                                df_clean[x_key].iloc[0],
                                df_clean[y_key].iloc[0],
                                color=unit.color,
                                s=50,
                                alpha=0.5,
                                label=label_txt + " (Fixed)",
                            )
                        else:
                            # When the trajectory is not a point
                            ax.plot(
                                df_clean[x_key],
                                df_clean[y_key],
                                "-",
                                color=unit.color,
                                linewidth=2,
                                alpha=0.7,
                                label=label_txt,
                            )
                            # Peaks
                            peaks = self.find_peaks(df_clean, x_key, y_key)
                            for px, py, ptype in peaks:
                                ax.scatter(px, py, color=unit.color, marker="x", s=60)

                # Current Point
                curr = unit.current_metrics
                cx, cy = curr.get(x_key, np.nan), curr.get(y_key, np.nan)
                if not (np.isinf(cx) or np.isinf(cy) or np.isnan(cx) or np.isnan(cy)):
                    ax.scatter(
                        [cx],
                        [cy],
                        color=unit.color,
                        s=150,
                        edgecolors="black",
                        zorder=10,
                    )
                    ax.annotate(
                        f"{unit.name}\n({cx:.2f}, {cy:.2f})",
                        (cx, cy),
                        xytext=(10, 10),
                        textcoords="offset points",
                        color=unit.color,
                        fontweight="bold",
                    )

            plot_unit(self.unit_A, "o")
            plot_unit(self.unit_B, "D")

            if self.w_detail.value:

                def get_stats_txt(unit):
                    m = unit.current_metrics
                    # Handle case when Fisher is not calculated (NaN)
                    fye_val = m.get("FYE", np.nan)
                    fye_str = f"{fye_val:.2f}" if pd.notna(fye_val) else "---"

                    return (
                        f"[{unit.name}]\n"
                        f"C={m['C']:,} ({unit.C_pct}%), W={m['W']:,} ({unit.W_pct}%)\n"
                        f"a={m['a']:,} (r={unit.r_pct}%)\n"
                        f"LLR: {m['LLR']:.2f}, PEARSONRESID: {m['PEARSONRESID']:.2f}\n"
                        f"LOGODDSRATIO: {m['LOGODDSRATIO']:.2f}, PMI: {m['PMI']:.2f}\n"
                        f"FYE: {fye_str}, DELTAPC2W: {m['DELTAPC2W']:.2f}\n"
                    )

                info_txt = (
                    f"N={self.w_N.value:,}\n\n"
                    + get_stats_txt(self.unit_A)
                    + "\n"
                    + get_stats_txt(self.unit_B)
                )
                props = dict(boxstyle="round", facecolor="white", alpha=0.9)
                ax.text(
                    1.02,
                    1.0,
                    info_txt,
                    transform=ax.transAxes,
                    fontsize=9,
                    verticalalignment="top",
                    bbox=props,
                )

            ax.set_xlabel(self.w_xaxis.label)
            ax.set_ylabel(self.w_yaxis.label)
            ax.grid(True, linestyle=":", alpha=0.6)
            ax.legend()
            ax.set_title("Collostructional Comparator: Two Data Points")
            plt.tight_layout()
            plt.show()

    def display(self):
        box_layout = widgets.Layout(
            border="1px solid #ddd",
            padding="10px",
            margin="5px",
            width="48%",
            min_width="300px",
        )

        box_A = widgets.VBox(
            [
                widgets.HTML("<b>Data A (Blue)</b>"),
                self.w_A_C,
                self.w_A_W,
                self.w_A_r,
                self.w_A_mode,
            ],
            layout=box_layout,
        )

        box_B = widgets.VBox(
            [
                widgets.HTML("<b>Data B (Red)</b>"),
                self.w_B_C,
                self.w_B_W,
                self.w_B_r,
                self.w_B_mode,
            ],
            layout=box_layout,
        )

        ui = widgets.VBox(
            [
                widgets.HTML("<h3>Collostructional Comparator (Fixed)</h3>"),
                self.w_N,
                widgets.HBox([box_A, box_B]),
                widgets.HBox(
                    [self.w_xaxis, self.w_yaxis, self.w_detail, self.w_signed]
                ),
                self.out_view,
            ]
        )
        display(ui)
        self.update_plot()


"""
# Execute
sim = CollostructionalComparisonSimulatorNCWr()
sim.display()
"""

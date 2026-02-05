import os
import sys
import ipywidgets as widgets
from IPython.display import display
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

try:
    # 1. repository structure with package
    from core.collostructional_analysis import CollostructionalAnalysisMain
except ImportError:
    # 2. Loading in Colab or flat structure (when files are in the same location)
    from collostructional_analysis import CollostructionalAnalysisMain

"""
Variables:
* K: Const Col 1 Total (a + c)
* M: Const Col 2 Total (b + d)
* r: Ratio (a / b)
* b: Value b
Trajectory Modes:
* Sweep Ratio r (fix b, vary r)
* Sweep b (fix r, vary b)
Memo:
* a: Co-occurrence Frequency
    a = r * b
* c: K - a
* d: M - b
* N: Total Corpus Size
    N = K + M
Note on rounding effects:
  Although a, b, c, and d are computed as floating-point values,
  they are rounded to integers for association metric calculations.
  Consequently, the resulting trajectories may deviate slightly
  from those computed using floating-point arithmetic.
"""

# ==========================================
# CollexemeUnitKMbr: 
# Data management for a single data point
# ==========================================

class CollexemeUnitKMbr:
    """
    Class: CollexemeUnitKMbr
    It maintains the state and calculation logic for a single word/data point.
    Given K, M, b, r, it calculates a, c, d, N and various association metrics.
    It also computes the trajectory data based on 
    the selected mode (Sweep Ratio r or Sweep b).
    Variables:
    * K: Const Col 1 Total (a + c)
    * M: Const Col 2 Total (b + d)
    * r: Ratio (a / b)
    * b: Value b
    Trajectory Modes:
    * Sweep Ratio r (fix b, vary r)
    * Sweep b (fix r, vary b)
    """

    def __init__(
        self, name, color, initial_K=100, initial_M=1000, initial_b=10, initial_r=1.0
    ):
        self.name = name
        self.color = color

        # Inputs
        self.K_val = initial_K  # col 1 total (a+c)
        self.M_val = initial_M  # col 2 total (b+d)
        self.b_val = initial_b
        self.r_val = initial_r  # ratio a/b

        self.traj_mode = "Sweep Ratio r"  # or 'Sweep b'

        # Cache
        self.current_metrics = {}
        self.trajectory_df = pd.DataFrame()

        # Internal calculated values
        self.a_calc = 0
        self.c_calc = 0
        self.d_calc = 0
        self.N_calc = 0

    def recalculate(self, signed_metrics=True, need_fisher=True):
        """Recalculate values based on the current N"""

        # --- Current Point Calculation ---
        self.N_calc = self.K_val + self.M_val

        # a = r * b
        raw_a = self.r_val * self.b_val

        # Constraints: a <= K, b <= M
        # Even if the UI restricts values, 
        # clip here to avoid inconsistencies in calculations
        # Ensure b does not exceed M
        valid_b = min(self.b_val, self.M_val)
        valid_a = min(raw_a, self.K_val)

        self.a_calc = valid_a
        self.b_calc = valid_b
        self.c_calc = self.K_val - self.a_calc
        self.d_calc = self.M_val - self.b_calc

        # Metrics calculation
        self.current_metrics = self._calc_single_point(
            self.a_calc,
            self.b_calc,
            self.c_calc,
            self.d_calc,
            self.N_calc,
            signed_metrics=signed_metrics,
            need_fisher=need_fisher,
        )

        # --- Trajectory Calculation ---
        self.trajectory_df = self._get_trajectory_data(
            signed_metrics=signed_metrics, need_fisher=need_fisher
        )

    def _calc_single_point(self, a, b, c, d, N, signed_metrics=True, need_fisher=False):
        # Convert a, b, c, d to int (for Fisher's exact test, etc.)
        # MEMO: Floating-point values may be used for smoothing purposes.
        # However, we round them to int here for safety.
        # Side effects: The trajectory may have slightly different values
        # due to rounding.
        ai, bi, ci, di = int(round(a)), int(round(b)), int(round(c)), int(round(d))
        Ni = ai + bi + ci + di

        # Avoid zero total corpus size
        if Ni == 0:
            Ni = 1

        res = AssociationStatsKernel.calculate_all_metrics(
            ai,
            bi,
            ci,
            di,
            total_corpus_size=Ni,
            signed_metrics=signed_metrics,
            include_fisher=need_fisher,
        )

        # Save parameters for axis usage
        res["a"] = a
        res["b"] = b
        res["c"] = c
        res["d"] = d
        res["K"] = a + c
        res["M"] = b + d
        res["r"] = (a / b) if b > 0 else 0
        return res

    def _get_trajectory_data(self, signed_metrics=True, need_fisher=False):
        results = []
        mode = self.traj_mode

        K = self.K_val
        M = self.M_val
        N = K + M

        if mode == "Sweep Ratio r":
            # Fix b, vary r
            b_fixed = min(self.b_val, M)
            # If b <= 0, set to small value to avoid division by zero
            if b_fixed <= 0:
                b_fixed = 1e-9

            max_r = K / b_fixed
            # make 200 steps
            r_values = np.linspace(0, max_r, 200)

            for r_i in r_values:
                a_i = r_i * b_fixed
                c_i = K - a_i
                d_i = M - b_fixed  # d is a constant

                results.append(
                    self._calc_single_point(
                        a_i,
                        b_fixed,
                        c_i,
                        d_i,
                        N,
                        signed_metrics=signed_metrics,
                        need_fisher=need_fisher,
                    )
                )

        elif mode == "Sweep b":
            # Fix r, vary b
            r_fixed = self.r_val
            if r_fixed < 0:
                r_fixed = 0

            limit_by_K = (K / r_fixed) if r_fixed > 1e-9 else M
            max_b = min(M, limit_by_K)

            # Make 200 steps (excluding 0)
            b_values = np.linspace(0.1, max_b, 200)

            for b_i in b_values:
                a_i = r_fixed * b_i
                c_i = K - a_i
                d_i = M - b_i

                results.append(
                    self._calc_single_point(
                        a_i,
                        b_i,
                        c_i,
                        d_i,
                        N,
                        signed_metrics=signed_metrics,
                        need_fisher=need_fisher,
                    )
                )

        return pd.DataFrame(results)


# ==========================================
# CollostructionalComparisonSimulatorKMbr:
# Simulator UI for comparing two data points
# ==========================================


class CollostructionalComparisonSimulatorKMbr:
    """
    Class: CollostructionalComparisonSimulatorKMbr

    If you want to change the initial values, modify the __init__ method.
    If you want to change the initial modes, modify 
    the self.w_A_mode and self.w_B_mode default values in create_widgets().
    """

    def __init__(self):
        # Initialize two data units
        self.unit_A = CollexemeUnitKMbr(
            "Data A",
            "#4477AA",
            initial_K=250,
            initial_M=750,
            initial_b=35,
            initial_r=3.00,
        )
        self.unit_B = CollexemeUnitKMbr(
            "Data B",
            "#EE6677",
            initial_K=250,
            initial_M=750,
            initial_b=70,
            initial_r=2.50,
        )

        # Axis options
        self.metric_options = {
            "Fisher (-log10 p)": "FYE",
            "Log Likelihood (G2)": "LLR",
            "Log Odds Ratio": "LOGODDSRATIO",
            "PMI": "PMI",
            "Pearson Residuals": "PEARSONRESID",
            "Delta P (C->W)": "DELTAPC2W",
            "Co-occurrence (a)": "a",
            "Value b": "b",
            "Ratio r (a/b)": "r",
            "Const K (a+c)": "K",
            "Const M (b+d)": "M",
        }
        self.traj_options = ["Sweep Ratio r", "Sweep b"]

        self.create_widgets()

    def create_widgets(self):
        style = {"description_width": "80px"}
        slider_layout = widgets.Layout(width="95%")
        input_layout = widgets.Layout(width="90%")

        # --- Unit A Widgets ---
        # K, M are IntText for easier input
        self.w_A_K = widgets.IntText(
            value=self.unit_A.K_val,
            description="K (a+c):",
            style=style,
            layout=input_layout,
        )
        self.w_A_M = widgets.IntText(
            value=self.unit_A.M_val,
            description="M (b+d):",
            style=style,
            layout=input_layout,
        )
        # b, r are FloatText and FloatSlider
        self.w_A_b = widgets.FloatText(
            value=self.unit_A.b_val, description="b:", 
            style=style, layout=input_layout
        )
        self.w_A_r = widgets.FloatSlider(
            value=self.unit_A.r_val,
            min=0.01,
            max=10.0,
            step=0.1,
            description="r (a/b):",
            style=style,
            layout=slider_layout,
        )
        self.w_A_mode = widgets.Dropdown(
            options=self.traj_options,
            value="Sweep Ratio r",
            description="Mode:",
            style=style,
            layout=input_layout,
        )

        # --- Unit B Widgets ---
        self.w_B_K = widgets.IntText(
            value=self.unit_B.K_val,
            description="K (a+c):",
            style=style,
            layout=input_layout,
        )
        self.w_B_M = widgets.IntText(
            value=self.unit_B.M_val,
            description="M (b+d):",
            style=style,
            layout=input_layout,
        )
        self.w_B_b = widgets.FloatText(
            value=self.unit_B.b_val, description="b:", style=style, layout=input_layout
        )
        self.w_B_r = widgets.FloatSlider(
            value=self.unit_B.r_val,
            min=0.01,
            max=10.0,
            step=0.1,
            description="r (a/b):",
            style=style,
            layout=slider_layout,
        )
        self.w_B_mode = widgets.Dropdown(
            options=self.traj_options,
            value="Sweep Ratio r", # it may be interesting to set "Sweep b" here
            description="Mode:",
            style=style,
            layout=input_layout,
        )

        # --- Axis Config ---
        self.w_xaxis = widgets.Dropdown(
            options=self.metric_options, value="PMI", 
            description="X Axis:", style=style
        )
        self.w_yaxis = widgets.Dropdown(
            options=self.metric_options, value="LOGODDSRATIO", 
            description="Y Axis:", style=style
        )
        self.w_detail = widgets.Checkbox(value=True, description="Show Stats Box")
        self.w_signed = widgets.Checkbox(value=True, description="Signed Metrics")

        self.out_view = widgets.Output()

        # Observers
        widget_list = [
            self.w_xaxis,
            self.w_yaxis,
            self.w_detail,
            self.w_signed,
            self.w_A_K,
            self.w_A_M,
            self.w_A_b,
            self.w_A_r,
            self.w_A_mode,
            self.w_B_K,
            self.w_B_M,
            self.w_B_b,
            self.w_B_r,
            self.w_B_mode,
        ]
        for w in widget_list:
            w.observe(self.update_plot, names="value")

    def update_state(self, signed_metrics=True, need_fisher=True):
        # Unit A
        self.unit_A.K_val = self.w_A_K.value
        self.unit_A.M_val = self.w_A_M.value
        self.unit_A.b_val = self.w_A_b.value
        self.unit_A.r_val = self.w_A_r.value
        self.unit_A.traj_mode = self.w_A_mode.value
        self.unit_A.recalculate(signed_metrics=signed_metrics, need_fisher=need_fisher)

        # Unit B
        self.unit_B.K_val = self.w_B_K.value
        self.unit_B.M_val = self.w_B_M.value
        self.unit_B.b_val = self.w_B_b.value
        self.unit_B.r_val = self.w_B_r.value
        self.unit_B.traj_mode = self.w_B_mode.value
        self.unit_B.recalculate(signed_metrics=signed_metrics, need_fisher=need_fisher)

    def find_peaks(self, df, x_col, y_col):
        if len(df) < 5:
            return []
        y = df[y_col].values
        x = df[x_col].values
        if np.var(y) < 1e-9:
            return []

        peaks = []
        grads = np.gradient(y)
        for i in range(1, len(grads) - 1):
            if grads[i - 1] > 0 and grads[i] < 0:
                peaks.append((x[i], y[i]))
            elif grads[i - 1] < 0 and grads[i] > 0:
                peaks.append((x[i], y[i]))
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
            fig, ax = plt.subplots(figsize=(10, 6))

            def plot_unit(unit, marker_curr):
                df = unit.trajectory_df
                if not df.empty:
                    # Clean data: remove NaN and Inf
                    df_clean = df.replace([np.inf, -np.inf], np.nan).dropna(
                        subset=[x_key, y_key]
                    )

                    if len(df_clean) > 0:
                        label_txt = f"{unit.name}: {unit.traj_mode}"
                        # Check if all points have the same value (switch to point plot)
                        is_single_point = (df_clean[x_key].nunique() <= 1) and (
                            df_clean[y_key].nunique() <= 1
                        )

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
                            for px, py in peaks:
                                ax.scatter(px, py, color=unit.color, marker="x", s=60)

                # Current Point
                curr = unit.current_metrics
                cx, cy = curr.get(x_key, np.nan), curr.get(y_key, np.nan)
                if (
                    pd.notna(cx)
                    and pd.notna(cy)
                    and not np.isinf(cx)
                    and not np.isinf(cy)
                ):
                    ax.scatter(
                        [cx],
                        [cy],
                        color=unit.color,
                        s=120,
                        edgecolors="black",
                        zorder=10,
                        marker=marker_curr,
                    )
                    ax.annotate(
                        f"{unit.name}",
                        (cx, cy),
                        xytext=(5, 5),
                        textcoords="offset points",
                        color=unit.color,
                        fontweight="bold",
                        fontsize=9,
                    )

            plot_unit(self.unit_A, "o")
            plot_unit(self.unit_B, "D")

            if self.w_detail.value:

                def get_stats_txt(unit):
                    m = unit.current_metrics

                    fye_val = m.get('FYE', np.nan)
                    fye_str = f"{fye_val:.2f}" if pd.notna(fye_val) else "---"

                    llr = m.get("LLR", 0)
                    pmi = m.get("PMI", 0)

                    return (
                        f"[{unit.name}]\n"
                        f"K={int(m['K']):,}, M={int(m['M']):,}\n"
                        f"b={m['b']:.1f}, r={m['r']:.2f}\n"
                        f"-> a={m['a']:.1f}, c={m['c']:.1f}, d={m['d']:.1f}\n"
                        f"-------------------\n"
                        f"LLR: {m['LLR']:.2f}, PEARSONRESID: {m['PEARSONRESID']:.2f}\n"
                        f"LOGODDSRATIO: {m['LOGODDSRATIO']:.2f}, PMI: {m['PMI']:.2f}\n"
                        f"FYE: {fye_str}, DELTAPC2W: {m['DELTAPC2W']:.2f}\n"
                    )

                info_txt = (
                    get_stats_txt(self.unit_A) + "\n\n" + get_stats_txt(self.unit_B)
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
            ax.legend(loc="upper left", bbox_to_anchor=(0, -0.15), ncol=2)
            ax.set_title(f"Collexeme Dynamics (K, M Fixed / Varying r or b)")
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
                self.w_A_K,
                self.w_A_M,
                self.w_A_b,
                self.w_A_r,
                self.w_A_mode,
            ],
            layout=box_layout,
        )

        box_B = widgets.VBox(
            [
                widgets.HTML("<b>Data B (Red)</b>"),
                self.w_B_K,
                self.w_B_M,
                self.w_B_b,
                self.w_B_r,
                self.w_B_mode,
            ],
            layout=box_layout,
        )

        ui = widgets.VBox(
            [
                widgets.HTML("<h3>Collostructional Simulator (K-M-b-r Model)</h3>"),
                widgets.HBox([box_A, box_B]),
                widgets.HTML("<hr>"),
                widgets.HBox(
                    [self.w_xaxis, self.w_yaxis, self.w_detail, self.w_signed]
                ),
                self.out_view,
            ]
        )
        display(ui)
        self.update_plot()


# Execute
# sim = CollostructionalComparisonSimulatorKMbr()
# sim.display()

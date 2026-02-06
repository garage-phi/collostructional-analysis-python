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
* N: Const Col 1 Total (a + b + c + d)
* s: Const col 2 Ratio (c / d)
* r: Ratio (a / b)
* b: Value b
Trajectory Modes:
* Sweep Ratio r (fix b, vary r)
* Sweep b (fix r, vary b)
Memo:
* a: Co-occurrence Frequency
    a = r * b
* d: {N - (1 + r) * b } / (1 + s)
* c: s * d
Note on rounding effects:
  Although a, b, c, and d are computed as floating-point values,
  they are rounded to integers for association metric calculations.
  Consequently, the resulting trajectories may deviate slightly
  from those computed using floating-point arithmetic.
"""

# ==========================================
# CollexemeUnitNsbr: 
# Data management for a single data point
# ==========================================

class CollexemeUnitNsbr:
    """
    Class: CollexemeUnitNsbr
    It maintains the state and calculation logic for a single word/data point.
    Given N, s, b, r, it calculates a, c, d, N and various association metrics.
    It also computes the trajectory data based on 
    the selected mode (Sweep Ratio r or Sweep b).
    Variables:
    * N: Const Col 1 Total (a + b + c + d)
    * s: Const col 2 Ratio (c / d)
    * r: Ratio (a / b)
    * b: Value b
    Trajectory Modes:
    * Sweep Ratio r (fix b, vary r)
    * Sweep b (fix r, vary b)
    """

    def __init__(
        self, name, color, initial_N=10000, initial_s=1.0, initial_b=10, initial_r=1.0
    ):
        self.name = name
        self.color = color

        # Inputs
        self.N_val = initial_N  # total corpus size (a + b + c + d)
        self.s_val = initial_s  # ratio (c / d)
        self.b_val = initial_b
        self.r_val = initial_r  # ratio a/b

        self.traj_mode = "Sweep b"  # or 'Sweep Ratio r'

        # Cache
        self.current_metrics = {}
        self.trajectory_df = pd.DataFrame()

        # Internal calculated values
        valid_a = max(0, self.r_val * self.b_val)
        valid_d = max(0, (self.N_val - (1 + self.r_val) * self.b_val) / (1 + self.s_val))
        valid_c = max(0, self.s_val * valid_d)

        self.a_calc = valid_a
        self.c_calc = valid_c
        self.d_calc = valid_d

    def recalculate(self, signed_metrics=True, need_fisher=True):
        """Recalculate values based on the current N"""

        # --- Current Point Calculation ---
        self.N_calc = self.N_val

        # Constraints: a, b, c, d >= 0
        valid_r = max(0, self.r_val)
        self.r_val = valid_r
        valid_s = max(0, self.s_val)
        self.s_val = valid_s

        # for c >= 0:
        valid_b = min(self.b_val, self.N_calc / (1 + self.r_val))
        valid_b = max(0, valid_b)
        self.b_calc = valid_b
        self.b_val = valid_b

        # a = r * b
        raw_a = self.r_val * self.b_calc
        valid_a = min(raw_a, self.N_calc)
        valid_a = max(0, valid_a)
        self.a_calc = valid_a

        # c >= 0  =>  N - (1 + r) * b  >= 0  =>  b <= N / (1 + r)
        raw_d = (self.N_calc - (1 + self.r_val) * self.b_val) / (1 + self.s_val)
        valid_d = min(raw_d, self.N_calc)
        valid_d = max(0, valid_d)
        self.d_calc = valid_d

        # d = s * c
        raw_c = self.s_val * self.d_calc
        valid_c = min(raw_c, self.N_calc)
        valid_c = max(0, valid_c)
        self.c_calc = valid_c


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
        res["N"] = N
        res["r"] = (a / b) if b > 0 else 0
        res["s"] = (c / d) if d > 0 else 0
        return res

    def _get_trajectory_data(self, signed_metrics=True, need_fisher=False):
        results = []
        mode = self.traj_mode

        a = self.b_val * self.r_val
        d = (self.N_val - (1 + self.r_val) * self.b_val) / (1 + self.s_val)
        c = self.s_val * d

        N = self.N_val
        K = a + c


        if mode == "Sweep Ratio r":
            # Fix b, vary r
            b_fixed = min(self.b_val, N)
            # If b <= 0, set to small value to avoid division by zero
            if b_fixed <= 0:
                b_fixed = 1e-9

            # a + b <= N  =>  a <= N - b  =>  r <= (N - b) / b
            max_r = (N / b_fixed) - 1
            # make 200 steps
            r_values = np.linspace(0, max_r, 200)

            for r_i in r_values:
                a_i = r_i * b_fixed
                d_i = (N - (1 + r_i) * b_fixed) / (1 + self.s_val)
                c_i = self.s_val * d_i

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

            # a + b <= N  =>  b <= N - a  =>  b <= N - r * b 
            #  =>  b * (1 + r) <= N  =>  b <= N / (1 + r)
            max_b = N / (1 + r_fixed)

            # Make 200 steps (excluding 0)
            b_values = np.linspace(0.1, max_b, 200)

            for b_i in b_values:
                a_i = r_fixed * b_i
                d_i = (N - (1 + r_fixed) * b_i) / (1 + self.s_val)
                c_i = self.s_val * d_i

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
# CollostructionalComparisonSimulatorNsbr:
# Simulator UI for comparing two data points
# ==========================================


class CollostructionalComparisonSimulatorNsbr:
    """
    Class: CollostructionalComparisonSimulatorNsbr

    If you want to change the initial values, modify the __init__ method.
    If you want to change the initial modes, modify 
    the self.w_A_mode and self.w_B_mode default values in create_widgets().
    """

    def __init__(self):
        # Initialize two data units
        self.unit_A = CollexemeUnitNsbr(
            "Data A",
            "#4477AA",
            initial_N=10000,
            initial_b=1500,
            initial_s=1.00,
            initial_r=3.00,
        )
        self.unit_B = CollexemeUnitNsbr(
            "Data B",
            "#EE6677",
            initial_N=10000,
            initial_b=1500,
            initial_s=1.00,
            initial_r=3.00,
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
        # N: total corpus size
        self.w_A_N = widgets.IntText(
            value=self.unit_A.N_val,
            description="N (Total):",
            style=style,
            layout=input_layout,
        )

        # s: FloatSlider for ratio c/d
        self.w_A_s = widgets.FloatSlider(
            value=self.unit_A.s_val,
            min=0.01,
            max=10.0,
            step=0.1,
            description="s (c/d):",
            style=style,
            layout=slider_layout,
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
        self.w_B_N = widgets.IntText(
            value=self.unit_B.N_val,
            description="N (Total):",
            style=style,
            layout=input_layout,
        )
        self.w_B_s = widgets.FloatSlider(
            value=self.unit_B.s_val,
            min=0.01,
            max=10.0,
            step=0.1,
            description="s (c/d):",
            style=style,
            layout=slider_layout,
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
            value="Sweep b", # it may be interesting to set "Sweep b" here
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
            self.w_A_N,
            self.w_A_s,
            self.w_A_b,
            self.w_A_r,
            self.w_A_mode,
            self.w_B_N,
            self.w_B_s,
            self.w_B_b,
            self.w_B_r,
            self.w_B_mode,
        ]
        for w in widget_list:
            w.observe(self.update_plot, names="value")

    def update_state(self, signed_metrics=True, need_fisher=True):
        # Unit A
        self.unit_A.N_val = self.w_A_N.value
        self.unit_A.s_val = self.w_A_s.value
        self.unit_A.b_val = self.w_A_b.value
        self.unit_A.r_val = self.w_A_r.value
        self.unit_A.traj_mode = self.w_A_mode.value
        self.unit_A.recalculate(signed_metrics=signed_metrics, need_fisher=need_fisher)

        # Unit B
        self.unit_B.N_val = self.w_B_N.value
        self.unit_B.s_val = self.w_B_s.value
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
                        f"N={int(m['N']):,}, s={m['s']:.2f}\n"
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
            ax.set_title(f"Collexeme Dynamics (N, s Fixed / Varying r or b)")
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
                self.w_A_N,
                self.w_A_s,
                self.w_A_b,
                self.w_A_r,
                self.w_A_mode,
            ],
            layout=box_layout,
        )

        box_B = widgets.VBox(
            [
                widgets.HTML("<b>Data B (Red)</b>"),
                self.w_B_N,
                self.w_B_s,
                self.w_B_b,
                self.w_B_r,
                self.w_B_mode,
            ],
            layout=box_layout,
        )

        ui = widgets.VBox(
            [
                widgets.HTML("<h3>Collostructional Simulator (N-s-b-r Model)</h3>"),
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
# sim = CollostructionalComparisonSimulatorNsbr()
# sim.display()

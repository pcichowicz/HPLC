# Main logic (Chromatogram class, etc.)
import pandas as pd
from pandas.core.array_algos.transforms import shift
from pandas.core.frame import DataFrame
from typing import Dict
import numpy as np
import tqdm
import scipy.signal
import scipy.optimize
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import warnings
import seaborn as sns
from datetime import datetime

from . import helpers

# Create Class object for Chromatogram

class Chromatogram:
    """
    Basic class for quantifying and processing of High-Performance Liquid Chromatograph

    Attributes
    ----------
    df : `pandas.core.frame.Dataframe`
        A pandas dataframe containing the chromatogram data, minimum columns
        with time and signal intensity
    peaks : `pandas.core.frame.Dataframe`
        A pandas dataframe containing the inferred properties;
        retention time, scale, skew, amplitude, and total area of
        individual peaks of the chromatogram.
    chromo_matrix : `numpy.ndarray`
        Matrix array where each row corresponds to time, each column to
        a probability density value for each individual peak. Used in the
        `plot` method.
    """

    def __init__(
            self,
            file: str | DataFrame,
            crop_window: bool = None,
            cols: Dict[str, str] = {"time": "time", "signal": "signal"},
    ) -> None:
        """
        Parameters
        ----------

        :param file: `str` or `pandas.core.frame.Dataframe`
            Path to csv file  or pandas dataframe of chromatogram to
            analyze. If None, a pandas Dataframe must be passed
            pandas Dataframe.
        dataframe: `pandas.core.frame.Dataframe`
            pandas Dataframe of chromatogram to analyze, if None, then
            a path to csv file must be passed.
        :param crop_window: `list` [start, end], optional
            The rentention time crop window of the chromatogram to be
            analyzed. If None, then the whole time range is used instead.
        :param cols:
        """

        # Check if file is dataframe
        if type(file) is not pd.core.frame.DataFrame:
            raise RuntimeError(
                f"Argument must be a Pandas Dataframe, given file is of type {type(file)}"
            )

        # Assign class variables
        self.time_col = cols["time"]
        self.signal_col = cols["signal"]

        # Load chromotagram and other components needed.
        chromatogram_dataframe = file.copy()
        self.df = chromatogram_dataframe

        # The average timestep in the chromatogram, ideally should be identical to
        # the values in the chromatogram data. Determine decimal place to ensure
        # there is no float-point precision issues

        self._timestep = np.mean(np.diff(chromatogram_dataframe[self.time_col].values))
        self._timestep_precision = int(np.abs(np.ceil(np.log10(self._timestep))))

        # Define variables that are used by other methods/functions

        self._baseline_corrected = None # removes neg values and subtracts noise
        self._crop_offset = 0

        self.peaks = None # df for peak properties
        self.known_peaks = None
        self._peak_indices = None

        if crop_window is not None:
            self.crop(crop_window)
        else:
            self.df = chromatogram_dataframe

    def __repr__(self):
        trange = f"(t: {self.df[self.time_col].values[0]} - {self.df[self.time_col].values[-1]})"

        rep = f"""                  Chromatogram REPORT                 """
        if self._crop_offset > 0:
            rep += f"\n\t Cropped {trange} \t\t\t\u2713"
        if self._baseline_corrected:
            rep += f"\n\t Baseline Subtracted\t\t\t\u2713"
        if self._peak_indices is not None:
            rep += f"\n\t {len(self._peak_indices)} Peak(s) Detected\t\t\t\u2713\n"
            # Use iloc for positional indexing
            peak_times = self.df.iloc[self._peak_indices][self.time_col].values
            rep += "\t\t\t→ Peak times: " + ", ".join(f"{t:.3f}" for t in peak_times)

        rep = "                  Chromatogram REPORT                 "
        if self._crop_offset > 0:
            rep += "\n\t{:<40} {:>5}".format(f"Cropped {trange}", "\u2713")
        if self._baseline_corrected:
            rep += "\n\t{:<40} {:>5}".format("Baseline Subtracted", "\u2713")
        if self._peak_indices is not None:
            rep += "\n\t{:<40} {:>5}".format(f"{len(self._peak_indices)} Peak(s) Detected", "\u2713")
            peak_times = self.df.iloc[self._peak_indices][self.time_col].values
            rep += "\n\t\t\t→ Peak times: " + ", ".join(f"{t:.3f}" for t in peak_times)

        with open("chromatogram_report.log", "a") as log_file:

            log_file.write(f"\n--- {datetime.now()} ---\n")
            # log_file.write(rep)

        return rep

    def crop(self,
            crop_window: list[float] | None,
            return_df: bool = False
    ) -> None | DataFrame:
        R"""

        :param crop_window:
        :param return_df:
        :return:
        """
        if self.peaks is not None:
            raise RuntimeError(
                """
                you are trying to crop when peaks are fit already.
                """
            )

        if len(crop_window) != 2:
            raise RuntimeError(
                f"`crop_window` must be of len 2, (start, end). Provided list is of length {len(crop_window)}."
            )

        if crop_window[0] >= crop_window[1]:
            raise RuntimeError(
                f"First index (start) is larger than second index (end)"
            )

        # Apply the crop
        self.df = self.df[
            (self.df[self.time_col] >= crop_window[0])
        & (self.df[self.time_col] <= crop_window[1])
        ]
        self._crop_offset = int(crop_window[0] / self._timestep)

        if return_df:
            return self.df

    def _assign_windows(self,
                         known_peaks=[],
                         tolerance=0.5,
                         prominence=0.01,
                         rel_height=1,
                         buffer=0,
                        peak_kwargs={}
    ) -> pd.DataFrame:

        if not (0 <= rel_height <= 1):
            raise ValueError("`rel_height` must be in [0, 1]")

        intensity = self.df[self.signal_col].values
        self.normint = helpers.normalize_signal(intensity)
        self._peak_indices = helpers._detect_peaks(self.normint, prominence, peak_kwargs)

        amps = np.sign(intensity[self._peak_indices])
        _widths = np.zeros_like(amps)
        _left = np.zeros_like(amps, dtype=int)
        _right = np.zeros_like(amps, dtype=int)

        for sign_mask, polarity in zip([amps > 0, amps < 0], [1, -1]):
            if np.any(sign_mask):
                sig = intensity if polarity == 1 else -intensity
                w_half, _, _ = helpers.calculate_peak_widths(sig, self._peak_indices, rel_height=0.5)
                w_full, l_ips, r_ips = helpers.calculate_peak_widths(sig, self._peak_indices, rel_height=rel_height)
                _widths[sign_mask] = w_half[sign_mask]
                _left[sign_mask] = l_ips[sign_mask]
                _right[sign_mask] = r_ips[sign_mask]

        ranges = helpers.build_peak_ranges(_left, _right, len(self.normint), buffer)
        ranges = helpers.remove_subset_ranges(ranges)
        self.ranges = ranges

        window_df = helpers.build_window_df(self, self.df, ranges)
        window_df = helpers.assign_background_windows(window_df)

        self.window_df = window_df
        self.window_props = helpers.extract_window_props(self, window_df, _widths)
        return window_df


    def deconvolve_peaks(
            self,
            verbose=True,
            param_bounds={},
            integration_window=[],
            max_iter=1_000_000,
            optimizer_kwargs={},
    ):
        if self.window_props is None:
            raise RuntimeError("Run `_assign_windows()` first.")

        iterator = (
            tqdm.tqdm(self.window_props.items(), desc="Deconvolving mixture")
            if verbose else self.window_props.items()
        )

        param_order = ["amplitude", "location", "scale", "skew"]
        t_range = helpers._generate_time_range(self.df, self.time_col, integration_window, self._timestep)

        peak_props = {}
        self._param_bounds = []
        self._p0 = []

        for k, v in iterator:
            if v["num_peaks"] == 0:
                continue

            window_dict = {}
            p0 = []
            bounds_lower, bounds_upper = [], []

            if v["num_peaks"] >= 10:
                warnings.warn(
                    f"\nToo many peaks ({v['num_peaks']}) between {v['time_range'].min()}–{v['time_range'].max()}.\n"
                    "This may take a while to finish... "
                )

            for i in range(v["num_peaks"]):
                peak_p0 = [v["amplitude"][i], v["location"][i], v["width"][i] / 2, 0]
                default_bounds = helpers._default_param_bounds(*peak_p0[:3], v["time_range"].min(), v["time_range"].max())

                # Optional tweaks from `param_bounds`
                if param_bounds:
                    custom = {
                        key: param_bounds.get(key, default_bounds[key])
                        for key in param_order
                    }
                    bounds = _adjust_param_bounds(peak_p0, custom, default_bounds, param_order)
                else:
                    bounds = {"lower": [default_bounds[k][0] for k in param_order],
                              "upper": [default_bounds[k][1] for k in param_order]}

                p0.extend(peak_p0)
                bounds_lower.extend(bounds["lower"])
                bounds_upper.extend(bounds["upper"])

            self._p0.append(p0)
            self._param_bounds.append((bounds_lower, bounds_upper))

            # Fit curves
            popt, _ = scipy.optimize.curve_fit(
                helpers._sum_skewnorms,
                v["time_range"],
                v["signal"],
                p0=p0,
                bounds=(bounds_lower, bounds_upper),
                maxfev=max_iter,
                **optimizer_kwargs
            )

            popt = np.reshape(popt, (v["num_peaks"], 4))
            for i, p in enumerate(popt):
                recon_signal = helpers._compute_skewnorm(t_range, *p)
                window_dict[f"peak_{i + 1}"] = {
                    "amplitude": p[0],
                    "retention_time": np.round(p[1], decimals=self._timestep_precision),
                    "scale": p[2],
                    "alpha": p[3],
                    "area": recon_signal.sum(),
                    "reconstructed_signal": recon_signal,
                    "signal_max": np.max(recon_signal),
                }

            peak_props[k] = window_dict

        self._peak_props = peak_props
        return peak_props


    def fit_peaks(
            self,
            tolerance: float = 0.5,
            prominence: float = 1e-2,
            rel_height: float = 1,
            approx_peak_width: float = 5,
            buffer: int = 0,
            param_bounds: Dict[str, list] = {},
            integration_window: list[float] = [],
            verbose: bool = True,
            return_peaks: bool = True,
            correct_baseline: bool = True,
            max_iter: int = 1000000,
            precision: int = 9,
            peak_kwargs: Dict = {},
            optimizer_kwargs: Dict = {},
    ) -> DataFrame:

        if correct_baseline and not self._baseline_corrected:
            self.correct_baseline(
                window=approx_peak_width,
                verbose=verbose,
                return_df=False
            )

        # Assign peak windows
        _ = self._assign_windows(
            tolerance=tolerance,
            prominence=prominence,
            rel_height=rel_height,
            buffer=buffer,
            peak_kwargs=peak_kwargs,
        )

        # Fit skew-normal peaks
        peak_props = self.deconvolve_peaks(
            verbose=verbose,
            param_bounds=param_bounds,
            max_iter=max_iter,
            integration_window=integration_window,
            **optimizer_kwargs,
        )

        # Build dataframe from fitted parameters
        rows = [
            {
                "retention_time": p["retention_time"],
                "scale": p["scale"],
                "skew": p["alpha"],
                "amplitude": p["amplitude"],
                "area": p["area"],
                "signal_maximum": p["signal_max"],
            }
            for window in peak_props.values()
            for p in window.values()
        ]
        peak_df = pd.DataFrame(rows).sort_values(by="retention_time")
        peak_df["peak_id"] = np.arange(1, len(peak_df) + 1).astype(int)
        self.peaks = peak_df

        # Reconstruct unmixed chromatogram matrix
        time = self.df[self.time_col].values
        out = np.zeros((len(time), len(peak_df)))
        i = 0
        for _, peaks in self._peak_props.items():
            for _, v in peaks.items():
                params = [v["amplitude"], v["retention_time"], v["scale"], v["alpha"]]
                out[:, i] = helpers._compute_skewnorm(time, *params)
                i += 1
        self.unmixed_chromatograms = np.round(out, decimals=precision)

        return peak_df if return_peaks else None

    def correct_baseline(
            self,
            window: float = 5,
            return_df: bool = False,
            verbose: bool = True,
            precision: int = 9,
    ) -> DataFrame | None:
        if self._baseline_corrected:
            warnings.warn(
                "Baseline has already been corrected. Rerunning on original signal..."
            )
            self.int_col = self.int_col.split("_corrected")[0]

        if (window / self._timestep) < 10:
            raise ValueError(
                f"""
    The approximate peak width ({window}) is too small relative to the time sampling interval ({self._timestep}).
    Either increase the width or set `correct_baseline=False` to skip this step.
    """
            )

        df = self.df
        signal = df[self.signal_col].copy()

        # Warning if signal has significant negative values, something wrong
        min_val, max_val = np.min(signal), np.max(signal)
        if min_val < 0 and (np.abs(min_val) / max_val) >= 0.1:
            warnings.warn(
                "\x1b[30m\x1b[43m\x1b[1m\n"
                "The chromatogram has significant negative signals. Background noise\n"
                "subtraction may not work as expected. Check results visually to deteremine if correct.\n"
                "\x1b[0m"
            )

        # Shift and clean/mask negative values
        shift = np.median(signal[signal < 0]) if (signal < 0).any() else 0
        signal -= shift
        signal *= np.heaviside(signal, 0)

        # Apply SNIP transformation (LLS)
        tform = np.log(np.log(np.sqrt(signal.values + 1) + 1) + 1)
        n_iter = int(((window / self._timestep) - 1) / 2)

        if verbose:
            self._bg_correction_progress_state = 1
            loop = tqdm.tqdm(range(1, n_iter + 1), desc="Performing baseline correction")
        else:
            self._bg_correction_progress_state = 0
            loop = range(1, n_iter + 1)

        for i in loop:
            tform_new = tform.copy()
            for j in range(i, len(tform) - i):
                tform_new[j] = min(tform[j], 0.5 * (tform[j + i] + tform[j - i]))
            tform = tform_new

        # Inverse transformation of LLS and subtraction
        inv_tform = (np.exp(np.exp(tform) - 1) - 1) ** 2 - 1
        corrected = np.round(signal - inv_tform, decimals=precision)
        estimated_bg = inv_tform + shift

        # Mark that the column has been corrected
        self.df = df.assign(
            **{
                f"{self.signal_col}_corrected": corrected,
                "estimated_background": estimated_bg,
            }
        )
        self._baseline_corrected = True
        self.signal_col = f"{self.signal_col}_corrected"

        return self.df if return_df else None

    def show(self, time_range: list[float] = []) -> list[Figure, Axes]:
        """
        Displays the chromatogram with mapped peaks and fitted signal.

        Parameters
        ----------
        time_range : list[float], optional
            [lower, upper] bounds for time axis zoom.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ax : matplotlib.axes._axes.Axes
        """
        sns.set()
        fig, ax = plt.subplots()
        time = self.df[self.time_col].values

        # Label setup
        ylabel_base = self.signal_col.split("_corrected")[0]

        ax.set_xlabel(self.time_col)
        ax.set_ylabel(f"{ylabel_base}")


        # Raw chromatogram
        ax.plot(time, self.df[self.signal_col], "k-", label="raw chromatogram")

        # Estimated background
        if "estimated_background" in self.df.columns:
            ax.plot(
                time,
                self.df["estimated_background"],
                color="blue",
                label="estimated background",
                zorder=10,
            )

        # Plot fitted peaks if available
        if self.peaks is not None:
            convolved = np.sum(self.unmixed_chromatograms, axis=1)
            ax.plot(time, convolved, "r--", label="inferred mixture")

            for _, peak in self.peaks.iterrows():
                peak_id = int(peak["peak_id"])
                peak_label = self._get_peak_label(peak_id)
                ax.fill_between(
                    time,
                    self.unmixed_chromatograms[:, peak_id - 1],
                    label=peak_label,
                    alpha=0.5,
                )

        # Zoom in on time range/ selected cropped region
        if len(time_range) == 2:
            ax.set_xlim(time_range)
            yvals = self.df[
                (self.df[self.time_col] >= time_range[0]) &
                (self.df[self.time_col] <= time_range[1])
                ][self.int_col]
            ax.set_ylim(ax.get_ylim()[0], 1.1 * yvals.max())

        ax.legend() #bbox_to_anchor=(1.5, 1) Can set to static position
        fig.patch.set_facecolor((0, 0, 0, 0))
        return [fig, ax]

    def _get_peak_label(self, peak_id: int) -> str:
        """
        Generate a label for the given peak_id.

        Parameters
        ----------
        peak_id : int
            The ID of the peak.

        Returns
        -------
        label : str
            Formatted label for plotting.
        """
        return f"peak {int(peak_id)}"

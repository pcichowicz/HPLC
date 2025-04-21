# Modular helper functions like compute_skewnorm
from typing import List, Dict, Tuple, Union
import numpy as np
import pandas as pd
import scipy.signal
import warnings

def normalize_signal(intensity: np.ndarray) -> np.ndarray:
    int_sign = np.sign(intensity)
    norm = (intensity - intensity.min()) / (intensity.max() - intensity.min())
    return int_sign * norm

def _detect_peaks(signal: np.ndarray, prominence: float, peak_kwargs: Dict) -> np.ndarray:
    peaks, _ = scipy.signal.find_peaks(signal, prominence=prominence, **peak_kwargs)
    return peaks

def calculate_peak_widths(intensity: np.ndarray, peak_indices: np.ndarray, rel_height: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=scipy.signal._peak_finding_utils.PeakPropertyWarning)
        widths, _, left_ips, right_ips = scipy.signal.peak_widths(intensity, peak_indices, rel_height=rel_height)
    return widths, left_ips.astype(int), right_ips.astype(int)

def enforce_known_peaks(self, known_peaks: Union[List, Dict], tolerance: float, _widths, _left, _right):
    """Insert known peaks, override auto-detected ones if within tolerance."""
    if isinstance(known_peaks, dict):
        enforced_times = list(known_peaks.keys())
    else:
        enforced_times = known_peaks

    enforced_inds = (np.array(enforced_times) / self._dt).astype(int) - self._crop_offset
    updated_times = np.round(self._dt * (enforced_inds + self._crop_offset), decimals=self._time_precision)

    if isinstance(known_peaks, dict):
        updated_known_peaks = {new: known_peaks[old] for new, old in zip(updated_times, enforced_times)}
    else:
        updated_known_peaks = updated_times.tolist()

    # Remove auto-detected peaks too close to enforced ones
    for loc in enforced_inds:
        close_inds = np.where(np.abs(self._peak_indices - loc) <= (tolerance / self._dt))[0]
        if close_inds.size:
            idx = close_inds[0]
            self._peak_indices = np.delete(self._peak_indices, idx)
            _widths = np.delete(_widths, idx)
            _left = np.delete(_left, idx)
            _right = np.delete(_right, idx)

    # Add enforced peaks
    for i, loc in enumerate(enforced_inds):
        self._peak_indices = np.append(self._peak_indices, loc)
        self._added_peaks = getattr(self, "_added_peaks", [])
        self._added_peaks.append((loc + self._crop_offset) * self._dt)

        peak_width = 1 / self._dt
        if isinstance(known_peaks, dict):
            props = updated_known_peaks[updated_times[i]]
            if isinstance(props, dict) and "width" in props:
                peak_width = props["width"] / self._dt

        _widths = np.append(_widths, peak_width)
        _left = np.append(_left, loc - peak_width)
        _right = np.append(_right, loc + peak_width)

    self._known_peaks = updated_known_peaks
    return _widths, _left, _right

def build_peak_ranges(_left, _right, norm_int_len: int, buffer: int) -> List[np.ndarray]:
    ranges = []
    for l, r in zip(_left, _right):
        rnge = np.arange(int(l - buffer), int(r + buffer))
        rnge = rnge[(rnge >= 0) & (rnge < norm_int_len)]
        ranges.append(rnge)
    return ranges

def remove_subset_ranges(ranges: List[np.ndarray]) -> List[np.ndarray]:
    valid = [True] * len(ranges)
    for i, r1 in enumerate(ranges):
        for j, r2 in enumerate(ranges):
            if i != j and set(r2).issubset(r1):
                valid[j] = False
    return [r for i, r in enumerate(ranges) if valid[i]]

def build_window_df(self, df: pd.DataFrame, ranges: List[np.ndarray]) -> pd.DataFrame:
    df = df.copy()
    df["time_idx"] = np.arange(len(df))
    df["window_id"] = 0
    df["window_type"] = "peak"

    for i, r in enumerate(ranges):
        df.loc[df["time_idx"].isin(r), "window_id"] = i + 1

    return df

def assign_background_windows(window_df: pd.DataFrame) -> pd.DataFrame:
    bg = window_df[window_df["window_id"] == 0]
    tidx = bg["time_idx"].values

    if not len(bg):
        return window_df

    diff = np.diff(tidx)
    split_inds = np.where(diff > 1)[0]

    if len(split_inds) == 0:
        window_df.loc[bg.index, ["window_id", "window_type"]] = [1, "interpeak"]
        return window_df

    split_inds = np.insert(split_inds + 1, 0, 0)
    split_inds = np.append(split_inds, len(tidx))

    for i, (start, end) in enumerate(zip(split_inds[:-1], split_inds[1:])):
        segment = tidx[start:end]
        if len(segment) >= 10:
            window_df.loc[window_df["time_idx"].isin(segment), "window_id"] = i + 1
            window_df.loc[window_df["time_idx"].isin(segment), "window_type"] = "interpeak"

    return window_df[window_df["window_id"] > 0]

def extract_window_props(self, window_df: pd.DataFrame, _widths) -> Dict:
    window_dict = {}
    for gid, group in window_df[window_df["window_type"] == "peak"].groupby("window_id"):
        peak_idxs = [i for i in self._peak_indices if i in group["time_idx"].values]
        peak_inds = [np.where(self._peak_indices == i)[0][0] for i in peak_idxs]
        window_dict[gid] = {
            "time_range": group[self.time_col].values,
            "signal": group[self.signal_col].values,
            "signal_area": group[self.signal_col].sum(),
            "num_peaks": len(peak_idxs),
            "amplitude": [group[group["time_idx"] == p][self.signal_col].iloc[0] for p in peak_idxs],
            "location": [np.round(group[group["time_idx"] == p][self.time_col].iloc[0], self._timestep_precision) for p in peak_idxs],
            "width": [_widths[i] * self._timestep for i in peak_inds],
        }
    return window_dict


def _compute_skewnorm(x, amplitude, loc, scale, alpha):
    _x = alpha * (x - loc) / scale
    norm = (1 / np.sqrt(2 * np.pi * scale**2)) * np.exp(-((x - loc) ** 2) / (2 * scale**2))
    cdf = 0.5 * (1 + scipy.special.erf(_x / np.sqrt(2)))
    return amplitude * 2 * norm * cdf

def _sum_skewnorms(x, *params):
    """
    Sum of skew-normal distributions for curve fitting.
    Each peak is represented by 4 parameters: amplitude, center, width, skew.
    """
    from scipy.stats import skewnorm
    n_params_per_peak = 4
    n_peaks = len(params) // n_params_per_peak
    y = np.zeros_like(x)

    for i in range(n_peaks):
        a, loc, scale, skew = params[i * 4:(i + 1) * 4]
        y += skewnorm.pdf(x, skew, loc, scale) * a

    return y


# def _sum_skewnorms(x, params):
#     params = np.reshape(params, (-1, 4))
#     return sum(_compute_skewnorm(x, *p) for p in params)


def _default_param_bounds(amplitude, location, width, time_min, time_max):
    return {
        "amplitude": np.sort([0.01 * amplitude, 100 * amplitude]),
        "location": [time_min, time_max],
        "scale": [0, (time_max - time_min) / 2],
        "skew": [-np.inf, np.inf],
    }


def _adjust_param_bounds(p0, custom_bounds, default_bounds, param_order):
    param_idx = {k: i for i, k in enumerate(param_order)}
    bounds = {"lower": [], "upper": []}
    for p in param_order:
        if p in custom_bounds:
            if p == "amplitude":
                default_bounds[p] = np.sort(custom_bounds[p])
            elif p == "location":
                default_bounds[p] = [p0[param_idx[p]] + delta for delta in custom_bounds[p]]
            else:
                if custom_bounds[p][0] <= p0[param_idx[p]] <= custom_bounds[p][1]:
                    default_bounds[p] = custom_bounds[p]
                else:
                    raise ValueError(f"Initial guess for '{p}' is outside bounds.")
        bounds["lower"].append(default_bounds[p][0])
        bounds["upper"].append(default_bounds[p][1])
    return bounds


def _generate_time_range(df, time_col, integration_window, timestep):
    if not integration_window:
        return df[time_col].values
    if len(integration_window) == 2:
        return np.arange(integration_window[0], integration_window[1], timestep)
    raise RuntimeError("Integration window must be empty or [start, stop].")

def _get_peak_label(self, peak_id: int) -> str:
    """
    Generate the appropriate label for a peak, using compound mapping if available.

    Parameters
    ----------
    peak_id : int
        The peak ID.

    Returns
    -------
    str
        Label for the peak.
    """
    if self._mapped_peaks is None or peak_id not in self._mapped_peaks:
        return f"peak {peak_id}"

    label = self._mapped_peaks[peak_id]
    try:
        d = self.quantified_peaks[self.quantified_peaks["compound"] == label]
        if "concentration" in d and not d["concentration"].isnull().all():
            label += f"\n[{d['concentration'].values[0]:.3g}"
            if "unit" in d:
                label += f" {d['unit'].values[0]}]"
            else:
                label += "]"
    except Exception:
        pass

    return label

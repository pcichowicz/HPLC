<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>








Subtracting this baseline estimate ($S'$) from the original signal ($S$) yields the "true signal" of the chromatogram. In the HPLC package, this process is wrapped in the .correct_baseline() method:


# Step 2 - Peak Detection

Notebook Code: [![License: MIT](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) Notebook Prose: [![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Baseline correction

Ideally, running a blank HPLC sample should result in a steady, zero or near-zero absorbance throughout the run. However, in reality, many factors disrupt this baseline. These disturbances introduce noise, which can impact the final analysis. Common sources include:

 -  Maintenance of the HPLC valves, flow cells, and UV detectors
 -  Temperature fluctuations of the laboratory during the day
 -  Mobile phase mixing
 -  Incomplete mobile phase degassing

To infer the baseline signal and subtract it from the observed signal, we apply a series of transformations:

- Log-transformation of the signal
- Iterative Minimum Filtering
- Inverse Transformation and Subtraction

### Log-transformation 

The first step in determining the baseline involves reducing the dynamic range of the signal data *S* using a _Linear Least Squares_ (LLS) approximation. This transformation prevents large peaks from dominating the filtering which can lead to the removal of smaller and still important peaks. The compression  $S \rightarrow S_{LLS}$ is achieved mathematically using 

$$
S_{LLS} = \ln \left[ \ln \left( \sqrt{S + 1} + 1 \right) + 1 \right]
$$

where the square root operation enhances smaller peaks, while the log operator compresses larger peaks — scaling the signal to a narrower, more balanced range.

### Iterative Minimum Filtering

Since peaks are considered "postive deviations" compared to the baseline, the local minimum would reflect the baseline. The iteration "erodes" the peaks gradually until only the baseline is left. This process is repeated several times until a smooth, baseline-only signal is obtained.

$$
S_{LLS_{m}}^{\prime} (t) = min \left[ S_{LLS_{m-1}} (t), \frac{S_{LLS_{m-1}} (t-m) + S_{LLS_{m-1}} (t+m)}{2} \right]
$$

### Inverse transformation and subtraction

Once we have the signal filtered and smoothed down to reflect the baseline, the $S_{LLS}$ can be then passed through the inverse LLS transformation to expand the range back to the scale of the original observed data.

$$
S' = \left( \exp \left[ \exp \left( S'_{LLS} \right) -1 \right\] -1 \right)^2 -1
$$

Subtraction of S - S' will remove the baseline signal and we are left with the "true signal" of the chromatogram. In the `HPLC` package, there is a method called `.correct_baseline()` that will perform all of these transformations.

```python
 chromatogram.correct_baseline()

 chromatogram.show()
```


## Detection of peaks

After the background noise is dealt with, the main functions of the `HPLC` package can perform their jobs. By calling `.fit_peaks()`, the deconvolution of the singal to separate peaks (and parameters/properties) is performed internally by the function. This involves the baseline correction, detecting the peaks, fitting of the peaks (using skew-normal distributions) and storing the properies of each peak and its respective window/ peak ID.

```python
  chromatogram.fit_peaks()

  chromatogram.show()
```

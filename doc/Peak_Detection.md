<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

# Step 2 - Peak Detection

Notebook Code: [![License: MIT](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) Notebook Prose: [![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Baseline correction

Ideally running a blank HPLC should result in a steady zero/near zero abosrbance during the entire HPLC run, but in reality there are many factors that disrupt this baseline. These factors cause "noise" that can impact the final analysis of HPLCs, such as but not limited to:
 -  Maintenance of the HPLC valves, flow cells, and UV detectors
 -  Temperature fluctuations of the laboratory during the day
 -  Mobile phase mixing
 -  Mobile phase degassing

To infer what the baseline signal within the data and subtract it from our observed signal, we need to employ few transformations;

- Log-transformation of the signal
- Iterative Minimum Filtering
- Inverse Transformation and Subtraction

### Log-transformation 

The first step in determining the baseline involves reducing the dynamic range of the signal data *S* using a _Linear Least Squares_ (LLS) approximation. This transformation prevents the large peaks from dominating the filtering which leads to the removal of smaller and still important peaks. The compression  $S \rightarrow S_{LLS}$ is achieved mathematically using 

$$
S_{LLS} = \ln \left[ \ln \left( \sqrt{S + 1} + 1 \right) + 1 \right]
$$

where the square root operation enhances the smaller peaks, and the log operator comrpesses the large peaks, scaling down to better match smaller peaks.

### Iterative Minimum Filtering

Since peaks are considered "postive deviations" compared to the baseline, the local minimum would reflect the baseline. The iteration "erodes" the peaks gradually until the baseline is left. This process repeatedly smoothes and filters the signal and keep the local minimum with each iteration.

$$
S_{LLS_{m}}^{\prime} (t) = min \left[ S_{LLS_{m-1}} (t), \frac{S_{LLS_{m-1}} (t-m) + S_{LLS_{m-1}} (t+m)}{2} \right]
$$

### Inverse transformation and subtraction

Once we have the signal filtered and smoothed down to a minimum value, the $S_{LLS}$ can be then passed through the inverse LLS operation to expand the range back to the scale of the original observed data.

$$
S' = \left( \exp \left[ \exp \left( S'_{LLS} \right) -1 \right\] -1 \right)^2 -1
$$

Subtraction of S - S' will remove the baseline signal and we are left with the "true signal" of the chromatogram. To do this in the `HPLC` package, there is a method called `.correct_baseline()` that will perform all of these transformations.

```python
 chromatogram.correct_baseline()

 chromatogram.show()
```

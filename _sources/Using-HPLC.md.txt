# High-Performace Liquid Chromatography

High-performance liquid chromatography (HPLC) is an analytical technique used to separate, identify, and quantify components in a chemical mixture. This technique employs high-pressure pumps to deliver the sample mixture, dissolved in a solvent (known as the mobile phase), through a tightly packed column filled with solid absorbent particles (the stationary phase). The different components in the mixture interact uniquely with the stationary phase, causing them to travel through the column at varying rates. This difference in migration allows for the separation of compounds.

As the separated components exit the column, a detector measures their presence and generates a signal. The output is visualized as a graph called a chromatogram, which plots signal intensity versus time. Each peak on the chromatogram corresponds to a distinct compound in the sample. The position of a peak, known as the retention time, reflects the time it took for that compound to travel through the column. Additionally, the area under each peak is proportional to the concentration of the corresponding compound in the mixture.


## Manual to Automated

There are many software tools capable of creating simple graphs and performing basic analysis; however, they are often limited in terms of customization. This is where this package—and many others like it—becomes advantageous, as it can be modified and extended to suit specific needs.


How  `HPLC` works

The peak detection and quantification  in HPLC involves the following steps.

1. Estimation of signal background using Statistical Nonlinear Iterative Peak (SNIP) estimation.
2. Using peak maxima with thresholds, automatic peak detection is performed
3. Division of the chromatogram into "peak windows", where each window ideally contains a single peak. If a region contains overlapping peaks, they are grouped into a single window.
4. For windows with multiple peaks, skew-normal distributions are fitted using initial guesses based on peak properties such as retention time, maximum signal, and width at half-maximum.
5. Each distribution’s expected signal is computed using the best-fit parameters, and the integrated area under each peak is calculated and stored.

The remainder of this notebook will walk through each of these steps in detail.

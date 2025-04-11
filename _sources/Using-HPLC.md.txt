# High-Performace Liquid Chromatography

High-performance liquid chromatography (HPLC) is an analytical technique used to separate, identify, and quantify components in a chemical mixture. This technique employs high-pressure pumps to deliver the sample mixture, dissolved in a solvent (known as the mobile phase), through a tightly packed column filled with solid absorbent particles (the stationary phase). The different components in the mixture interact uniquely with the stationary phase, causing them to travel through the column at varying rates. This difference in migration allows for the separation of compounds.

As the separated components exit the column, a detector measures their presence and generates a signal. The output is visualized as a graph called a chromatogram, which plots signal intensity versus time. Each peak on the chromatogram corresponds to a distinct compound in the sample. The position of a peak, known as the retention time, reflects the time it took for that compound to travel through the column. Additionally, the area under each peak is proportional to the concentration of the corresponding compound in the mixture.


## Manual to Automated

There are many softwares that are capable of creating simple graphs and analysis, however they are often limited with customizatiom. This is where this package and many others available can be advantageous and be modified/changedto suit the desired needs.


How  `HPLC` works

The peak detection and quantification  in HPLC involves the following steps.

1. Estimation of signal background using Statistical Nonlinear Iterative Peak (SNIP) estimation.
2. Using peak maxima with thresholds, automatic peak detection is done
3. Chromatogram is divided into "peak windows" where one peak is present, if a window is heavily overlapped with singal then the peaks are grouped into a single winodw.
4. Windoes with N peaks, skew-normal distributions are inferred with peak properties (location, maximum value, and width at half-maximum) as initial guesses. Each window in the chromatogram is inferred in this way.
5. Each distributions expected signal is computed with best-fit parameters, then the integrated signal over the entire peak is computed and stored.

The rest of the notebook will go into more detail for each part.

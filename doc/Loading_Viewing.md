# Step 1 - Loading & viewing the data

Notebook Code: [![License: MIT](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) Notebook Prose: [![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Loading Chromatograms

Part of the `HPLC` package includes the data loading method called `load_chromatogram` that takes in text files of chromatogram data and reads it into a pandas DataFrame. The text file should be similar to something like this:

Metadata line 1  
Metadata line 2  

time,signal  
0.0,0.00037  
0.00833,-0.00135  
0.01666,-0.00244  
0.02499,0.00208  
0.03332,2e-05  
...,...   
19.95868,0.00044  
19.96701,-0.00057  
19.97534,-0.00095  
19.98367,-0.00078  
19.992,0.00062  

Loading the text file is done with one line of command using `hplc.io.load_chromatogram()` that will create the required DataFrame object that will be used with the `Chromatogram` class.

```python
  df = hplc.io.load_chromatogram("data.csv",
                                ['time', 'signal'])
```

| time  |  signal |
| ---- | ------ |
| 0.0  | 0.00037|
| 0.00833 | -0.00135 |
| 0.01666 | -0.00244 |
| 0.02499 | 0.00208 |
| 0.03332 | 2e-05 |
| 0.04165 | -1e-05 |
| ... | ... |
| 19.95868 | 0.00044 |
| 19.96701 | -0.00057 |
| 19.97534 | -0.00095 |
| 19.98367 | -0.00078 |
| 19.992 | 0.00062 |

## Viewing the data

The `Chromatogram` class has methods for viewing, cropping and quantification of the raw chromatogram data. Using the `Chromatogram` class to initialize a new class object, followed by `.show()` on the newly created object you can view the raw chromatogram.

```python
  chromatogram = Chromatogram(df)

  chromatogram.show()
```

![](https://raw.githubusercontent.com/pcichowicz/HPLC/main/doc/plots/Chrom_raw.png)


It is also possible to crop the plot to the desired window of interest for better resolution if there are multiple peaks overlapping. In this instance that is not the case, however lets crop and focus on the third peak in the chromatogram.

```python
  chromatogram.crop([10,20])

  chromatogram.show()
```

![][(https://github.com](https://raw.githubusercontent.com)/pcichowicz/HPLC/blob/main/doc/plots/Chrom_raw_crop.png)


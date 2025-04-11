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

To view the raw chromatogram data, you load the DataFrame into the `Chromatogram` object that has simple methods for viewing and cropping the plot.

```python
  chromatogram = Chromatogram(df)

  chromatogram.show()
```


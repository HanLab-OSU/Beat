# Beat - Base Editing Analysis Tool
Beat determines editing efficiency after subtracting
the background noise without the need to normalize to control samples. 
It finds the noise and filters the outliers from the noise
using the Median Absolute Deviation (MAD) method

## The Executable files
To download the .exe file, please use this link:
```bash
https://drive.google.com/open?id=13NvLL70i7sTlNw04FF_n4b-RxFndlovw
```

To download the macOS executable, please click this link:
```bash
https://github.com/HanLab-OSU/Beat/releases/download/macOS/Beat_OSX
```

## Usage
The usage:
```bash
python Beat.py [--data-path] [--filename|all] [--spacer] [-base-change-position] [--change-pattern] 
```
or for the batch analysis:
```bash
python Beat.py [--csv-path]
```
## Example
To run batch analysis:
```bash
python Beat.py ./data all GAGTATGAGGCATAGACTGC 5 AG
```
to analyze all sequencing files if they have the same spacer.
And if they contain different spacers, just run as
```bash
python Beat.py ./data/template1.csv 
```

To analyze individual sequencing file, run
```bash
python Beat.py ./data Site2_PCR_420.ab1 GAGTATGAGGCATAGACTGC 5 AG
```
## Dependency
To install the necesssary packages, run:
```bash
pip install biopython pandas numpy scipy openpyxl
```

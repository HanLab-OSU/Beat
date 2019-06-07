# Beat - Base Editing Analysis Tool
Beat8 determines editing effiency after substracting
the background noise without normalization to controls. 
It finds the noise and filters the outliers from the noise
using the Median Absolute Deviation (MAD) method

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
You can run the following code to install the necesssary packages:
```bash
pip install biopython pandas numpy scipy openpyxl
```

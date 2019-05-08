# Beat8 - Base Editing Analysis Tool version 8
Beat8 determines editing effiency after substracting
the background noise without normalization to controls. 
It finds the noise and filters the outliers from the noise
using the Median Absolute Deviation (MAD) method

## Usage example
To run batch analysis:
```bash
python Beat8.py foldername all GAGTATGAGGCATAGACTGC 5 AG
```
to analyze all sequencing files if they have the same spacer.
And if they contain different spacers, just run as
```bash
python Beat8.py ./example/template1.csv 
```

To analyze individual sequencing file, run
```bash
python Beat8.py foldername Site2_PCR_420.ab1 GAGTATGAGGCATAGACTGC 5 AG
```

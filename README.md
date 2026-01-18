# Dimuon Scounting analysis

## Running the analysis


#### Merging the histograms
```py
usage: mergeHistos.py [-h] [-i INPUTFILES] [-s SEARCHSTRING] [-t TAG] [-o OUTFILE]

options:
  -h, --help            show this help message and exit
  -i INPUTFILES, --inputFiles INPUTFILES
                        inputFile
  -s SEARCHSTRING, --searchString SEARCHSTRING
                        search string
  -o OUTFILE, --outFile OUTFILE
                        Output json file

## @ will be replaced by * 
python3 python/mergeHistos.py  -s results/analysis/v0/MinBias/parts/out_data_highptMuMu_MinBias_@.json -o results/analysis/v0/MinBias/out_data_highptMuMu_MinBias.json
```

### Printing the summary


```bash
usage: printSummary.py [-h] [-i INPUTFILE] [-l LUMI] [-x XS] [-n HNAME] [-p PROC] [--xmin XMIN] [--xmax XMAX] [-o OUTFILE]
                       [--printAllHistograms]

options:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        inputFile
  -l LUMI, --lumi LUMI  Target luminosity
  -x XS, --xs XS        Input cross-section
  -n HNAME, --hname HNAME
                        hitogram to check selections
  -p PROC, --proc PROC  Process name to get xs from
  --xmin XMIN           min bound in histogram
  --xmax XMAX           max bound in histogram
  -o OUTFILE, --outFile OUTFILE
                        Output fileout
  --printAllHistograms  print the details of all the available histograms

## Example
python3 python/printSummary.py -i results/analysis/v0/MinBias/out_data_highptMuMu_MinBias.json -p MinBias

```

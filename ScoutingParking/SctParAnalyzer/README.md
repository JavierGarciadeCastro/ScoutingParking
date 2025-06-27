# nTuplizer and Analyzer

This code currently takes 2024 central samples as inputs and outputs a flat TTree of variables to facilitate the plotting process. For the moment the data used is only from scouting, and the plotting is very preliminar (also there is a bug in the muon and SV ttrees, they are not being filled correctly).

## HOW TO RUN:
1. Set up a CMSSW release:
```bash
cmsrel CMSSW_15_0_2
cd CMSSW_15_0_2/src
cmsenv
```
2. Clone this repository in you local machine:
```bash
git clone https://github.com/JavierGarciadeCastro/ScoutingParking.git
```
3. Run the Ntuplizer:
```bash
cd ScoutingParking/SctParAnalyzer
cmsRun test/run_dark_shower.py
```
4. This should output a root file called output.root, optionally you can plot certain variables written in this root file with:
```bash
python3 make_plots.py
```

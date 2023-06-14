# G4EPP
A Python wrapper for processing and plotting GEANT4 energetic particle precipitation (EPP) model outputs.

To install, run 

```python 
python pip install -m G4EPP
```
## Notes
- The data requires about 0.5 GB free disk space.

The data directory doesn't get cleaned correctly if you run 
```python
python pip uninstall G4EPP
```
and you'll need to remove it by hand. You can do that with a line like

rm -rf ~/anaconda3/envs/py3/lib/python3.7/site-packages/G4EPP/data

if you're using anaconda, or you can run

```bash
find ~/ -name G4data_mono_discretePAD_0degLat.pkl
```

then remove that file. Email me at grant.berland@colorado.edu if you have any questions!


# G4EPP
A Python wrapper for processing and plotting GEANT4 energetic particle precipitation (EPP) model outputs.

To install, run 

```python 
python pip install -m G4EPP

## Example Usage

The following examples are simple usages of G4EPP with Juypter Notebooks. 

[Example 1](G4EPP/G4EPP/examples/Ionization_profile_examples.ipynb) is an example of various input electron spectra and the resulting ionization profiles they would generate.

[Example 2](G4EPP/G4EPP/examples/X_ray_analysis_examples.ipynb) is a simple X-ray spectral analysis example.



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


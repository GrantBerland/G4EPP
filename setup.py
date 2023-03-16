from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'GEANT4 EPP Simulation wrapper'
LONG_DESCRIPTION = 'Python wrapper to process and plot GEANT4 EPP simulation outputs'

setup(
        name="G4EPP", 
        version=VERSION,
        author="Grant Berland",
        author_email="grant.berland@colorado.edu",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["matplotlib.pyplot", "pandas", "numpy", "scipy", "seaborn", "collections"],  
        keywords=['python', 'EPP', 'UQ'],
        classifiers= []
)

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from collections import OrderedDict

import numpy as np
import pandas as pd
import scipy as sci


import seaborn as sns
sns.set_theme(font_scale=1.5);

from scipy.optimize import curve_fit

from utils import FileReader, EnergyDistributions

# Globals
Earray = np.round(np.logspace(np.log10(0.250), np.log10(10000), 100), 2)
h      = np.linspace(0, 499, 500);
runList = np.array([10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
energyDistList = ["mono", "exp"]
PAlist   = np.arange(0, 70+5, 5);


class api(FileReader):
    def __init__(self):

        # Altitude array abscissa
        self.h = h 
        # Energy array abscissa
        self.Earray = Earray 

        
        self.binWidth = np.hstack([np.diff(self.Earray), np.diff(self.Earray)[-1]])


        self.runList        = runList
        self.energyDistList = energyDistList
        self.PAlist         = PAlist


        super().__init__(self.Earray, self.runList, self.PAlist, self.energyDistList)

        self.X_energy, self.Y_altitude = np.meshgrid(self.Earray, self.h)
        
    def plot_ionization_profile(self, energyDist, energy, flux, pitchAngle, particle=None):

        self._validateInputs(energyDist, energy)

        mean, std = self._readInData(energyDist, energy, "ioni", pitchAngle, particle)

        if energyDist == "exp":
            t = "Exponential"
        elif energyDist == "mono":
            t = "Monoenergetic"
        
        plt.semilogx(flux * mean, self.h)
        plt.fill_betweenx(self.h, flux * (mean - std), 
                                  flux * (mean + std), 
                                  alpha=0.5, 
                                  label='$\pm \sigma$')
    
        plt.xlabel("Deposited Energy [keV km$^{-1}$]");
        plt.ylabel("Altitude [km]");

        plt.ylim(0, 300);
        plt.grid(True, which='both')
        
        plt.title("Ionization Profile\n" + str(energy) + " keV " + t + " Distribution")
        
        return mean

    def plot_spectral_profile(self, energyDist, energy, flux, pitchAngle, particle=None):

        self._validateInputs(energyDist, energy)
        
        sp_mean, sp_std = self._readInData(energyDist, energy, "spectra", pitchAngle, particle)
        
        if energyDist == "exp":
            t = "Exponential"
        elif energyDist == "mono":
            t = "Monoenergetic"
        
        p = plt.pcolormesh(self.X_energy, self.Y_altitude, flux * sp_mean, norm=LogNorm())

        plt.colorbar(p, label='Flux [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')

        plt.xlabel('Energy [keV]')
        plt.ylabel('Altitude [km]')
        
        plt.xscale('log')
        plt.ylim(0, 300)

        plt.grid(True, which='both')

        plt.title("Spectral Profile\n" + str(energy) + " keV " + t + " Distribution")

        return sp_mean 

    def plot_example(self, energyDist, energy, flux, pitchAngle, particle=None):

        self._validateInputs(energyDist, energy)
        
        mean, std       = self._readInData(energyDist, energy, "ioni", pitchAngle, particle)
        sp_mean, sp_std = self._readInData(energyDist, energy, "spectra", pitchAngle, particle)
        
        if energyDist == "exp":
            t = "Exponential"
        elif energyDist == "mono":
            t = "Monoenergetic"
       
        plt.figure(figsize=(14,8))
        plt.subplot(1,2,1)
        plt.semilogx(mean * flux, self.h)
        plt.fill_betweenx(self.h, flux * (mean - std), flux * (mean + std), 
                                  alpha=0.5, 
                                  label='$\pm \sigma$')
    
        plt.xlabel("Deposited Energy [keV km$^{-1}$]");
        plt.ylabel("Altitude [km]");

        plt.ylim(0, 300);
        plt.grid(True, which='both')
        
        plt.title("Ionization Profile\n" + str(energy) + " keV " + t + " Distribution")

        ######

        plt.subplot(1,2,2)
        p = plt.pcolormesh(self.X_energy, self.Y_altitude, flux * sp_mean, norm=LogNorm())

        plt.colorbar(p, label='Flux [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$]')

        plt.xlabel('Energy [keV]')
        
        plt.xscale('log')
        plt.ylim(0, 300)

        plt.grid(True, which='both')

        plt.title("Spectral Profile\n" + str(energy) + " keV " + t + " Distribution")
    
    def plot_all_ionization_profiles(self, energyDist, flux, pitchAngle, particle=None):

        for item in self.runList:
            mean, std = self._readInData(energyDist, item, "ioni", pitchAngle, particle)

            plt.semilogx(flux * mean / item,  self.h, label='%.0f keV' % item)

    
        plt.xlabel("Deposited Energy [keV km$^{-1}$ keV$^{-1}$ / (cm$^{-2}$ s$^{-1}$ sr$^{-1}$)]");
        plt.ylabel("Altitude [km]");
        plt.legend()

        plt.ylim(0, 300);
        plt.grid(True, which='both')
        
    # Getter methods
    def get_ionization_profile(self, energyDist, energy, flux, pitchAngle, particle=None):
      
        self._validateInputs(energyDist, energy)

        return flux * np.array(self._readInData(energyDist, energy, "ioni", pitchAngle, particle))

    def get_all_ionization_profiles(self):
        return self._get_ionization_table()


    def get_spectral_profile(self, energyDist, energy, flux, pitchAngle, particle=None):

        self._validateInputs(energyDist, energy)

        return flux * np.array(self._readInData(energyDist, energy, "spectra", pitchAngle, particle))
    
    def get_altitude_array(self):
        return self.h

    def get_spectral_abscissa(self):
        return self.X_energy, self.Y_altitude

    


    def get_energy_array(self):
        return self.Earray

    def get_bin_width(self):
        return self.binWidth

    def get_run_list(self):
        return self.runList

    def get_PA_list(self):
        return self.PAlist


class XrayAnalysis(FileReader):
    def __init__(self):
        
        # Import globals
        self.Earray         = Earray
        self.runList        = runList
        self.energyDistList = energyDistList

        # Send them up
        super().__init__(self.Earray, self.runList, self.energyDistList)

    def getXraysAtAltitude(self, energyDist, energy, flux, altitude):
        
        mean, std = self._readInData(energyDist, energy, "spectra", "photon")

        # Integrates over a 10 km bin
        spectrum = np.sum(flux * mean[altitude-5:altitude+5, :], axis=0)

        plt.step(Earray, spectrum, label="%.0f km" % altitude)

        plt.yscale('log')
        plt.xscale('log')
        
        plt.xlabel('Energy [keV]')
        plt.ylabel('Flux [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$]')
        
        plt.grid(True, which='both')

        plt.title('X-ray Spectrum')

        return spectrum

    def invertToElectronSpectrum(self, xray_spectrum, altitude=400):
        
         
        electron_spectrum = xray_spectrum

        self.beams = OrderedDict()
        # Get monoenergetic spectral profiles
        for item in self.runList:
            mean, std = self._readInData("mono", item, "spectra", "photon")

            # Integrate at correct altitude and store in dictionary
            self.beams[item] = np.sum(mean[altitude-5:altitude+5, :], axis=0)


        popt, pcov = curve_fit(self._exponentialDistribution, self.Earray, xray_spectrum) 

        # TODO: implement cost-fnc fit method used in /AEPEX/inversion_analysis

        return electron_spectrum

    def _exponentialDistribution(self, x, a, b):
        return a * np.exp(-x / b)



class RadiationAnalysis:
    def __init__(self, material):
        
        if material is "human":
            
            filepath_pre  = "../data/radiationData/"
            filepath_post = "_human_dose_conversion.csv"

            # Read in data
            # [energy in keV, conversion factor in Sv-cm^2]
            self.el_conv = pd.read_csv(filepath_pre + "electron" + filepath_post, names=['E', 'C']) 
            self.ph_conv = pd.read_csv(filepath_pre + "photon"   + filepath_post, names=['E', 'C'])

    def calculate_radiation_dose(self, el_spectra, ph_spectra):
    
        dosage_alt_array = np.zeros([len(h)]);
   
        for i in range(0, len(h)):
            dosage_alt_array[i] = np.sum(el_spectra[i,:] * self.el_conv.C + \
                                         ph_spectra[i,:] * self.ph_conv.C)

        return dosage_alt_array;

    def _calculate_differential_radiation_dose(self, el_spectra, ph_spectra):
    
        dosage_alt_array = np.zeros([len(h), 100]);
    
        for i in range(0, len(h)):
            dosage_alt_array[i,:] = el_spectra[i,:] * self.el_conv.C + \
                                    ph_spectra[i,:] * self.ph_conv.C

        return dosage_alt_array;



'''
TODO: loop in Jovian EPP model
class JUNO(DataHandler):
    def __init__(self):
        pass
'''

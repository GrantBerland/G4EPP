
import numpy as np
import os 

import pickle


class EPP_Exception(Exception):
    def __init__(self, issue, message="Error: "):
        self.issue   = issue
        self.message = message

        super().__init__(self.message)

class EPP_Exception_Handler(EPP_Exception):
    def __init__(self, runList):

        self.runList = runList

    def _validateInputs(self, energyDist, energy):

        if energy not in self.runList:

            raise EPP_Exception(energy, "Error: %.1f not in %s" % (energy, self.runList))


class FileReader(EPP_Exception_Handler):
    def __init__(self, Earray, runList, PAlist):

        self.data_path = os.path.dirname(__file__) + "/data/"

        super().__init__(runList)

        # Load in pkl data table 
        self.D = pickle.load(open(self.data_path + "G4data_mono_discretePAD_0degLat.pkl", "rb"))

        self.runList        = runList
        self.PAlist         = PAlist



    def _get_ionization_table(self):

        table = np.zeros([500, len(self.runList), len(self.PAlist)]);

        for ind1, ene in enumerate(self.runList):
            for ind2, pa in enumerate(self.PAlist):
                table[:, ind1, ind2] = self.D[('electron', 'ioni', ene, pa)][0] + \
                                       self.D[('photon', 'ioni', ene, pa)][0] / 100



        return table

    def _get_all_data(self):
        return self.D

    def _formGreensFunctionSpectrum(self,
                                    energyDistribution,
                                    pitchAngleDistribution,
                                    flux,
                                    dataType,
                                    particle=None):

        testArray = np.hstack([energyDistribution, pitchAngleDistribution])

        if (np.isnan(testArray)).any():
            raise ValueError("Nan(s) in inputs!")

        if (np.isinf(testArray)).any():
            raise ValueError("Inf(s) in inputs!")

        # Energy array in eV for convienience
        energyAbsc = self.runList * 1e3

        # Normalize energy distribution
        energyDistribution     /= np.trapz(energyDistribution, x=energyAbsc)

        # Normalize pitch angle distribution
        pitchAngleDistribution /= np.trapz(pitchAngleDistribution, x=np.deg2rad(self.PAlist))


        # Compute energy flux input 
        # Solid angle calculation
        int1 = 2 * np.pi * np.trapz(pitchAngleDistribution * np.sin(np.deg2rad(self.PAlist)),
                                    x=np.deg2rad(self.PAlist))

        # First moment of energy distribution
        int2 = np.trapz(energyDistribution * energyAbsc,
                        x=energyAbsc)

        norm = flux * int1 * int2

        data = self._get_all_data()

        if dataType == 'ioni':
            result     = np.zeros(500)
            multFactor = 0.035

        elif dataType == "spectra":
            result     = np.zeros([500, 100])
            multFactor = 1

        # TODO: make this a matrix multiplication for SPEED
        # See documentation for explanation of divisive factors in below loop
        for ind1, ene in enumerate(self.runList):
            for ind2, pa in enumerate(self.PAlist):

                weight = energyDistribution[ind1] * pitchAngleDistribution[ind2]

                angFactor = (1e5 *  2 * np.pi * np.cos(np.deg2rad(pa)) * multFactor)

                if particle is None:

                    result += weight * (data[('electron', dataType, ene, pa)][0] + \
                           data[('photon', dataType, ene, pa)][0]/100) / angFactor

                if particle == 'electron':

                    result += weight * data[('electron', dataType, ene, pa)][0] / angFactor

                if particle == 'photon':

                    result += weight * data[('photon', dataType, ene, pa)][0]/100 / angFactor


        # (ioni ~ cm^-3 s^-1, spectra ~ keV cm^-2 s^-1 sr^-1 keV^-1)
        # norm ~ eV / cm^2 / sec
        return result, norm

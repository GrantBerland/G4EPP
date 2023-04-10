import numpy as np
import pandas as pd

import pickle

from scipy.signal import savgol_filter


class EPP_Exception(Exception):
    def __init__(self, issue, message="Error: "):
        self.issue   = issue
        self.message = message

        super().__init__(self.message)

class EPP_Exception_Handler(EPP_Exception):
    def __init__(self, runList, energyDistList):

        self.runList = runList

        self.energyDistList = energyDistList

    def _validateInputs(self, energyDist, energy):

        if energy not in self.runList:
            
            raise EPP_Exception(energy, "Error: %.1f not in %s" % (energy, self.runList))
        
        elif energyDist not in self.energyDistList:

            raise EPP_Exception(energy, "Error: Distribution %s not in %s" % (energyDist, self.energyDistList))



class FileReader(EPP_Exception_Handler):
    def __init__(self, Earray, runList, PAlist, energyDistList):
        
        self.data_path = "../data/"

        super().__init__(runList, energyDistList)

        # Load in pkl data table 
        self.D = pickle.load(open(self.data_path + "G4data_mono_discretePAD_test.pkl", "rb"))


        self.runList        = runList
        self.PAlist         = PAlist
        self.energyDistList = energyDistList


    def _readInData(self, energyDist, energy, dataType, pitchAngle, particle=None):
        
        # Select if electron, photon, or both ionization is desired
        if particle is None:
    
            return self.D[('electron', dataType, energy, pitchAngle)][0] + \
                   self.D[('photon', dataType, energy, pitchAngle)][0], \
                   self.D[('electron', dataType, energy, pitchAngle)][1] + \
                   self.D[('photon', dataType, energy, pitchAngle)][1]

        elif particle is "electron" or particle is "photon":
            
            return self.D[(particle, dataType, energy, pitchAngle)][0], \
                   self.D[(particle, dataType, energy, pitchAngle)][1]
        

        else:
            raise ValueError("%s not an option; please select from 'electron', 'photon', or leave empty" % particle)


    def _get_ionization_table(self):
        
        table = np.zeros([500, len(self.runList), len(self.PAlist)]);

        for ind1, ene in enumerate(self.runList):
            for ind2, pa in enumerate(self.PAlist):
                table[:, ind1, ind2] = self.D[('electron', 'ioni', ene, pa)][0] + \
                                       self.D[('photon', 'ioni', ene, pa)][0] / 100



        return table 

    def _get_all_data(self): 
        return self.D


class EnergyDistributions:
    def __init__(self, Nsamples):
        
        # Modified Bessel function for relativistic Maxwellian
        from scipy.special import kn 
        
        # Distribution PDFs
        self.f_powerLaw       = lambda E, alpha, Emin: (alpha-1) / Emin *  (E / Emin)**-alpha
        
        self.f_exponential    = lambda E, E0: 1/E0 * np.exp(-E / E0)

        self.f_doubMaxwellian = lambda E, T1, T2: np.sqrt(E/np.pi) * (1/T1)**(3/2) * np.exp(-E/T1) + np.sqrt(E/np.pi) * (1/T2)**(3/2) * np.exp(-E/T2) 
       
        m_el = 511 # Electron mass , [keV/c^2]
        self.f_relMaxwell     = lambda E, T: (1 + E/m_el)**2 * np.sqrt(1 - 1/(1 + E/m_el)**2) / ( m_el * (T/m_el) * kn(2, m_el/T) ) * np.exp(-(1 + E/m_el) / (T/m_el))

        # Inverse CDF sampling formulas
        self.powerLawSampler    = lambda U, alpha, Emin: np.power(U, 1/(1-alpha)) * Emin
        self.exponentialSampler = lambda U, E0: -E0 * np.log(U)

        # Save input data
        self.Nsamples = Nsamples
        self.samples  = np.zeros([Nsamples])


    def _doubleMaxwellianSampler(self, T1, T2):

        i = 0
        c = 1
        while i < self.Nsamples:
    
            # Helper distribution (exponential with E0 = (T1 + T2)/2
            Y = -(T1+T2)/2 * np.log(np.random.rand());

            # Accept/reject algorithm
            if np.random.rand() < self.f_doubleMaxwell(Y, T1, T2)/(c * self.f_exponential(Y, (T1 + T1)/2)):
        
                # If loop entered, accept sample point 
                self.samples[i] = Y;
        
                i += 1;

        return self.samples
   
    def _relativisticMaxwellianSampler(self, T):
        
        i = 0
        c = 2
        while i < self.Nsamples:
    
            # Helper distribution (exponential with E0 = T on [1, inf)
            Y = -T * np.log(np.random.rand());

            # Accept/reject algorithm
            if np.random.rand() < self.f_relMaxwell(Y, T)/(c * self.f_exponential(Y, T)):
        
                # If loop entered, accept sample point 
                self.samples[i] = Y;
        
                i += 1;
        
        return self.samples

    def generate_mixture_spectrum(self, Earray, weights, coefficients):
        
        # Power law 
        alpha = coefficients[0]
        Emin  = coefficients[1]

        # Exponential
        E0    = coefficients[2]

        # Double Maxwellian
        T1    = coefficients[3]
        T2    = coefficients[4]

        # Relativistic Maxwellian
        T     = coefficients[5]

        return weights[0] * f_powerlaw(Earray, alpha, Emin) +  \
               weights[1] * f_exponential(Earray, E0) +        \
               weights[2] * f_doubMaxwellian(Earray, T1, T2) + \
               weights[3] * f_relMaxwellian(Earray, T)

    def sample_from_mixture_spectrum(self, weights, coefficients):

        U = np.random.rand(Nsamples, 4)

        return weights[0] * self.powerlawSampler(U[:,0], alpha, Emin) +    \
               weights[1] * self.exponentialSampler(U[:,1], E0) +          \
               weights[2] * self._doubMaxwellianSampler(U[:,2], T1, T2) +  \
               weights[3] * self._relativisticMaxwellianSampler(U[:,3], T)





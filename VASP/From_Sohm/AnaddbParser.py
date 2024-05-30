"""
Parse ANADDB output file
Currently able to parse: EO tensor, Raman tensor, polarities.
"""
from itertools import tee
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

class AnaddbParser:
    def __init__(self, filename, natom=45):
        with open(filename, 'r') as file:
            self.data = file.readlines()
        self.natom = natom


    def get_EO(self):
        EO_tensors = pd.DataFrame(columns=["Mode", "Freq", "EO Tensor"])
        index = self.data.index(" Output of the EO tensor (pm/V) in Voigt notations\n")+6

        for mode in np.arange(4, self.natom*3+1):
            freq = float(re.findall("\s\d+.\d+", self.data[index])[0])

            EO = np.zeros((6,3))
            for i in range(6):
                EO[i] = np.array(self.data[index+i+1].split(), dtype=float)

            EO_tensors.loc[len(EO_tensors.index)] = [mode, freq, EO]
            index += 8
        
        return EO_tensors


    def get_Raman(self, NAC=False):
        Raman_tensors = pd.DataFrame(columns=["Mode", "Freq", "Raman susceptibility"])

        if NAC:
            index = self.data.index(" Raman susceptibility of zone-center phonons, with non-analyticity in the\n")+23
        else:
            index = self.data.index("  Raman susceptibilities of transverse zone-center phonon modes\n")+22

        for mode in np.arange(4, self.natom*3+1):
            freq = float(re.findall("\s\d+.\d+", self.data[index])[0])

            Raman = np.zeros((3,3))
            for i in range(3):
                Raman[i] = np.array(self.data[index+i+1][1:].split(), dtype=float)

            Raman_tensors.loc[len(Raman_tensors.index)] = [mode, freq, Raman]
            index += 6

        return Raman_tensors


    def get_polarity(self):
        polarity = pd.DataFrame(columns=["Mode", "px", "py", "pz", "|p|"])
        index = self.data.index(" Mode effective charges \n")+5
        
        for mode in np.arange(4, self.natom*3+1):
            polarity.loc[len(polarity.index)] = self.data[index][1:].split()
            index += 1

        return polarity


    def plot_r_33(self, filename):
        EO_tensors = self.get_EO()
        r_33 = [i[2,2] for i in EO_tensors["EO Tensor"]]
        # inv_freq_sq = 1 / np.array(EO_tensors["Freq"])**2

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        fig.set_figwidth(15)
        fig.set_figheight(5)

        ax1.plot(EO_tensors["Mode"], r_33, 'b-')
        ax2.plot(EO_tensors["Mode"], EO_tensors["Freq"], 'r-')
        ax1.set_xlabel("phonon mode number")
        ax1.set_ylabel("r_33 (pm/V)", color="b")
        ax2.set_ylabel("mode frequency (cm-1)", color="r")
        ax1.set_title("Mode Decomposed EO Tensor")
        plt.savefig(filename)


    def plot_Raman(self, filename, NAC=False):
        Raman_tensors = self.get_Raman()
        EO_tensors = self.get_EO()
        r_33 = [i[2,2] for i in EO_tensors["EO Tensor"]]
        Raman_33 = [np.abs(i[2,2]) for i in Raman_tensors["Raman susceptibility"]]

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        fig.set_figwidth(15)
        fig.set_figheight(5)

        ax1.plot(Raman_tensors["Mode"], r_33, 'b-')
        ax2.plot(Raman_tensors["Mode"], Raman_33, 'r-')
        # ax1.vlines(x=41, ymin=0, ymax=8)
        ax1.set_xlabel("phonon mode number")
        ax1.set_ylabel("r_33 (pm/V)", color="b")
        ax2.set_ylabel("|Raman_33|", color="r")
        ax1.set_title("Mode Decomposed EO Tensor")
        plt.savefig(filename)


    def plot_polarity(self, filename):
        polarity = self.get_polarity()
        EO_tensors = self.get_EO()
        r_33 = [i[2,2] for i in EO_tensors["EO Tensor"]]
        p = [float(i) for i in polarity["|p|"]]
        modes = [float(i) for i in polarity["Mode"]]

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        fig.set_figwidth(15)
        fig.set_figheight(5)

        ax1.plot(modes, r_33, 'b-')
        ax2.plot(modes, p, 'r-')
        # ax1.vlines(x=41, ymin=0, ymax=8)
        ax1.set_xlabel("phonon mode number")
        ax1.set_ylabel("r_33 (pm/V)", color="b")
        ax2.set_ylabel("|polarity|", color="r")
        ax1.set_title("Mode Decomposed EO Tensor")
        plt.savefig(filename)


        


if __name__ == "__main__":
    SBN = AnaddbParser("SBN_anaddb.abo")
    SBN.plot_Raman("Raman_33.pdf")
    
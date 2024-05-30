"""
Calculates the direct bandgap from a VASP calculation.
Can find the k-point which has the minimum direct bandgap,
or can calulate the bandgap at a specified k-point.
"""
import numpy as np
import pyprocar


def direct_band_gap(outcar="OUTCAR", procar="PROCAR", kpoint=None):
    outcarparser = pyprocar.UtilsProcar()
    fermi = outcarparser.FermiOutcar(outcar)


    procar_file = pyprocar.ProcarParser()
    procar_file.readFile(procar=procar)

    bands = np.array(procar_file.bands)
    sub_bands = np.subtract(bands, fermi)
    valance_top = np.where(sub_bands<0, sub_bands, -np.inf).max(axis=1)
    conduction_bottom = np.where(sub_bands>0, sub_bands, np.inf).min(axis=1)

    bandgaps = conduction_bottom - valance_top
    if kpoint==None:
        return (np.min(bandgaps), np.argmin(bandgaps))
    else:
        return (bandgaps[kpoint], kpoint)



if __name__ == "__main__":
    gap,kpoint = direct_band_gap()
    print("Bandgap:", gap, "eV\t@ K-point:", kpoint)

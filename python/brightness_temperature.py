import numpy as np

from astropy import constants, units

def T_b(wavelength, I_lambda):
        """
        Method to calculate brightness temperature

        Parameters
        ----------
        I_lambda : n-D array (float)
            Specific intensity per wavelength, si units

        wavelength : n-D array (float)
            wavelength in m

        Returns
        -------
        T_b : n-D array (float)
            brightness temperature
        """

        h = constants.h
        c = constants.c
        k_B = constants.k_B

        T_b = h*c/(k_B*wavelength)*np.log(1 + 2*h*c**2/(I_lambda*wavelength**5))**(-1)

        return T_b.to("K")

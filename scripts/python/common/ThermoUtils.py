#!/usr/bin/env python3
#
# Script to generate moist thermodynamic variables from temperature, pressure
# and vapor mixing ratio variables.
#

import numpy as np

####################################################################################################
# Global constants
####################################################################################################
Rd = 287  # Gas constant, dry air (J/Kg/K)
Rv = 461  # Gas constant, vapor (J/Kg/K)

Cp = 1004 # specific heat of dry air (J/Kg/K)

L0 = 2.555e6 # latent heat of vaporization (J/kg)

P0 = 1.0e5 # reference pressure (Pa)
           # 1000 mb => 100000 Pa

T0 = 298 # temperature at P0 in atmosphere (K)
         #  25 C (mean temp at 1000mb in my sims) => 298 K


####################################################################################################
# DryAirDensity
#
# This function calculates the dry air density for a given temperature (K),
# and pressure (Pa).
# 
def DryAirDensity(Temp, Press):
    global Rd

    Rho = Press / ( Rd * Temp)
    
    return Rho

####################################################################################################
# SatVaporPress
#
# This function calculates the saturation vapor pressure for a given temperature (K).
# Based on formula in Tetens (1930).
# 
def SatVaporPress(Temp):

    Esat = 610.78 * np.exp(17.269 * (Temp - 273.15) / (Temp - 35.86))
    
    return Esat

####################################################################################################
# SatVaporMixRatio
#
# This function calculates the saturation vapor mixing ratio for a given temperature (K)
# and pressure (Pa). 
#
def SatVaporMixRatio(Temp, Press):

    Esat = SatVaporPress(Temp)
    Rvapsat = 0.622 * Esat / (Press - Esat)

    return Rvapsat


####################################################################################################
# PotTemp
#
# This function calculates the potential temperature for a given temperature (K)
# and pressure (Pa)
#
def PotTemp(Temp, Press):
    global Rd
    global Cp
    global P0

    Theta = Temp * np.power((P0/Press), (Rd/Cp))

    return Theta

####################################################################################################
# RelHum
#
# This function calculates the relative humidity for a given temperature (K),
# pressure (Pa), and vapor mixing ratio (kg/kg).
#
def RelHum(Temp, Press, Rvap):

    Rvapsat = SatVaporMixRatio(Temp, Press)
    RH = 100 * (Rvap/Rvapsat)

    return RH

####################################################################################################
# EquivPotTemp
#
# This function calculates the equivalent potential temperature for a given temperature (K),
# pressure (Pa), and vapor mixing ratio (kg/kg). Formula is from Bryan (2008).
#
def EquivPotTemp(Temp, Press, Rvap):
    global Rv
    global Cp
    global L0

    Theta = PotTemp(Temp, Press)
    Rvapsat = SatVaporMixRatio(Temp, Press)
    H = Rvap / Rvapsat

    ThetaE = Theta * (H ** (-(Rv * Rvap) / Cp)) * np.exp((L0 * Rvap) / (Cp * Temp))

    return ThetaE

####################################################################################################
# SatEquivPotTemp
#
# This function calculates the saturation equivalent potential temperature for a given
# temperature (K) and pressure (Pa). Formula is adapted from the general theta-e formula
# from Bryan (2008).
#
def SatEquivPotTemp(Temp, Press):
    global Cp
    global L0

    Theta = PotTemp(Temp, Press)
    Rvapsat = SatVaporMixRatio(Temp, Press)

    ThetaEs = Theta * np.exp((L0 * Rvapsat) / (Cp * Temp))

    return ThetaEs

####################################################################################################
# Entropy
#
# This function calculates the moist entropy for a given temperature (K),
# pressure (Pa), and vapor mixing ratio (kg/kg). Formula is from Bryan (2008).
#
def Entropy(Temp, Press, Rvap):
    global Rd
    global Rv
    global Cp
    global L0

    Rvapsat = SatVaporMixRatio(Temp, Press)
    H = Rvap / Rvapsat

    S = (Cp*np.log(Temp/T0)) - (Rd*np.log(Press/P0)) + ((L0 * Rvap)/Temp) - (Rv*Rvap*np.log(H))

    return S

####################################################################################################
# SatEntropy
#
# This function calculates the saturation moist entropy for a given temperature (K),
# pressure (Pa) Formula is adapted from Bryan (2008).
#
def SatEntropy(Temp, Press):
    global Rd
    global Rv
    global Cp
    global L0

    Rvapsat = SatVaporMixRatio(Temp, Press)

    Ssat = (Cp*np.log(Temp/T0)) - (Rd*np.log(Press/P0)) + ((L0 * Rvapsat)/Temp)

    return Ssat


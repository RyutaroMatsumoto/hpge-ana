

"""
Ge_Egap(T::Real)
bandgap energy (eV) of germanium as a function of temperature (K)
"""
function Ge_Egap(T::Real)
    return 0.744 - (4.774*1e-4 * T^2) / (T + 235)
end
# Ge_Egap(77)

"""
Ge_Energy_per_eholePair(T::Real)
energy (eV) required to create an electron-hole pair in germanium as a function of temperature (K)
"""
function Ge_Energy_per_eholePair(T::Real)
    return 2.2 * Ge_Egap(T) + 1.99 * Ge_Egap(T)^(3/2) * exp(4.75 * Ge_Egap(T) / T)
end
# Ge_Energy_per_eholePair(90)

"""
Ge_NumberChargeCarrier(E::Real, T::Real)
number of charge carriers created by an energy E (eV) in germanium at temperature T (K)
"""
function Ge_NumberChargeCarrier(E::Real, T::Real)
return E / Ge_Energy_per_eholePair(T)
end 
# Ge_NumberChargeCarrier(1e6, 77)/1e5


electron_charge = 1.60217662e-19 # C
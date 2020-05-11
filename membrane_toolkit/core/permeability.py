# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Permeability and permeance - related methods
"""

# '''
#                 COMPUTE THE PERMEABILITY COEFFICIENTS
#                 '''
#                 # compute the average slope of the volume change rate
#                 # if the slope is 0 (indicating no data), exclude from averaging
#                 if slope_C != 0:
#                     if slope_D != 0:
#                         avg_slope = (abs(slope_C) + abs(slope_D))/2
#                     else:
#                         avg_slope = abs(slope_C)
#                 else:
#                     avg_slope = abs(slope_D)
#                 # compute water permeance in Liters / m2 / hr / bar
#                 water_permeance = avg_slope / 1000 * 3600 / area / osmotic_pressure.magnitude
#                 # compute the water flow in m3/s
#                 water_flow = avg_slope / 10 ** 6
#                 #compute the salt permeability in m**2/s
#                 salt_permeability = -1 * slope_LHS * thickness * volume / 2 / area

#                 '''
#                 COMPUTE THE PERMEABILITY COEFFICIENTS, ACCOUNTING FOR WATER FLUX
#                 '''
#                 # define the expression for dilute concentration to be solved as a function of permeability
#                 def RHS(permeability):
#                     # calculate the "u" factor
#                     # if u is complex (negative inside the sqrt), then only return the real part
#                     u = scipy.sqrt(permeability **2 * area ** 2 + water_flow **2 * thickness **2)
#                     u = u.real

#                     df['calc'] = -initial_concentration * area * permeability / u *                                (scipy.exp((2*permeability*area-water_flow*thickness-u)*df['time (s)']/(2*volume*thickness))                                 -scipy.exp((2*permeability*area-water_flow*thickness+u)*df['time (s)']/(2*volume*thickness)))
                    
#                     df['diff']=(df['Dilute conc. (M)']-df['calc'])**2
                    
#                 # define a function that iteratively solves permeability to minimize the sum of squared errors between
#                 # measured and calculated dilute concentration
#                 def solve(permeability):

#                     RHS(permeability)
#                     # Scale the output by 1e6 to eliminate machine precision issues with the solver
#                     return df['diff'].sum() * 1e6

#                 # use a nonlinear solver to solve permeability
#                 from scipy.optimize import minimize
#                 result = minimize(solve,salt_permeability,method='Nelder-Mead',tol=1e-8,options={'maxiter':None})
#                 print(result.message)

#                 salt_permeability_osmosis = result.x
#                 # call RHS one more time to make sure the dataframe gets updated
#                 RHS(salt_permeability_osmosis)

# # calculate the water flux
# # 1 bar = 100000 Pa
# results['water flux'] = results['P_w^H (m2.s-1.Pa-1)'] / results['Thickness (um)'] * results['Delta pi, bar'] * 100000

# # calculate the diffusion coefficient of water, Dw
# def dw(PwH,T,phi):
#     return PwH * 8.314 * (T+273.15) / 0.0000182 / phi * (1-phi)
# def apply_dw(row):
#     return dw(row['P_w^H (m2.s-1.Pa-1)'],row['Temperature (degC)'],row['Vol_frac'])
# results['D_w, m2/s'] = results.apply(apply_dw,axis=1)
# #results['D_w, m2/s'] = results['P_w^H (m2.s-1.Pa-1)'] * 8.314 * (results['Temperature (degC)'] + 273.15) / 0.0000182 / results['Vol_frac'] * (1-results['Vol_frac'])

# # calculate the thermo-corrected diffusion coefficient of water with Flory-Huggins
# def dwFH(PwH,T,phi,chi):
#     return PwH * 8.314 * (T+273.15) / 0.0000182 * (1- phi) **2  * (1-2*chi*phi)/phi
# def apply_dwFH(row):
#     return dwFH(row['P_w^H (m2.s-1.Pa-1)'],row['Temperature (degC)'],row['Vol_frac'],row['chi'])
# results['Dw,FH'] = results.apply(apply_dwFH,axis=1)
# #results['Dw,FH'] = results['P_w^H (m2.s-1.Pa-1)'] * 8.314 * (results['Temperature (degC)'] + 273.15) / 0.0000182 * (1-results['Vol_frac'])**2 * (1-2*results['chi']*results['Vol_frac']) / results['Vol_frac']

# # calculate water permeability, Pw
# def Pw(Dw,Kw):
#     return Dw*Kw
# def apply_Pw(row):
#     return Pw(row['D_w, m2/s'],row['K_w'])
# results['P_w (m2/s)']=results.apply(apply_Pw,axis=1)
# #results['P_w (m2/s)'] = results['D_w, m2/s'] * results['K_w']

# # calculate the thermo-corrected water permeability, Pw with Flory-Huggins
# def PwFH(DwFH,Kw):
#     return DwFH*Kw
# def apply_PwFH(row):
#     return PwFH(row['Dw,FH'],row['K_w'])
# results['P_w,FH (m2/s)'] = results.apply(apply_PwFH,axis=1)
# #results['P_w,FH (m2/s)'] = results['Dw,FH'] * results['K_w']

# # calculate the water:salt permeability ratio
# def PwPs(Pw,Ps):
#     return Pw/Ps
# def apply_PwPs(row):
#     return PwPs(row['P_w (m2/s)'],row['Salt Permeability DK w/osmosis, m2/s'])
# results['P_w/P_s'] = results.apply(apply_PwPs,axis=1)
# #results['P_w/P_s'] = results['P_w (m2/s)']/results['Salt Permeability DK w/osmosis, m2/s']

# # calculate the water:salt permeability ratio
# def PwFHPs(PwFH,Ps):
#     return PwFH/Ps
# def apply_PwFHPs(row):
#     return PwFHPs(row['P_w,FH (m2/s)'],row['Salt Permeability DK w/osmosis, m2/s'])
# results['P_w,FH/P_s'] = results.apply(apply_PwFHPs,axis=1)
# #results['P_w,FH/P_s'] = results['P_w,FH (m2/s)']/results['Salt Permeability DK w/osmosis, m2/s']


def water_permeability_hydraulic():
    pass


def water_permeability_diffusive():
    pass


def water_permeance():
    pass


def salt_permeability():
    pass


def salt_permeance():
    pass

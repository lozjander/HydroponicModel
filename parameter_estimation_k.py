# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 17:21:38 2019

@author: ljander
"""

def residuals(p, data_K, t, root_dm, total_dm, PAR, DSR):
    """Calculates the difference between the simulated and measured data (percentage of composition[DM/DM]).

    Parameters
    ----------
    P : list
      vector of parameters that need to be estimated
    data_K :
      the literature data of potassium content
    t : list
      timestep of the simulation
    root_dm :
      dynamic dry weight of root obtained from simulation (vector)
    total_dm :
      dynamic dry weight of the whole tomato plant (vector)
    PAR : list
      photsynthesis active radation
    DSR : list
      Daily Solar Radiation

    Variables
    ---------
    RSA : list
      Root Surface Area
    tc : array
      an array of 't[-1]+1' numbers, equally spaced between '0' - 't[-1]'
    uptake_K :

    cuml_K :

    percentage_K :

    sim_K :

    data_K :

    e_K :


    Constants
    ---------
    k_rsa : float
       coefficient which is used to convert root dry matter to root surface area [0.096]
    EC : float
      electrical conductivity [1.5]

    Returns
    -------
    The difference between the simulated and measured values of uptake of potassium

    Notes
    -----
    # root_dm, total_dm are obtained from simulation of dry matter production model
    # is p a list?
    # spacing around commas in functions
    # are all the parameters just lists???????

        """

    # root_dm, total_dm are obtained from simulation of dry matter production model
    k_rsa = 0.096 # coefficient which convert root dry matter to root surface area
    # a,b,c,d are parameters that need to be estimated
    a = p[0]
    b = p[1]
    c = p[2]
    d = p[3]
    
    # compute the dynamic root surface area (vector)
    RSA = root_dm*k_rsa

    # create time vector for simulation
    tc = np.linspace(0, t[-1], t[-1]+1)

    # parameters for uptake (irrelevent)
    EC = 1.5 # electrical conductivity
    # calculate uptake of potassium
    uptake_K = np.zeros_like(tc)
    for j,k in enumerate(DSR):
     uptake_K[j] = K_uptake(j,RSA[j],PAR[j],EC,a,b,c,d,DSR[j])
    # calculate cumulative uptake of potassium
    cuml_K = np.cumsum(uptake_K)
    # calculate potassium content
    percentage_K = cuml_K/1000/total_dm*100
    # extract data point from simulated
    sim_K = np.zeros_like(data_K)
    for i,l in enumerate(t):
     sim_K[i] = percentage_K[l]
    # delete the last element (artificial data for the last time step to ensure the size of the vector are
    # correct
    sim_K = np.delete(sim_K,29)
    data_K = np.delete(data_K,29)
    # compute difference between data and simulated result (residuals)
    e_K = abs(data_K-sim_K)
    self = e_K
    return self

if __name__ == '__main__':
 # load the ratio between tomato fruit and whole plant (DM) for computing the mineral content of whole
 # plant based on mineral content of leaves and fruits
 r_f = np.load('r_f.npy')
# construct dynamic content data based on literature
 data_K = (np.array([5.55,6.75,6.85,6.75,6.75,6.9,4.91,4.69,4.53,3.68,3.96,3.91,4.02,3.72,3.87,
 5.11*(1-r_f[24*56])+3.76*r_f[24*56],5.51*(1-r_f[24*56])+3.76*r_f[24*56],
 5.53*(1-r_f[24*56])+3.76*r_f[24*56],5.25*(1-r_f[24*56])+3.76*r_f[24*56],
 4.6*(1-r_f[24*56])+3.76*r_f[24*56],4.76*(1-r_f[24*56])+3.76*r_f[24*56],
 4.78*(1-r_f[24*84])+3.76*r_f[24*84],4.78*(1-r_f[24*84])+3.76*r_f[24*84],
 5*(1-r_f[24*84])+3.76*r_f[24*84],4.75*(1-r_f[24*84])+3.76*r_f[24*84],
 5.51*(1-r_f[24*84])+3.76*r_f[24*84],4.2*(1-r_f[24*84])+3.76*r_f[24*84],
 3.72*(1-r_f[24*112])+3.76*r_f[24*112],4.07*(1-r_f[24*112])+3.76*r_f[24*112],
 4.07*(1-r_f[2879])+3.76*r_f[2879]]))/5.5

# time vector of the data point
t = np.array([5*24,10*24,14*24,18*24,20*24,22*24,21*24,21*24,21*24,
 28*24,28*24,28*24,28*24,28*24,28*24,
 56*24,56*24,56*24,56*24,56*24,56*24,
 84*24,84*24,84*24,84*24,84*24,84*24,
 112*24,112*24,2879])
 # initial guess
 p0 = [0.0007,241.491,0.006,2198.3]
 # -- least_squares --
 print('--- least_squares function results --- \n')
 # The known parameters can be passed as args
 lsresult = least_squares(residuals, p0, args=( data_K ,t,root_dm,
 total_dm,PAR,DSR),
 method='trf', bounds=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf]))
 print('Estimated parameters \n', lsresult['x'],'\n')
 # The function returns jac, which seems to be X matrix
 lsX = lsresult['jac']
 lscovx = inv(np.dot(lsX.T,lsX))
 print('Calculated inv(X.T * X) with X=returned jac = \n',lscovx,'\n')
 # The residuals are returned as fun in lsresult
 lsres = lsresult['fun']
 # The calculated variance of residuals is
 lsvarres = 1/(lsX.shape[0]-lsX.shape[1]) * np.dot(lsres.T,lsres)
 # The covariance matrix of parameters is
 lscovp = lsvarres*lscovx
 print('Covariance matrix of parameters \n',lscovp,'\n')
 # The standard deviations of the parameters are:
 lssd = np.sqrt(np.diag(lscovp))
 print('Standard deviations of parameters \n',lssd,'\n')
 print('Objective function value \n',lsresult['cost'])
 print('-'*25,'\n')
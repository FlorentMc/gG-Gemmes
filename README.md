# gG-Gemmes
All codes of the Gemmes project

Instead of using the TwoSect_GoodwinKeenInf and TwoSect_GoodwinKeenInfTay models, we advise you to use directly 
the TwoSect_GoodwinKeenInfTayInv model and if you would like to remove some features, just put the right parameter(s) to zero :
-> sigma = 0; blocks the allocation of investment dynamics (d_theta_i)
-> istar = 0; phi_T = 0; blocks the interest rate r at the rstar value
-> eta_i = 0; blocks the price dynamics

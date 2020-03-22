import numpy as np 
from scipy.constants import speed_of_light, epsilon_0, mu_0

eta_0 = np.sqrt(mu_0 / epsilon_0)

class PyPanel():

    fIni = 1e3
    fEnd = 1000e9
    N = 1001
    omega = np.linspace(fIni,fEnd,N) * 2 * np.pi

    def __init__(self, epsilon_r = 1.0 , mu_r = 1.0, sigma = 0.0, d = 100e-6):
        self.sigma = sigma
        self.d = d
        self.mu = mu_r * mu_0
        self.epsilon_c = -np.complex(0,1)*sigma/PyPanel.omega + \
                         epsilon_0 * epsilon_r
        self.eta = np.sqrt(self.mu / self.epsilon_c)
        self.phi = self._phi()
        self.denom = self._denom()

    def _phi(self):
        gamma = np.complex(0,1) * PyPanel.omega * np.sqrt(self.mu * self.epsilon_c)
        gd = gamma * self.d
        return np.array([[np.cosh(gd),self.eta*np.sinh(gd)], \
                         [np.sinh(gd)/self.eta,np.cosh(gd)]])

    def _denom(self):
        return self.phi[0,0] * eta_0 + self.phi[0,1] + \
               self.phi[1,0] * eta_0**2 + self.phi[1,1] * eta_0

    def T(self):
        return 2 * eta_0 / self.denom

    def R(self):
        num = self.phi[0,0] * eta_0 + self.phi[0,1] -\
              self.phi[1,0] * eta_0**2 - self.phi[1,1] * eta_0 
        return num / self.denom
        

        

class MyClass(object):
    def __init__(self, value):
        self.value = value
        
    def squared(self):
        return self.value ** 2
    
    def __str__(self):
        return "Value = %s" % self.squared()

m = MyClass(10)
print m

"""
class Cosmology(object):
    def __init__(self, OmegaM=0.3, h=0.7):
        self.OmegaK = 0
        self.OmegaM = OmegaM
        self.OmegaL = 1. - OmegaM
        self.H0 = h * 100
        self.DH = C / self.H0

    def _Einv(self, z):
        return 1./np.sqrt(self.OmegaM*((1+z)**3) + self.OmegaK*((1+z)**2) + self.OmegaL)

    def DC(self, z):
        return self.DH*integrate.quad(self._Einv, 0, z)[0] # is this self._Einv the correct syntax?

    def DM(self, z):
        if self.OmegaK == 0:
            DM = self.DC(z)
        else:
            DM = (self.DH/np.sqrt(np.abs(self.OmegaK)))*np.sinh((np.sqrt(np.abs(self.OmegaK))/self.DH)*self.DC(z))
        return DM
            
    def DA(self, z):
        return self.DM(z)/(1+z)

    def DL(self, z):
        return self.DM(z)*(1+z)

    def mu(self, z):
		return 5*np.log10(self.DL(z)*(1e6)/10)

"""



#phivalues = np.zeros([len(xpoints), L], float)
#phivalues [0,:] = -g*eta # append initial conditions at z = 0


"""
phihat = np.fft.rfft(-g*eta)
phihat = phihat*karray*np.sin(omegak*t)/omegak*np.exp(karray*zpoints) # Iuno about this
phiprime = np.fft.irfft(phihat)
phivalues[i,:] = phiprime
"""
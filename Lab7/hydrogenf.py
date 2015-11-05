""" This script solves the hydrogen problem using scipy routine
    integrate.odeint to integrate Schroedinger differential equation
    It coputes the electron density for heaviear atoms in the
    non-interacting limit.
"""
from scipy import *
from scipy import integrate, optimize
from pylab import *


def Schrod_deriv(y, r, l, E, Z):
    """ Routine which is called by odeint to solve the
    Schroedinger equation of the Hydrogen atom.
    We are solving the Eq:
        u'' + (E -l(l+1)/r^2 + 2Z/r ) * u = 0
    which is rewritten as
        y     == [u, u']
        dy/dr == [u', (l(l+1)/r^2 - 2Z/r -E ) u ]
    """
    du2 = y[0]*(l*(l+1)/r**2 - 2*Z/r - E)
    return [y[1], du2]

def SolveSchroedinger(E,R,l,Z):
    """ Integrates Schroedinger equation given a energy E, and angular momentum
    eigenvalue l, and charge Z.
    For stability, we integrate Schroedinger equation starting at infinity
    and integrating down to zero.
    We use initial condition at infinity : u(infinity)=0, and u'(infinity)=small.
    In order that solution is of the order of unity, we normalize it to (roughly) the peak
    value of the function.
    """
    R0 = R[::-1]
    y0=[0, -1e-4] # [value, derivative]
    y = integrate.odeint(Schrod_deriv, y0, R0, args=(l,E,Z) )

    norm = max(y[:len(R)*0.9,0])
    u = y[::-1,0]/norm
    return u
    

def Shoot(E, R, l, Z):
    """ Given a energy E (momentum eigenvalue l, and charge Z), and initial condition
        at infinity (u(infinity)=0, u'(infinity)=small) it finds the value of u(r=0).
        For the state to be bound state, we have the following condition: u(r=0)==0.
    """
    u = SolveSchroedinger(E,R,l,Z)
    return u[0] + (u[1]-u[0])*(0.0-R[0])/(R[1]-R[0])
    

def FindBoundStates(R,lmax,Z,nmax,MaxSteps):
    """ For each momentum eigenvalue l=[0,...lmax], it find first nmax-l bound-states.
    It first brackets roots of function Shoot using a logarithmic mesh dense at zero.
    It later refinds the roots by optimize.brentq method.
    """
    Eb=[]
    for l in range(lmax+1):
        E0 = -1.2*Z**2
        dE = 0.1*Z**2
        uo = Shoot(E0, R, l, Z)
        nfound=0
        for i in range(MaxSteps):
            E0 += dE
            un = Shoot(E0, R, l, Z)
            if uo*un < 0.0:
                Eb0 = optimize.brentq(Shoot, E0-dE, E0, args=(R, l, Z),xtol=1e-17)
                Eb.append([l,Eb0])
                nfound+=1
                print 'Found bound state at ', Eb0
                if (nfound>=nmax-l): break
            uo = un
            dE /= 1.091
    return Eb

def cmpb(x,y):
    """This subroutine is used for sorting bound states.
    It sorts them according to eigenvalue, but if two eigenstates are
    degenerate, it sorts according to angular momentum l.
    """
    if abs(x[1]-y[1])>1e-4:
        return cmp(x[1],y[1])
    else:
        return cmp(x[0],y[0])

R = logspace(-6,2.2,500) # Logarithmic mesh for our integration: [10^-6,...10^2.2]
Z=1                      # Charge of the nucleous
nmax = 4                 # Maximum principal state
lmax = 3                 # Maximum angular momentum
MaxSteps=1000            # Maximum number of steps in searching for the bound state
Ntot = 36                 # nuber of all electrons we put into these hydrogen-like states.

# Computes eigenstates requested by lmax and nmax
bound_states = FindBoundStates(R,lmax,Z,nmax,MaxSteps)
# Sorts eigenestates according to their energy
bound_states.sort(cmpb)

N=0                             # Total number of electrons in the charge density
rho = zeros(len(R),dtype=float) # Total charge density

for state in bound_states:
    l = state[0]    # l of current state
    Ei = state[1]   # Energy of current state

    u = SolveSchroedinger(Ei,R,l,Z)  # u(r) for this bound state
    u2 = u*u                         # normalizing bound state
    norm = integrate.simps(u2,R)     # normalization factor
    u2 *= 1/norm

    dN = 2*(2*l+1)                   # degeneracy of this particular bound state
    if N+dN<Ntot:                    
        ferm = 1                     # This bound state will be fully filled
    else:
        ferm = (Ntot-N)/(dN+0.0)     # This bound state needs only partial filling up to N=Ntot
        
    rho += u2*dN*ferm                # Updating charge density
    N += dN                          # ... and number of electrons
    print 'adding', state, 'ferm=', ferm, N, Ntot
    if N>=Ntot: break
    
    
plot(R,rho,'o-')
a = axis()
axis([0.0,30.,a[2],a[3]])
show()
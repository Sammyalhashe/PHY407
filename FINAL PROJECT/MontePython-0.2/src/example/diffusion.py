#!/usr/bin/env python

import pypar
import sys, math, Numeric, MLab, os, pickle, numarray
from MontePython import MontePython

try:
    ifile = sys.argv[1]
except:
    ifile = '' # let MontePython instance make its own file

d = MontePython(pypar,ifile)

# 1 particle and alpha=0.5 yields a harmonic oscillator with
# energy == 3/2:
d.params[0] = 0.5 
d.reset_params()
#d.silent = True # to avoid too much noise
if not os.access('checkpoint.dat',os.F_OK):
    d.blocks = 5
    d.tau_list = [[0.001,100,d.blocks],]
    if os.access('results_'+d.outfile,os.F_OK):
        os.unlink('results_'+d.outfile)
else:
    d.blocks = d.tau_list[0][2]

blocks = d.blocks

loc_tmp_nrg = Numeric.empty(3,'f')
tmp_nrg = Numeric.empty(3*d.numproc,'f')
tmp_nrg.shape = (d.numproc,3)


def mk_chkpt_file(steps,termalization):
    file_str = 'checkpoint.dat'
    tmp = []
    for i in d.params:
        tmp += [i]
    d.params = tmp
    f = open(file_str,'w')
    init = """set wave type       = %s  #type of wave function
set numerical       = %d      #numerical
set metropolis      = %d      #metropolis=1,linear=0,quadric=2
set save pos        = %d      #save_pos=1
set no. of walkers  = %d      #number of walkers
set desired walkers = %d      #desired number of walkers
set particles       = %d      #particles
set dimensions      = %d      #dimensions
set steps           = %d      #time steps
set termalization   = %d      #thermalization
set metropolis termalization = %d      #metropolis_termalization
set step length     = %15.10e #step length
set idum            = %d      #idum
set tau             = %15.10e #tau
set D               = %d      #D
set e_trial         = %15.10e #trial energy
# params are (in order): initial param no. 1, dummy, particle size,
# dummy, healing length, lambda, kappa, 
# do Jackson-Feenberg
set params          = %s
        
# parameters only used in python version:
set scale	    = %d #scale used in calculation of energy
set outfile         = '%s'

#extensions from checkpointing from diffusion2.py
set tau list        = %s
"""%(d.wave_type, d.numerical, d.metropolis, d.save_pos,
     d.no_of_walkers, d.desired_walkers, d.particles,
     d.dimensions, steps, termalization,
     0, d.step_length, d.ran.get_idum(),
     d.tau, d.D, d.e_trial, d.params, d.scale,
     d.outfile, d.tau_list)
    f.write(init)
    f.close()

def checkpoint(steps,termalization):
    mk_chkpt_file(steps,termalization)
    d.save_walkers()

def timestep(i_step):
    M = d.no_of_walkers;

    d.monte_carlo_steps()
    
    d.update = False

    #bring out your dead
    for i in xrange(len(d.w_block)-1,-1,-1):
        if d.w_block.isDead(i):
            del(d.w_block[i])  #removing walker
        else:
            while d.w_block.tooAlive(i):
                d.w_block += d.w_block[i]
                d.w_block.calmWalker(len(d.w_block)-1)
                d.w_block.madeWalker(i)

    d.load_balancing()

    nrg = 0.; pot_nrg = 0.; vort_nrg = 0
    loc_nrg = 0.; loc_pot_nrg = 0.; loc_vort_nrg = 0

    d.time_step_counter += 1

    loc_nrg     = d.w_block.warray.getLocalEnergies(*d.num_args)
    loc_pot_nrg = d.w_block.warray.getOtherEnergies(d.pot_e)
    if hasattr(d,'vort_e'):
        loc_vort_nrg = d.w_block.warray.getOtherEnergies(d.vort_e)

    loc_nrg      /= float(d.no_of_walkers*d.particles*d.scale)
    loc_pot_nrg  /= float(d.no_of_walkers*d.particles*d.scale)
    loc_vort_nrg /= float(d.no_of_walkers*d.particles*d.scale)
    loc_tmp_nrg[0] = loc_nrg
    loc_tmp_nrg[1] = loc_pot_nrg
    loc_tmp_nrg[2] = loc_vort_nrg
    d.pypar.gather(loc_tmp_nrg,
                   d.master_rank,
                   buffer=tmp_nrg)

    if d.master:
        nrg=float(Numeric.sum(tmp_nrg[:,0]))
        pot_nrg=float(Numeric.sum(tmp_nrg[:,1]))
        vort_nrg=float(Numeric.sum(tmp_nrg[:,2]))
        d.energy  += nrg
        d.energy2 += nrg*nrg

        if d.time_step_counter > d.termalization:
            d.observables[i_step,0] = nrg
            d.observables[i_step,1] = pot_nrg
            d.observables[i_step,2] = vort_nrg


        if not hasattr(d,'silent'):
            print "Step no.", i_step+1, "\tnrg", nrg, "\tpot_nrg",\
                  pot_nrg, "\tvort_nrg", vort_nrg, "\tno. of walkers",\
                  d.no_of_walkers

    # adjust trial energy (and no. of walkers)
    nrg = -.5*math.log(float(d.no_of_walkers)/float(M))/float(d.tau)
    d.e_trial += nrg

    d.num_args[-1] = d.update


#set initial walker positions:
d.uni_dist()

for i in xrange(d.metropolis_termalization):
    d.metropolis_steps()

nrg = 0.

nrg = d.w_block.warray.getLocalEnergies(*d.num_args)

d.e_trial = d.all_reduce(nrg/d.no_of_walkers)

if d.master:
    print d.e_trial/d.particles/d.scale
    print d.accept/float(d.trials)


first_tau=True
while d.tau_list: #loop over taus
    tau = d.tau_list[0]
    d.tau   = tau[0]
    d.steps = tau[1]
    d.blocks = tau[2]
    if not os.access('observables_tau%g.dat'%d.tau,os.F_OK):
        if d.master:
            observables = Numeric.zeros((3,d.blocks),'float')
    elif d.master:
        f = open('observables_tau%g.dat'%d.tau,'r')
        observables = pickle.load(f)
        f.close()

    # resetting walkers if not first tau:
    if first_tau:
        first_tau=False
    else:
        d.uni_dist()
        for i in xrange(d.metropolis_termalization):
            d.metropolis_steps()
        
    while tau[2] > 0: #loop over blocks
        d.reset()
        d.time_step_counter = 0
        for i in range(d.termalization):
            checkpoint(d.steps,d.termalization)
            timestep(i)

        d.energy = 0; d.energy2 = 0

        for i in range(d.steps):
            d.tau_list[0][1] = d.steps-i
            checkpoint(d.steps-i,0)
            timestep(i)

        if d.master:
            d.energy  /= float(d.steps)
            d.energy2 /= float(d.steps)
            d.energy2 -= d.energy*d.energy

            observables[0,d.blocks-tau[2]] = d.energy

            f = open('observables_tau%g.dat'%d.tau,'w')
            pickle.dump(observables,f)
            f.close()
            
            print 'tau',d.tau,'block',d.blocks-tau[2],\
                  'no. of walkers',d.no_of_walkers
            print "energy = %g +/- %g"%(d.energy,\
                                        math.sqrt(d.energy2/\
                                                  float(d.steps)))
            print "sigma = %g"%d.energy2

        tau[2] -= 1

    
    if d.master:
        print 'tau',d.tau,'no. of walkers',d.no_of_walkers
        print "energy = %g +/- %g"%(MLab.mean(observables[0]),\
                                    math.sqrt(MLab.std(observables[0])/\
                                              float(len(observables[0]))))
        print "sigma = %g"%MLab.std(observables[0])
        os.remove('observables_tau%g.dat'%d.tau)
        f = open('results_'+d.outfile,'a')
        f.write('tau %d no. of walkers %d'%(d.tau,d.no_of_walkers))
        f.write("energy = %g +/- %g"%(MLab.mean(observables[0]),\
                                    math.sqrt(MLab.std(observables[0])/\
                                              float(len(observables[0])))))
        f.write("sigma = %g\n"%MLab.std(observables[0]))
        f.close()

    del d.tau_list[0]
    
if d.master: os.remove('checkpoint.dat')
os.remove('walkers%d_'%d.myrank+d.outfile)
pypar.Finalize()

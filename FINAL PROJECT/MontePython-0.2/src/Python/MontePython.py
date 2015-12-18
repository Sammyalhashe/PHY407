#!/usr/bin/env python

from MontePython_cxx import *
import re,math,Numeric,pickle,os,WalkerArray,copy,time


class MontePython:
    """
Parallellized Diffusion Monte Carlo Simulation
Example on usage:
>>> import pypar
Pypar (version 1.9.1) initialised MPI OK with 1 processors
>>> import sys, math, Numeric, MLab, os, pickle
>>> from MontePython import MontePython
>>> d = MontePython(pypar)
>>> d.reset_params()
>>> d.silent = True # to avoid too much noise
>>> d.blocks = 5
>>> d.tau_list = [[0.001,100,d.blocks],]
>>> blocks = d.blocks
>>> loc_tmp_nrg = Numeric.empty(3,'f')
>>> tmp_nrg = Numeric.empty(3*d.numproc,'f')
>>> tmp_nrg.shape = (d.numproc,3)
>>> def timestep(i_step):
...     M = d.no_of_walkers;
...     d.monte_carlo_steps()
...     d.update = False
...     #bring out your dead
...     for i in xrange(len(d.w_block)-1,-1,-1):
...         if d.w_block.isDead(i):
...             del(d.w_block[i])  #removing walker
...         else:
...             while d.w_block.tooAlive(i):
...                 d.w_block += d.w_block[i]
...                 d.w_block.calmWalker(len(d.w_block)-1)
...                 d.w_block.madeWalker(i)
...     d.load_balancing()
...     nrg = 0.; pot_nrg = 0.; vort_nrg = 0
...     loc_nrg = 0.; loc_pot_nrg = 0.; loc_vort_nrg = 0
...     d.time_step_counter += 1
...     loc_nrg     = d.w_block.warray.getLocalEnergies(*d.num_args)
...     loc_pot_nrg = d.w_block.warray.getOtherEnergies(d.pot_e)
...     if hasattr(d,'vort_e'):
...         loc_vort_nrg = d.w_block.warray.getOtherEnergies(d.vort_e)
...     loc_nrg      /= float(d.no_of_walkers*d.particles*d.scale)
...     loc_pot_nrg  /= float(d.no_of_walkers*d.particles*d.scale)
...     loc_vort_nrg /= float(d.no_of_walkers*d.particles*d.scale)
...     loc_tmp_nrg[0] = loc_nrg
...     loc_tmp_nrg[1] = loc_pot_nrg
...     loc_tmp_nrg[2] = loc_vort_nrg
...     d.pypar.gather(loc_tmp_nrg,
...                    d.master_rank,
...                    buffer=tmp_nrg)
...     if d.master:
...         nrg=float(Numeric.sum(tmp_nrg[:,0]))
...         pot_nrg=float(Numeric.sum(tmp_nrg[:,1]))
...         vort_nrg=float(Numeric.sum(tmp_nrg[:,2]))
...         d.energy  += nrg
...         d.energy2 += nrg*nrg
...         if d.time_step_counter > d.termalization:
...             d.observables[i_step,0] = nrg
...             d.observables[i_step,1] = pot_nrg
...             d.observables[i_step,2] = vort_nrg
...         if not hasattr(d,'silent'):
...             print "Step no.", i_step+1, "\tnrg", nrg, "\tpot_nrg",\
                   pot_nrg, "\tvort_nrg", vort_nrg, "\tno. of walkers",\
                   d.no_of_walkers
...     # adjust trial energy (and no. of walkers)
...     nrg = -.5*math.log(float(d.no_of_walkers)/float(M))/float(d.tau)
...     d.e_trial += nrg
...     d.num_args[-1] = d.update
...
>>> d.uni_dist()
>>> for i in xrange(d.metropolis_termalization):
...     d.metropolis_steps()
...
>>> nrg = 0.
>>> nrg = d.w_block.warray.getLocalEnergies(*d.num_args)
>>> d.e_trial = d.all_reduce(nrg/d.no_of_walkers)
>>> if d.master:
...     print d.e_trial/d.particles/d.scale
...     print d.accept/float(d.trials)
...
2.51834776963
0.526733333333
>>> first_tau=True
>>> while d.tau_list: #loop over taus
...     tau = d.tau_list[0]
...     d.tau   = tau[0]
...     d.steps = tau[1]
...     d.blocks = tau[2]
...     if not os.access('observables_tau%g.dat'%d.tau,os.F_OK):
...         if d.master:
...             observables = Numeric.zeros((3,d.blocks),'float')
...     elif d.master:
...         f = open('observables_tau%g.dat'%d.tau,'r')
...         observables = pickle.load(f)
...         f.close()
...     # resetting walkers if not first tau:
...     if first_tau:
...         first_tau=False
...     else:
...         d.uni_dist()
...         for i in xrange(d.metropolis_termalization):
...             d.metropolis_steps()
...     while tau[2] > 0: #loop over blocks
...         d.reset()
...         d.time_step_counter = 0
...         for i in range(d.termalization):
...             timestep(i)
...         d.energy = 0; d.energy2 = 0
...         for i in range(d.steps):
...             d.tau_list[0][1] = d.steps-i
...             timestep(i)
...         if d.master:
...             d.energy  /= float(d.steps)
...             d.energy2 /= float(d.steps)
...             d.energy2 -= d.energy*d.energy
...             observables[0,d.blocks-tau[2]] = d.energy
...             f = open('observables_tau%g.dat'%d.tau,'w')
...             pickle.dump(observables,f)
...             f.close()
...             print 'tau',d.tau,'block',d.blocks-tau[2],\
                   'no. of walkers',d.no_of_walkers
...             print "energy = %g +/- %g"%(d.energy,\
                                         math.sqrt(d.energy2/\
                                                   float(d.steps)))
...             print "sigma = %g"%d.energy2
...         tau[2] -= 1
...     if d.master:
...         print 'tau',d.tau,'no. of walkers',d.no_of_walkers
...         print "energy = %g +/- %g"%(MLab.mean(observables[0]),\
                                     math.sqrt(MLab.std(observables[0])/\
                                               float(len(observables[0]))))
...         print "sigma = %g"%MLab.std(observables[0])
...     del d.tau_list[0]
...
tau 0.001 block 0 no. of walkers 100
energy = 2.50824 +/- 0.000804816
sigma = 6.47729e-05
tau 0.001 block 1 no. of walkers 100
energy = 2.51304 +/- 0.000524061
sigma = 2.7464e-05
tau 0.001 block 2 no. of walkers 100
energy = 2.53004 +/- 0.000766155
sigma = 5.86993e-05
tau 0.001 block 3 no. of walkers 100
energy = 2.5335 +/- 0.00120665
sigma = 0.0001456
tau 0.001 block 4 no. of walkers 100
energy = 2.52183 +/- 0.000481917
sigma = 2.32244e-05
tau 0.001 no. of walkers 100
energy = 2.52133 +/- 0.0464206
sigma = 0.0107744
>>> pypar.Finalize()
>>>
    """

#####################init functions:#########################

    def make_infile(self):
        file_str = 'initialize.dat'
        f = open(file_str,'w')
        init = """set wave type       = DuBois  #type of wave function
set numerical       = 0       #numerical
set metropolis      = 1       #metropolis=1,linear=0,quadric=2
set save pos        = 1       #save_pos=1
set no. of walkers  = 100     #number of walkers
set desired walkers = 100     #desired number of walkers
set particles       = 2       #particles
set dimensions      = 3       #dimensions
set steps           = 100     #time steps
set termalization   = 50      #thermalization
set metropolis termalization = 150     #metropolis_termalization
set step length     = 1.5     #step length
set idum            = 68758   #idum
set tau             = 5e-05   #tau
set D               = 1       #D
set e_trial         = 6.01042 #trial energy
# params are (in order): initial param no. 1, dummy, particle size,
dummy, healing length, lambda, kappa, 
do Jackson-Feenberg
set params          = [0.53, 0, 0.15155, 0, 200, 2.82842712475, 1, 0]
        
# parameters only used in python version:
set scale	    = 2  #scale used in calculation of energy
set outfile         = 'output.dat'
"""
        f.write(init)
        f.close()
        return file_str

    def __init__(self, pypar, file_str = ''):
        """Initialization of MontePython"""
        self.pypar = pypar
        self.file_str = file_str
        if not self.file_str: 
            self.file_str = self.make_infile() # if no file, make one
        try: f = open('checkpoint.dat')
        except IOError:
            self.file_str = self.file_str
        else:
            if self.pypar.rank() == 0:
                print 'found checkpoint file, will start from checkpoint..'
            self.file_str = 'checkpoint.dat'
        self.__init_file(self.file_str)
        self.nabla  = Nabla() # gradient
        self.nabla2 = Nabla2() # laplacian
        self.move = DuBoisMove()
        self.__set_params()
        self.__set_move_params()
        self.__init_mpi()
        self.__set_warray()
        self.__set_work_data()
        self.t0 = time.time()
        self.first_balance = True

    def __set_params(self):
        """Function setting parameters for functors"""
        # convert params to QickArray
        conv = Convert_QickArray()
        qparams = conv.py2my(self.params)
        # and set the params
        self.local_e.setParams(qparams)
        self.pot_e.setParams(qparams)
        self.wf.setParams(qparams)
        self.wf_all.setParams(qparams)
        self.q_force.setParams(qparams)
        if hasattr(self,'vort_e'):
            self.vort_e.setParams(qparams)


    def __set_move_params(self):
        """Function setting parameters for movement of walkers,
        This is separated from set_params as tau is likely to
        change during simulation"""
        self.move_params = [self.params[2],\
                            math.sqrt(2.*self.D*self.tau)]
        conv = Convert_QickArray()
        qmove_params = conv.py2my(self.move_params)
        self.move.setParams(qmove_params)


    def __set_w_p_node(self):
        self.loc_walkers = [self.no_of_walkers/self.numproc]*self.numproc
        for i in xrange(self.numproc):
            for j in range(self.no_of_walkers%self.numproc):
                if(i==j):
                    self.loc_walkers[self.numproc-1-i] += 1


    def __init_mpi(self):
        """mpi specific initialization"""
        self.numproc = self.pypar.size()
        self.myrank  = self.pypar.rank()
        self.master_rank = 0
        self.master = self.myrank == self.master_rank
        self.slave  = not self.master
        #each node must have different idum
        self.ran = Random(self.idum)#*(self.myrank+1)) 
        # no. of local walkers:
        self.__set_w_p_node()


    def __set_w_block(self,buffer=None):
        #function for setting local walker array
        #if buffer is a Numeric array, this should be a walker array
        #it will be used instead of creating new walkers
        args = [self.loc_walkers[self.myrank],
                self.particles,self.dimensions]
        if buffer:
            if isinstance(buffer,Numeric.ArrayType):
                args += [buffer]
                sizeof_walker = args[1]*args[2]*9+9
                self.loc_walkers[self.myrank] = len(buffer)/sizeof_walker
                args[0] = self.loc_walkers[self.myrank]
            else:
                raise 'wrong buffer!'
        self.w_block = WalkerArray.WalkerArray(*args)

    def __set_warray(self):
        self.__set_w_block()
    

    def __set_work_data(self):
        """function setting work arrays etc."""
        self.accept        = 0
        self.trials        = 0
        self.accept_walker = 0
        self.energy        = 0.
        self.energy2       = 0. #energy squared
        self.update        = True
        if self.numerical:
            self.num_args = [self.wf_all,self.nabla,self.nabla2,self.pot_e,\
                             True]
        else: self.num_args = [self.local_e,True]
        #total, potential and vortex energy:
        self.observables = Numeric.zeros((self.steps,3),'d') 
        
    def __wave_type2funcs(self, text):
        code  = 'self.local_e = '+text+'LocalEnergy();'
	code += 'self.pot_e = '+text.split('Vortex')[0]+'PotentialEnergy();'
        code += 'self.wf = '+text+'Wave();'
        code += 'self.wf_all = '+text+'WaveAll();'
        code += 'self.q_force = '+text+'QForce();'
        if text.find('Vortex') > 0:
            code += 'self.vort_e = '+text+'PotentialEnergy();'
        return code


    def __init_file(self, file_str):
        ifile = open(file_str)

        line_re = re.compile(r'set (.*?)=(.*)') 
        value_re = re.compile(r'(.*)#')

        for line in ifile:
            m = line_re.search(line)
            if m:
                variable = m.group(1).strip()
                value = m.group(2).strip()
            v = value_re.search(value)
            if v:
                value = v.group(1).strip()
            # replacing illegal chars in variable names with _
            variable = variable.replace(' ','_')
            # (the variables should stay in this class)
            variable = variable.replace('.','')
            if variable == 'wave_type':
                code = self.__wave_type2funcs(value)
                code += "self." + variable + " = \'" + value + "\'"
            else:
                code = "self." + variable + " = " + value
            exec code

    def __mk_chkpt_file(self,steps,termalization):
        file_str = 'checkpoint.dat'
        tmp = []
        for i in self.params:
            tmp += [i]
        self.params = tmp
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
"""%(self.wave_type, self.numerical, self.metropolis, self.save_pos,
     self.no_of_walkers, self.desired_walkers, self.particles,
     self.dimensions, steps, termalization,
     0, self.step_length, self.idum,
     self.tau, self.D, self.e_trial, self.params, self.scale,
     self.outfile)
        f.write(init)
        f.close()
        return file_str

######################user functions:########################

    def uni_dist(self,reset=False):
        """Function distributing walkers uniformly
        (or loading checkpointed walkers)"""
        if (self.file_str == 'checkpoint.dat') and \
           os.access(('walkers%d_'%self.myrank)+
                     self.outfile,os.F_OK) and (not reset):
            self.load_walkers()
        else:
            self.w_block.warray.UniDist(self.ran)
            

    def reset_params(self):
        """Function allowing users to change parameters sat in
        the initialization file. Convenient if self has made a
        default file. Note that it is safer to edit initialize.dat
        and start over again."""
        self.__set_params()

    def reset_work_data(self):
        """Function allowing blocking without deleting everything."""
        self.__set_work_data()

    def reset(self):
        self.reset_params()
        self.reset_work_data()


    def save_walkers(self):
        if not hasattr(self,'outfile'):
            raise 'no outfile defined; no walkers saved'
        f = open(('walkers%d_'%self.myrank)+self.outfile+'tmp','w')
        pickle.dump(self.w_block[:],f,
                    protocol=pickle.HIGHEST_PROTOCOL)
        f.close()
        os.system('mv -f walkers%d_'%self.myrank+\
                  self.outfile+'tmp walkers%d_'%self.myrank+self.outfile)

	if hasattr(self,'densdir'):
            #if not os.access(self.densdir,os.F_OK):
            #    os.mkdir(self.densdir)
            os.system("mkdir -p %s"%self.densdir)
            f = open(os.path.join(self.densdir,
                                  'dens%d_'%self.myrank+self.outfile+'tmp'),'w')
            pickle.dump(self.loc_dens_array,f,protocol=pickle.HIGHEST_PROTOCOL)
            f.close()
            os.system('mv -f '+
                      os.path.join(self.densdir,
                                   'dens%d_'%self.myrank+self.outfile+'tmp')+' '+
                      os.path.join(self.densdir,
                                   'dens%d_'%self.myrank+self.outfile))
            f = open(os.path.join(self.densdir,
                                  'dens_pure%d_'%self.myrank+self.outfile+'tmp'),'w')
            pickle.dump(self.loc_dens_pure,f,protocol=pickle.HIGHEST_PROTOCOL)
            f.close()
            os.system('mv -f '+
                      os.path.join(self.densdir,
                                   'dens_pure%d_'%self.myrank+self.outfile+'tmp')+' '+
                      os.path.join(self.densdir,
                                   'dens_pure%d_'%self.myrank+self.outfile))

    def load_walkers(self):
        if not hasattr(self,'outfile'):
            raise 'no outfile defined; no walkers loaded'
        if not os.access('walkers%d_'%self.myrank+self.outfile,os.F_OK):
            print 'Can\'t find file '+'walkers%d_'%self.myrank\
                  +self.outfile+', ignored'
        else:
            f = open('walkers%d_'%self.myrank+self.outfile,'r')
            walkers = pickle.load(f)
            self.__set_w_block(walkers)
            f.close()
            if hasattr(self,'densdir'):
                f = open(os.path.join(self.densdir,
                                      'dens%d_'%self.myrank+self.outfile),'r')
                self.loc_dens_array = pickle.load(f)
                f.close()
                f = open(os.path.join(self.densdir,
                                      'dens_pure%d_'%self.myrank+self.outfile),'r')
                self.loc_dens_pure = pickle.load(f)
                f.close()
            

        self.no_of_walkers = int(self.all_reduce(len(self.w_block)))

        self.__set_w_p_node()
        self.update = False


    def checkpoint(self,steps,termalization):
        self.__mk_chkpt_file(steps,termalization)
        self.save_walkers()


    def all_reduce(self,x):
        xs = self.pypar.gather(Numeric.array([x]),self.master_rank)
        return self.pypar.broadcast(float(Numeric.sum(xs)),self.master_rank)


    def metropolis_steps(self):
        """Wrap around MCpp().MetropolisStep()"""
        conv = Convert_QickArray()
        counters = Numeric.array([self.trials, self.accept])
        counters = conv.py2my(counters)
        args = [self.ran, self.move, self.wf,\
                conv.py2my(self.params), counters,\
                self.step_length]
        self.w_block.warray.MetropolisSteps(*args)
        counters = conv.my2py(counters)
        self.trials = int(counters[0])
        self.accept = int(counters[1])


    def monte_carlo_steps(self):
        """Wrap around MCpp().MonteCarloStep()"""
        conv = Convert_QickArray()
        counters = Numeric.array([self.trials, self.accept,
                                  self.accept_walker])
        counters = conv.py2my(counters)
        qparams  = conv.py2my(self.params)
        args = [self.ran, self.wf, self.wf_all, self.local_e,
                self.q_force, self.move, self.pot_e, self.nabla,
                self.nabla2, qparams, counters,
                self.metropolis, self.D, self.tau, 
                self.e_trial, self.update, self.numerical]
        self.w_block.warray.MonteCarloSteps(*args)
        counters = conv.my2py(counters)
        self.params = conv.my2py(qparams)
        self.trials        = int(counters[0])
        self.accept        = int(counters[1])
        self.accept_walker = int(counters[2])
        

    def NumericFloat2IntList(self,x):
        return Numeric.array(x.tolist(),'i').tolist()

    def __find_opt_w_p_node(self):
        """Help function for load_balancing()
        Finds and returns the optimal number of walkers
        per node for a possibly non-uniform set of nodes
        """
        timings = self.pypar.gather(Numeric.array([abs(self.t1-self.t0)]),
                                    self.master_rank)
        tmp_timings = copy.deepcopy(timings)
        timings = self.pypar.broadcast(tmp_timings,
                                       self.master_rank)

        C = self.no_of_walkers/sum(1./timings)

        tmp_loc_walkers = C/timings

        self.loc_walkers = self.NumericFloat2IntList(tmp_loc_walkers)
        remainders = tmp_loc_walkers-self.loc_walkers

        while sum(self.loc_walkers) < self.no_of_walkers:
            maxarg = Numeric.argmax(remainders)
            self.loc_walkers[maxarg] += 1
            remainders[maxarg] = 0

        if self.master and self.first_balance and\
               not hasattr(self,'silent'):
            print timings
            print self.loc_walkers

        return self.loc_walkers
    

    def load_balancing(self):
        """Function balancing load between nodes"""
        self.t1 = time.time()
        w_numbers = self.pypar.gather(Numeric.array([len(self.w_block)]),
                                      self.master_rank)
        tmp_w_numbers = copy.deepcopy(w_numbers)
        w_numbers = self.pypar.broadcast(tmp_w_numbers,self.master_rank)

        self.no_of_walkers = Numeric.sum(w_numbers)

        self.__find_opt_w_p_node()
        #self.__set_w_p_node()
        self.first_balance = False
        balanced = Numeric.array(self.loc_walkers)

        difference = w_numbers-balanced

        diff_sort = Numeric.argsort(difference)
        prev_i_min = diff_sort[0]
        
        while sum(abs(difference))!=0:
            diff_sort = Numeric.argsort(difference)
            i_max = diff_sort[-1]
            i_min = diff_sort[0]
            
            if i_min == prev_i_min and i_max!=diff_sort[1]:
                i_min = diff_sort[1]
            
            if self.myrank==i_max:
                self.pypar.send(self.w_block[balanced[i_max]:],i_min)
                args = [balanced[i_max],
                        self.particles,
                        self.dimensions,
                        self.w_block[0:balanced[i_max]]]
                self.w_block = WalkerArray.WalkerArray(*args)
            elif self.myrank==i_min:
                recv_buff = self.pypar.receive(i_max)
                args = [len(self.w_block)+difference[i_max],
                        self.particles,
                        self.dimensions,
                        Numeric.concatenate((self.w_block[:],recv_buff))]
                self.w_block = WalkerArray.WalkerArray(*args)
            difference[i_min]+=difference[i_max]
            difference[i_max]=0
            prev_i_min = i_min
        #self.t0 = self.t1


def _test():
    import doctest, MontePython
    return doctest.testmod(MontePython)

if __name__ == '__main__':
    _test()

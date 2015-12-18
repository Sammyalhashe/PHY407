#!/usr/bin/env python

import math,Numeric,pickle,MontePython_cxx,copy

class WalkerArray(object):
    """
    Class for walker array book keeping,
    walker info is stored in a Numeric array,
    
    Example of usage:
    >>> from WalkerArray import WalkerArray
    >>> wa = WalkerArray(20,1,3)
    >>> from MontePython_cxx import DuBoisLocalEnergy
    >>> local_e = DuBoisLocalEnergy()
    >>> params = [0.53, 0, 0.15155, 0, 200, 1, 1, 0]
    >>> from MontePython_cxx import Convert_QickArray
    >>> conv = Convert_QickArray()
    >>> local_e.setParams(conv.py2my(params))
    >>> from MontePython_cxx import Random
    >>> ran = Random(42)
    >>>
    """

    def __warray(self):
        """function returning an array with info on walkers"""
        #arrays that Walker will handle itself
        skip_size = self.particles*self.dimension*9 

        w = Numeric.zeros(self.sizeof_walker*self.size,'f')
        last_i=0
        stride=self.sizeof_walker
        # set 
        for i in Numeric.arange(0,len(w),stride):
            index = i+skip_size
            w[index] = 0 # e_l
            index += 1
            w[index] = 0 # wave_func
            index += 1
            w[index] = self.particles
            index += 1
            w[index] = self.dimension
            index += 1
            w[index] = 0 # particles_moved
            index += 1
            w[index] = 0 # replicate
            index += 1
            w[index] = self.sizeof_walker
            index += 1
            w[index] = 0 # dead
            index += 1
            w[index] = 0 # old_eq_new
            index += 1

        return w

    def __init__(self, size, particles, dimension, walkerchunk=None):
        self.particles     = particles; self.dimension = dimension
        self.positionsize  = self.particles*self.dimension
        self.sizeof_walker = self.positionsize*9+9
        self.size          = size
        self.arraysize     = size
        if isinstance(walkerchunk,Numeric.ArrayType):
            if len(walkerchunk)==self.size*self.sizeof_walker:
                self.walkers = walkerchunk
            else:
                raise "wrong size of walkerchunk; should be %d but is %d"\
                      %(self.size*self.sizeof_walker,len(walkerchunk))
        else:
            self.walkers       = self.__warray()
        self.warray = MontePython_cxx.Warray()
        self.walkers = self.warray.pyInitialize(size, particles, dimension,
                                                self.walkers)
        # set up var index:
        self.skip_arrays   = self.sizeof_walker-9
        self.vars = {'new_pos'            : self.positionsize,
                     'o_pure'             : self.positionsize*5,
                     'o_pure_new'         : self.positionsize*6,
                     'o_pure_further'     : self.positionsize*7,
                     'o_pure_further_new' : self.positionsize*8,
                     'e_l'                : self.skip_arrays,
                     'wave_func'          : self.skip_arrays+1,
                     'particles'          : self.skip_arrays+2,
                     'dimension'          : self.skip_arrays+3,
                     'particles_moved'    : self.skip_arrays+4,
                     'replicate'          : self.skip_arrays+5,
                     'size'               : self.skip_arrays+6,
                     'dead'               : self.skip_arrays+7,
                     'old_eq_new'         : self.skip_arrays+8,
                     }

    def __len__(self):
        return self.size

    def _slicer(self,slic):
        start=0;step=1;stop=0
        if slic.stop:
            if slic.stop < len(self)+1: stop = slic.stop
            else: stop = len(self)
        else: stop = len(self)
        if slic.start: start = slic.start
        else: start = 0
        if slic.step: step = slic.step
        else: step = 1
        if step>1:
            raise "Step size %d > 1!"%step\
                  +"Only continous blocks of walkers are implemented"
        size = len(range(start,stop,step))
        return start,stop,step,size
        

    def __getitem__(self,w_nr):
        """getting walker chunk, if slice, size of walker is defined
        by slice as len(range(start,stop,step))
        """
        if isinstance(w_nr,slice):
            start,stop,step,size = self._slicer(w_nr)
            start *= self.sizeof_walker
            stop  *= self.sizeof_walker
            #step  *= self.sizeof_walker
            return self.walkers[start:stop:step]
        elif abs(w_nr) < len(self):
            start = w_nr*self.sizeof_walker
            stop  = (w_nr+1)*self.sizeof_walker
            return self.walkers[start:stop]

    def __setitem__(self,w_nr,walkerchunk):
        """setting walker. if w_nr is a slice, walker must be
        of length defined by slice as len(range(start,stop,step))
        """
        if not isinstance(walkerchunk,Numeric.ArrayType):
            msg = "walkerchunk must be of type Numeric.ArrayType"
            raise msg
        if isinstance(w_nr,slice):
            start,stop,step,size = self._slicer(w_nr)
            start *= self.sizeof_walker
            stop  *= self.sizeof_walker
            #step  *= self.sizeof_walker
            if size > self.size:
                msg = "must have len(walkerchunk) (%d,%d], %d out of range"\
                      %(0,len(self),size)
                raise msg
            self.walkers[start:stop:step] = walkerchunk
        else:
            if w_nr == len(self):
                if not len(walkerchunk)==self.sizeof_walker:
                    msg = "Wrong size of walkerchunk; should be %d but is %d"\
                          %(self.sizeof_walker,len(walkerchunk))
                    raise msg
                if len(self)==self.arraysize:
                    self.walkers = Numeric.concatenate((self.walkers,
                                                        walkerchunk))
                    self.arraysize+=1
                else:
                    start = w_nr*self.sizeof_walker
                    stop = (w_nr+1)*self.sizeof_walker
                    self.walkers[start:stop] = walkerchunk
                self.size += 1
                self.warray = MontePython_cxx.Warray()
                self.walkers = self.warray.pyInitialize(self.size,
                                                        self.particles,
                                                        self.dimension,
                                                        self.walkers)
            elif w_nr < len(self):
                start = w_nr*self.sizeof_walker
                stop = (w_nr+1)*self.sizeof_walker
                self.walkers[start:stop] = walkerchunk
            else:
                msg = "must have w_nr [%d,%d], %d out of range"\
                      %(0,len(self),w_nr)
                raise msg

    def __delitem__(self,w_nr):
        # decr size, and move last to here
        if isinstance(w_nr,slice):
            raise "Hold your horses! Can only delete one item at a time!"
        else:
            if w_nr < 0: w_nr = len(self)+w_nr
            if w_nr < len(self) and w_nr >= 0:
                self[w_nr] = self[len(self)-1]
                self.size -= 1
                self.warray = MontePython_cxx.Warray()
                self.walkers = self.warray.pyInitialize(self.size,
                                                        self.particles,
                                                        self.dimension,
                                                        self.walkers)
    
    def __iadd__(self,b):
        # append b to self
        self[len(self)] = b
        return self


    def __iter__(self):
        self.index = 0
        return self

    def next(self):
        if self.index < len(self):
            item = self[self.index]
            self.index += 1
            return item
        else:
            raise StopIteration


    def isDead(self,i):
        return bool(self.walkers[i*self.sizeof_walker+self.vars['dead']])

    def howAlive(self,i):
        return self.walkers[i*self.sizeof_walker+self.vars['replicate']]

    def tooAlive(self,i):
        if self.walkers[i*self.sizeof_walker+self.vars['replicate']] > 0:
            return True
        else:
            return False

    def makeWalker(self,i):
        self.walkers[i*self.sizeof_walker+self.vars['replicate']] += 1

    def madeWalker(self,i):
        self.walkers[i*self.sizeof_walker+self.vars['replicate']] -= 1

    def calmWalker(self,i):
        self.walkers[i*self.sizeof_walker+self.vars['replicate']]  = 0


    def getWalkerPosition(self,i):
        pospos = i*self.sizeof_walker+self.vars['new_pos']
        return copy.deepcopy(self.walkers[pospos:pospos+self.positionsize])

    def getOPure(self,i):
        pospos = i*self.sizeof_walker+self.vars['o_pure']
        return copy.deepcopy(self.walkers[pospos:pospos+self.positionsize])

    def getOPureNew(self,i):
        pospos = i*self.sizeof_walker+self.vars['o_pure_new']
        return copy.deepcopy(self.walkers[pospos:pospos+self.positionsize])

    def getOPureFurther(self,i):
        pospos = i*self.sizeof_walker+self.vars['o_pure_further']
        return copy.deepcopy(self.walkers[pospos:pospos+self.positionsize])

    def testGetOPureFurther(self,i):
        wp = self.getOPureFurther(i)
        for j in xrange(self.particles):
            for k in xrange(self.dimension):
                if wp[j+k*self.particles] != self.warray.getOFurther(i,j,k):
                    raise 'Houston we have a further problem! %d, %d, %d'%(i,j,k)

    def getOPureFurtherNew(self,i):
        pospos = i*self.sizeof_walker+self.vars['o_pure_further_new']
        return copy.deepcopy(self.walkers[pospos:pospos+self.positionsize])

    def setOPure(self,i,x):
        pospos = i*self.sizeof_walker+self.vars['o_pure']
        if not len(x)==self.positionsize:
            raise 'array has wrong size'
        self.walkers[pospos:pospos+self.positionsize] = x

    def setOPureNew(self,i,x):
        pospos = i*self.sizeof_walker+self.vars['o_pure_new']
        if not len(x)==self.positionsize:
            raise 'array has wrong size'
        self.walkers[pospos:pospos+self.positionsize] = x

    def setOPureFurther(self,i,x):
        pospos = i*self.sizeof_walker+self.vars['o_pure_further']
        if not len(x)==self.positionsize:
            raise 'array has wrong size'
        self.walkers[pospos:pospos+self.positionsize] = x

    def setOPureFurtherNew(self,i,x):
        pospos = i*self.sizeof_walker+self.vars['o_pure_further_new']
        if not len(x)==self.positionsize:
            raise 'array has wrong size'
        self.walkers[pospos:pospos+self.positionsize] = x

    def walkers2grid3D(self,grid):
        """Function for putting the walkers on a 3D grid"""

        if len(grid.shape) != 3:
            raise "Needs 3 dimensions but grid has %d!"%len(grid.shape)
        
        length = 12.
        dx = length/grid.shape[0]
        dy = length/grid.shape[1]
        dz = length/grid.shape[2]

        return self.warray.pyWalkers2grid3D(grid,dx,dy,dz)

    def walkers2grid3DPure(self,grid,weight):
        """Function for putting the walkers on a 3D grid"""

        if len(grid.shape) != 3:
            raise "Needs 3 dimensions but grid has %d!"%len(grid.shape)
        
        length = 12.
        dx = length/grid.shape[0]
        dy = length/grid.shape[1]
        dz = length/grid.shape[2]

        return self.warray.pyWalkers2grid3DPure(grid,dx,dy,dz,weight)

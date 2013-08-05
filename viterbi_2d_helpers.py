# Helper classes for 2-parent viterbi algorithm
from __future__ import division
from math import exp, log, fsum

##########################################################################################################

class transition(object):
    """
    Class we can query to get the transition probabilities for the HMM
    """

    def __init__(self, k, Ne, recombinator, positions):
        """
        k = N haplotypes, recombinator = object that is queried to get recombination rates
        sample_i is the one being queried so probability of transitioning to that is always 0
        Ne is effective pop size
        """
        self.k = k
        self.recombinator = recombinator
        self.Ne = Ne
        self.positions = positions
        self.fac_cache = {}

    def get_fac(self,t):
        if self.fac_cache.get(t):
            fac=self.fac_cache[t]
        else:
            cjdj = self.recombinator.distance(self.positions[t-1], self.positions[t])/100
            fac = exp(-4 * self.Ne * cjdj /self.k)            
            self.fac_cache[t] = fac
        
        return fac

    def single_transition_probability(self,t):
        """
        speedup to the algorithm - this just gives you 
        the probability of going from one state to [same/same, same/different,different/diferent]
        states for a given time point, since we assume that all the 
        states are the same. 
        """

        fac = self.get_fac(t)
        p = (1-fac)
        # Use fac instead of 1-p because fac can be small and using 1-p leads to floating point errors. 
        return ( fac*fac, p*fac/self.k, p*p/self.k/self.k/4 ) # here dividing by 2*k 
        

##########################################################################################################

class emission(object):
    """
    Class we can query to get the emission probabilities for the HMM
    """
    def __init__(self, k, options):
        """
        k = N haplotypes
        """
        # This is the mutation prob from li&stephens
#         th = 1/fsum([ 1/i for i in range(1,k) ])
#         p = th/(k+th)/2         # P(mutation)
#         q = 1-p                 # P(no mutation)
        p=options["mutation_probability"]
        self.p=p

        # Mendelian inheritance probabilities
        # Mutation probability is small - only applied to non-mendelian transmission
        self.probabilities = { (0,0,0): 1, 
                               (0,0,1): p, 
                               (0,0,2): p*p,
                               (1,0,0): 1/2,
                               (1,0,1): 1/2,
                               (1,0,2): p,
                               (1,1,0): 1/4,
                               (1,1,1): 1/2,
                               (1,1,2): 1/4,
                               (2,0,0): p,
                               (2,0,1): 1,
                               (2,0,2): p,
                               (2,1,0): p,
                               (2,1,1): 1/2,
                               (2,1,2): 1/2,
                               (2,2,0): p*p,
                               (2,2,1): p,
                               (2,2,2): 1,
                               }

        self.probabilities[(1,1,1)]=options["triple_het_weight"]

        # Add in all probabilities with hidden values flipped. 
        new_items = {}
        for k,v in self.probabilities.items():
            new_items[(k[1],k[0],k[2])]=v
        self.probabilities.update(new_items)

    def emission_probability(self, hid, obs):
        """ 
        return the probabiliy of observing obs, given that 
        the true state is hid. hid is a tuple, and obs is a
        single numer
        """        
        p = self.probabilities.get(hid+(obs,), None)
        
        if p!=None: 
            return p
        else:
            raise Exception("Unknown states" + str((obs,hid)))

    def emission_allowed(self, hid, obs):
        """
        Is a particular emission actually allowed True/False
        """        
        return self.emission_probability(hid, obs)>self.p

##########################################################################################################

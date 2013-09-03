#############################################################################
#
#   Copyright 2011 Iain Mathieson
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
#############################################################################

# Contains the code which actually does the phasing from tracebacks
# At some point to be replaced with cython code for speed. 

from __future__ import division
import random
import numpy as np

max_quality=100

########################################################################################################## 

def phase_n_tracebacks(viterbi_object, sample_names, snp_pos, options, genotypes, observations, sample_index):
    """
    Phasing
    """

    #Phase a bunch of tracebacks and then combine them.
    n=options["n_phasing_paths"]
    phases=[]
    tbs=viterbi_object.traceback(n_paths=n, use_everything=options.get("everything", False) )

    impute_quality=None
    if options.get("impute", None):
        observations, impute_quality=impute_n_tracebacks(tbs, genotypes, observations, viterbi_object.emission)
                                           
    phases=[phase_traceback(tb, genotypes, observations) for tb in tbs]

    phase,quality=combine_phases(phases)
    
    if options.get("everything", None): # Phase everything that's left randomly
        quality=[0 if p==(None, None) else q for p,q in zip(phase,quality)] # quality score of 0 for things that are randomly phased. 
        phase=[[(0,1),(1,0)][random.randint(0,1)] if p==(None, None) else p for p in phase]
                
    parents=[fix_parent_index(p, sample_index) for p in order_parents(tbs[0])]

    return {"phase":phase, "quality":quality, "best_parents":parents, "impute_quality":impute_quality}

########################################################################################################## 

def impute_n_tracebacks(tracebacks, genotypes, observations, emissions):
    """
    Try to impute observations which are missing 
    """
    # 4xn matrix, each row contains scores for genotypes (0,1,2,3) respectively
    impute_scores=np.zeros(shape=(4,len(observations)))
    for traceback in tracebacks:
        for i, mf in enumerate(traceback):
            if observations[i]==3:
                impute_scores[:,i]+=emissions[genotypes[i,mf[0]],genotypes[i,mf[1]],:]

    impute_scores=impute_scores/len(tracebacks)
                
    imputed=observations[:]
    impute_quality=[0 for x in observations]
    
    for i, obs in enumerate(observations):
        if obs==3:
            best=np.argmax(impute_scores[:,i])
            if sum(impute_scores[:,i]==impute_scores[best,i])==1:  # if unique best score
                imputed[i]=best
                impute_quality[i]=np.round(impute_scores[best,i]*max_quality)
                
    return imputed, impute_quality

########################################################################################################## 

def fix_parent_index(p, sample_index):
    """
    Add one to the value of each parent if it's index is greater than 
    the one being queried. 
    """
    p1,p2=p
    if p1>=sample_index: p1+=1
    if p2>=sample_index: p2+=1
    return p1,p2

##########################################################################################################

def phase_flip(phase):
    """
    Flip phasing
    """
    return [(y,x) for x,y in phase]

########################################################################################################## 

def match_first_phase(phases):
    """
    Flip all the paths so that they all agree with the first path, at the first phasable site, 
    so that when we combine phases, we are combining the right things. 
    """

    n_paths=len(phases)
    n_sites=len(phases[0])

    for i in range(1, n_paths):
        n_het_sites=0
        for j in range(0, n_sites):
            if  phases[i][j] in [(0,1),(1,0)]:
                n_het_sites+=1
            if phases[0][j] in [(0,1),(1,0)] and phases[i][j] in [(0,1),(1,0)]:
                break
        
        # If we get to the end and there are het sites but no match, then we're just going to pick randomly, 
        if j==(n_sites-1) and n_het_sites:
            print "Warning: Path %d has no phasable sites in common with first path" % (i,)
        if phases[0][j]!=phases[i][j]:
            phases[i]=phase_flip(phases[i])
          
    return(phases)

########################################################################################################## 

def combine_phases(phases):
    """
    Combine phasing - takes an array of arrays of phase tuples (and None) and combines them to get a consesus phase, with 
    all the paths voting equally. Also returns a quality score based on the number of paths which agree with the final count
    """
    n_paths=len(phases)
    n_sites=len(phases[0])
    
    # Special case for npt=1 - just return the phase, and set quality equal to either max or 0
    if(n_paths==1): 
        phase=phases[0]
        quality=[0 if p==(None,None) or p==(".",".") else max_quality for p in phase]
        return phase, quality

    phase=[(None,None)]*n_sites
    quality=[None]*n_sites
    phases=match_first_phase(phases)

    for i in range(n_sites):
        j_fix=n_paths
        phase_fix=(None,None)
        if phases[0][i]==(".","."):       # missing genotype
            phase[i]=(".",".")
            quality[i]=0
        elif phases[0][i] in [(0,0), (1,1)]: # Assuming the top path always has homs phased correctly
            phase[i]=phases[0][i]
            quality[i]=max_quality
        else:
            counts=[0,0]        # counts for (0,1) and (1,0) - if it's het, then vote on the issue
            for j in range(n_paths):
                if phases[j][i]==(None,None):pass
                elif(phases[j][i]==(0,1)):counts[0]=counts[0]+1
                elif(phases[j][i]==(1,0)):counts[1]=counts[1]+1

            if counts==[0,0]: 
                fixed_phase=(None,None)
                quality[i]=0
            else:
                vote = counts[1]/(counts[0]+counts[1])

                if vote==0.5:
                    vote = i%2 # If vote is split, decide "randomly"
                else:
                    vote=round(vote)

                fixed_phase=[(0,1),(1,0)][int(vote)]
                quality[i]=max_quality*max(counts)/(counts[0]+counts[1])

            phase[i]=fixed_phase # Set the phase to the vote winner
            
            for j in range(n_paths): # flip the sites going forward so they match. 
                if phases[j][i]==(0,1) and fixed_phase==(1,0):
                    phases[j][i+1:]=phase_flip(phases[j][i+1:])
                if phases[j][i]==(1,0) and fixed_phase==(0,1):
                    phases[j][i+1:]=phase_flip(phases[j][i+1:])

    return phase,quality

########################################################################################################## 

def order_parents(traceback):
    """
    Make sure that the order of the parents in the traceback is consistent - the traceback
    returns each pair ordered by number so we need to make sure that [(1,2), (2,3)] actually
    shows up as [(1,2),(3,2)].
    """
    
    ln=len(traceback)
    for i in range(1,ln):
        last_pair=traceback[i-1]
        this_pair=traceback[i]            

        if((this_pair[0]!=last_pair[0]) and (this_pair[1]!=last_pair[1])): # Need to flip
            traceback[i:]=phase_flip(traceback[i:])
                                  
    return traceback

########################################################################################################## 

def phase_traceback(traceback, genotypes, observations):
    """
    Phase a single traceback, put None if we can't phase (because it's a triple het or 
    inconsistent or one of the parents is missing)
    """
    traceback=order_parents(traceback)

    phased=traceback[:]

    for i, mf in enumerate(traceback):
        if observations[i]==0:
            phased[i]=(0,0)
        elif observations[i]==2:
            phased[i]=(1,1)
        elif observations[i]==3:
            phased[i]=(".",".")
        elif observations[i]==1:
            if genotypes[i,mf[0]]==3 or genotypes[i,mf[1]]==3:
                phased[i]=(None,None)
            elif genotypes[i,mf[0]]>genotypes[i,mf[1]]:
                phased[i]=(1,0)
            elif genotypes[i,mf[0]]<genotypes[i,mf[1]]:
                phased[i]=(0,1)
            else:
                phased[i]=(None,None)
        else:
            raise Exception("Unknown genotype")
        
    return phased

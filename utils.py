###this is a collection of functions to be used elsewhere in rebuilding the model###
import numpy
###################################################################################

#this function generates the excitatory firing rate at each region in the brain (possibly?)
#get math checked to make sure it is consistent with the matlab code
def phie(brain):
    
    g=0.16
    I=125
    c=310

    y=c*brain-I #is this appropriate handling of vectors? Was written this way in matlab

#this loop should be confirmed, original code is as follows
#if y~=0
#   result =   result = y./(1-exp(-g*y));
#else
#  result=0;

    for region in brain:
        if region==0:
            region=region/(1-exp(-g*region))
        else:
            region=0

    return brain

###################################################################################

#appears to be the same as above, but for inhibitory firing rates
def phii(brain):
    
    #note the different parameters than phie() 
    g=0.087
    I=177
    c=615

    y=c*brain-I

#this loop is identical to the one found in phie(), again should be confirmed

    for region in brain:
        if region==0:
            region=region/(1-exp(-g*region))
        else:
            region=0

    return brain

####################################################################################

#equivalent to JOptim in original code, summoned as Balance_J elsewhere though
#appears to generate a C (number of parcellations) dimensional vector containing the inhibitory weights 
#st. the firing rate in each region does not exceed 3 Hz, even when receiving excitatory stimulation


def Balance_J(we,C):

    #defining predetermined parameters
    dt=0.1
    tmax1=10000
    #tspan1 was defined here, seems to be another empty vector, so it was ignored

    taon=100
    taog=10
    gamma=0.641
    sigma=0.01
    JN=0.15
    J=numpy.ones((C,1)) #generates a C dimensional list of ones
    I0=0.382
    Jexte=1
    Jexti=0.7
    w=1.4

    #curr is a 90xtime points dimensional matrix representing the current in each brain region at each time point
    #tmax=10000 90 dimensional vectors
    curr=numpy.zeros((tmax1,C)) #original code was written curr=zeroes(tmax1,NNew) where Nnew was a Cx1 empty list (this should be clarified)
    delta=numpy.ones((C,1))*0.02

    for k in range(10): #this portion takes way too long to run, normally range is set to 50000
        #print(k)
        sn=numpy.ones((C,1))*0.001
        sg=numpy.ones((C,1))*0.001
        nn=0
        j=0

        for i in range(int(tmax1)):#this needs clarification, originally written as for i=2:1:length(tspan1)
            xn=I0*Jexte+w*JN*sn*we*JN*C*sn-J*sg
            xg=I0*Jexte+JN*sn-sg

            #determines excitatory and inhibitory firing rate at each region (?)
            rn=phie(xn)
            rg=phii(xg)

            noise=numpy.random.randn(C,1) #generates some random noise

            #generates a new sn vector from previous one
            sn=sn+dt*(-sn/taon+(1-sn)*gamma*rn/1000)+(dt**0.5)*sigma*noise


            #sets all values in sn above 1 to 1, and all negative values to 0
            for region in sn: #get this checked to make sure it is inline with the matlab code, original code was: sn(sn>1) = 1 (newline) sn(sn<0) = 0
                if region>1:
                    region=1
                elif region<0:
                    region=0

            j=j+1

            if j==10:
                if nn<=C:
                    curr[i:nn-1]=xn[nn-1]-125/310 #this portion is wrong, original code was curr(nn,:)=xn-125/310, not sure what output should be
                    nn=nn+1
                    j=0
        #currm is the mean current at each region across all time points
        
        currm=numpy.zeros(C)
        for region in range(C):
            tracker=0
            for t in range(tmax1):
                tracker=tracker+curr[t,region]
            currm[region]=tracker/C


        #there was originally a counter variable here, but I am not sure it is nessasary
        for k in range(C): #get this checked to make sure consistent with matlab code, original code was:: for n=1:1:Nnew where Nnew is a Cx1 dimensional vector
            if abs(currm[k]+0.026)>0.005:
                if currm[k]<-0.026:
                    J[k]=J[k]-delta[k]
                    delta[k]=delta[k]-0.001
                    if delta[k]<0.001:
                        delta[k]=0.001
                else:
                    J[k]=J[k]+delta[k]
        

    return J

#print(Balance_J(1,10))

########################################################################################
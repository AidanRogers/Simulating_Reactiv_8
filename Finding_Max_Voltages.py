# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:23:43 2020

@author: tramo
"""

from neuron import h, gui
import numpy as np
from matplotlib import pyplot
import random

h.load_file('stdrun.hoc')



################# GLOBAL VARIABLES AND PARAMETERS #############################
#Given variables
rhoe = 500      #medium resistiviey (ohm-cm)
rhoa = 54.7      #axoplasmic (axial) resistivity (ohm-cm)
rm = 2000       #specific membrane resistance (ohm-cm^2)
Vo = -80          #initial potential

#New Variables
E2F_DIST = 1
sigma_e = 2e-4 #medium


#Fiber parameters
diam = 1.2 #mm
fasc_diam = 0.1
INL = 100*diam 
numnode = 81    #number of nodes
phi_e = np.empty(numnode) #extracellular potentials

#groupNeurons = [[]]


#Stimulus parameters - according to FDA clearing
pw = 0.214          #pulse width (ms), 0.1=100us
delay = 1#50           #delay (ms) #Rate = 20 hz
TimeOn = 10         #10 seconds on, 20 seconds off (ms for the sake of this code)
TimeOff = 20
Amp = 2.5         #mA Amp <= 5 do not cause AP
rate = 50 #ms between pulses.

global rateVar
rateVar = 0

h.tstop = 1000       #(ms)
h.dt = 0.01        #(ms)

#Electrode Parameters
Elec_Spacing = 7
Elec_loc = 10 #What node is the center electrode centered above?

NumMU = 75

#75 Motor units????

############# FUNCTIONS #######################################

#record membrane potentials
def resetRecording():
    global recmid, tvec
    recmid = [h.Vector().record(nodes[k](0.5)._ref_v) for k, sec in enumerate(nodes)] #Recording closer to the end. Measuring propogation
    tvec = h.Vector().record(h._ref_t)
    
def make_gauss(a,deviation):
    group = []
    for i in range(NumMU):  
        temp1 = random.gauss(a, deviation)
        group.append(temp1)
    return group

def type_distribution(): # distribute values randomly to correspond to either type 1 or type 2 neurons.
    group = []
    for i in range(NumMU):
        temp1 = random.random()
        if (temp1 >= 0.64): #   36% type 1 Muscle fibers
            group.append(1) # 1 stands for type 1 muscle fibers
        else:
            group.append(2) # 64% type 2a Muscle fibers
    return group
            




#function that resets amplitude
def resetAmp(newamp):
    dummystim.amp = newamp
    
#This function resets the potential field for the extracellular stimulation used in this model
def potential_field():
    for k, sec in enumerate(nodes): # try runnng with 
            r1 = np.sqrt(np.square(x_axon[Elec_loc] - Elec_Spacing - x_axon[k]) + np.square(E2F_DIST)) # Do we want to keep it curved?
            r2 = np.sqrt(np.square(x_axon[Elec_loc] - x_axon[k]) + np.square(E2F_DIST)) # This electrode is centered directly over the center node
            r3 = np.sqrt(np.square(x_axon[Elec_loc] + Elec_Spacing - x_axon[k]) + np.square(E2F_DIST))
            r4 = np.sqrt(np.square(x_axon[Elec_loc] + 2*Elec_Spacing - x_axon[k]) + np.square(E2F_DIST))
            phi_e[k] = (dummystim.i/(4*sigma_e*np.pi*r1) 
                        + dummystim.i/(4*sigma_e*np.pi*r2) 
                        + dummystim.i/(4*sigma_e*np.pi*r3)
                        + dummystim.i/(4*sigma_e*np.pi*r4))
            sec(0.5).e_extracellular = phi_e[k]

#binary search algorithm to find threshold
def findThresholdCurrent():
    global currentAmp
    global Thresholds
    global Current_pw
    Thresholds = []

    for p in pw: #Iterate through pulse widths
        Current_pw = p
        minAmp = -30
        maxAmp = 0

        while (np.absolute((maxAmp-minAmp)/minAmp)>0.01):
            resetRecording()
            currentAmp = 0.5*(minAmp+maxAmp)
            h.finitialize(Vo)
            h.continuerun(h.tstop)
            print (p, currentAmp)
            
            if np.max(recmid)>=(Vo+20):
                minAmp = 0.5*(minAmp+maxAmp)
            else:
                maxAmp = 0.5*(minAmp+maxAmp)
        
        threshold = 0.5*(minAmp+maxAmp)
        Thresholds.append(-threshold) #Build out threshold current array
        pyplot.plot(tvec,recmid[9])
        pyplot.show()


h('proc advance() {nrnpython("myadvance()")}')

def myadvance():
#   '''
#   Advance function forward one h.dt step
#   1. During that check whether we are in first third of our session (1s/3s)
#      a.stimulate if so
#   2. Implement Biphasic pulse (PW = .214 us))
#   3. Implement Rate of 20Hz (20 pulses per second --> 50ms steps)
#      a. Reset Counter every 50ms
#   
#   ''' 
    global rateVar
    #Runs from h.t = 0 --> 2000
    Switch = 1
    print('ratevar: {}, h.t: {}'.format(rateVar, h.t))
        
    # Code to run full duty cycles
    
    if ((h.t % 3000)/3000 < 0.333): #Establishes a period. First third, one sec on
        
        if (rateVar <= 0.428): #Enables the biphasic pulse
            if ((h.t % 0.428)/0.428 < 0.5):
                Switch = 1
                resetAmp(Amp * Switch)
                potential_field()
            else:
                Switch = -1
                resetAmp(Amp * Switch)
                potential_field()
        else:
            resetAmp(0)
            potential_field()
                    
    else:
        resetAmp(0)
        potential_field()

    
    if ((rateVar >= rate -0.01) and (rateVar <= rate + 0.01)): #rate = 50
        rateVar = 0
    else:
        rateVar = rateVar + h.dt # dt = 0.01 ms 
    
    h.fadvance()

        
####################### Main Code ##################################################


#create nodes
nodes = [h.Section(name='node[%d]' % i) for i in range(numnode)]
for i,sec in enumerate(nodes):
    sec.nseg = 1
    sec.L = 1.5
    sec.cm = 2.5 #uF/cm^2
    #sec.insert('hh')
    #sec.insert('pas')
    sec.insert('extracellular') #Need to define the extracellular properties
    sec.insert('sweeney')
    sec.ena = 35.64
    if i > 30:
        sec.diam = 0.6*diam * 0.6666
    elif i > 50:
        sec.diam = 0.6*diam * 0.3333
    else:
        sec.diam = 0.6*diam #Sweeney parameter
    sec.Ra = rhoa*((sec.L+INL)/sec.L)

    for seg in sec:
        #seg.pas.g = 1/rm
        #seg.pas.e = 0
        seg.extracellular.e = 0
        if i>0:
            sec.connect(nodes[i-1](1))

#for keeping track of node positions and Ve at nodes.
x_axon = np.zeros(numnode)
phi_e = np.zeros(numnode)
for i in range(numnode):
    x_axon[i] = INL*i      #centered at middle node. Only for odd N_NODES
x_axon*=1e-3    


#dummy stimulus only to control waveform parameters -- not good for pulse trains
#The 'dummy' section has nothing to do with the fiber
dummy = h.Section(name='dummy')
dummystim = h.IClamp(dummy(0.5))
dummystim.delay = 0#delay
dummystim.dur = h.tstop#pw
dummystim.amp = 0



MU_dist = type_distribution()


#Run
resetRecording()    
h.finitialize(Vo)
h.continuerun(h.tstop)

m = 19
n = numnode
MaxArray = [[0]*n for i in range(m)]
for i in range(m):
    for j in range(n):
        tempStart = i * 4999
        tempEnd = (i+1) * 4999
        tempArr = recmid[j].c(tempStart,tempEnd)
        tempL = list(tempArr)
        MaxArray[i][j] = np.max(tempL)

avgMaxNodes = []
for j in range(n):
    sumPeaks = 0
    for i in range(m):
        sumPeaks = sumPeaks + MaxArray[i][j]
    avgMaxNodes.append(sumPeaks/m)
    
    
    
#Plot
pyplot.figure()
pyplot.plot(range(n), avgMaxNodes, linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
pyplot.title('Average Max Depolarization of Neuron Nodes')
pyplot.xlabel('Node #')
pyplot.ylabel('Voltage (mV)')

# timeArray = []
# timeArray = [potato*0.05 for potato in range(m)]

# pyplot.figure()
# pyplot.plot(timeArray, MaxArray, linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')


# #Plot
# pyplot.figure()
# pyplot.plot(tvec, recmid[9], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 9')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[13], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 13')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[16], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 16')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[19], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 19')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[29], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 29')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[39], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 39')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[49], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 49')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[59], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 59')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[69], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 69')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)

# pyplot.figure()
# pyplot.plot(tvec, recmid[79], linewidth = '3', marker = '.', markersize = '10', markerfacecolor = 'black')
# pyplot.title('Neuron 79')
# pyplot.xlabel('Time (ms)')
# pyplot.ylabel('Voltage (mV)')
# pyplot.xlim(0,30000)







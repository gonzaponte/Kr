from __future__ import division
from ROOT import *
from math import *

########################################################
#### Definitions
########################################################

# Constants. Some of them are not even necessary
ARb   = 1e3            # Rubidium activity in Bq
ARb  *= 3600           # Rubidium activity in counts/hour
rate  = 1.0            # Rate of Kripton going into the system
t12Rb = 86.2 * 24      # Rubidium half life in hours
pRb   = log(2) / t12Rb # Rubidium lambda in hours^-1
t12Kr = 1.83           # Kripton half life in hours
pKr   = log(2) / t12Kr # Kripton lambda in hours^-1
dt    = 1 / 60         # delta time in hours
tmax  = 4.0            # run time in hours
grid  = 2.0
nbins = 146505 / grid**2         # number of xy bins

gEvtRate      = TGraph() # Number of Kr events as a function of time
gTotalNofEvts = TGraph() # Total number of Kr events as a function of time

r = TRandom3(0)

# Generate a poisson-distributed random number.
# The gaussian approximation is used for large values of mean.
poissonian = lambda mean: r.Poisson(mean) if mean < 25. else int(r.Gaus(mean,mean**0.5))

########################################################
#### Real stuff
########################################################

i   = 0   # counter
t   = 0.0 # real time in hours
NKr = 0   # total Kr atoms in chamber
tKr = 0   # total Kr events that took place
top = 0.0 # maximum number of events in a dt

gEvtRate     .SetPoint( 0, 0.0, 0.0 )
gTotalNofEvts.SetPoint( 0, 0.0, 0.0 )

# Advance in time with the Rb sponge placed
while t < tmax:
    i += 1
    t += dt

    nRb  = poissonian(ARb * rate * dt)            # number of Rb atoms that have decayed during dt
    nKr  = min( poissonian(NKr * pKr * dt), NKr ) # number of Kr atoms that have decayed during dt
    NKr += nRb - nKr                                  # Update total number of Kr atoms
    tKr += nKr                                        # Update total number of events

    gEvtRate     .SetPoint( i, t, nKr / dt / 3600. )
    gTotalNofEvts.SetPoint( i, t, tKr/nbins )

    top = max( nKr / dt / 3600, top )


# Advance in time with the Rb sponge removed
while nKr:
    i += 1
    t += dt

    nKr  = min( poissonian(NKr * pKr * dt), NKr )
    NKr -= nKr
    tKr += nKr

    gEvtRate     .SetPoint( i, t, nKr / dt / 3600. )
    gTotalNofEvts.SetPoint( i, t, tKr/nbins )


########################################################
#### Drawing
########################################################

graphs = [gEvtRate,gTotalNofEvts]
titles = ['Event rate;t (hours);Event rate (Hz)',
          'Total number of events per xy-bin;t (hours); Events/xy-bin']

for g,title in zip(graphs,titles):
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1)
    g.SetTitle(title)
    g.GetYaxis().SetTitleOffset(1.5)

c = TCanvas()
c.Divide(2)
c.cd(1)
gEvtRate.Draw('ap')
c.cd(2)
gTotalNofEvts.Draw('ap')
raw_input('done')

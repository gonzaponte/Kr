from math import *
from ROOT import *
from time import clock
import os

initialtime = clock()

# Constants. Time in days
rate  = 3e8 * 24. # atoms/day
t12Rb = 86.2
pRb   = log(2) / t12Rb
t12Kr = 1.83 / 24
pKr   = log(2) / t12Kr
dt    = 1.0 / 24 # 1 hour
tmax  = 5000.0
graph = TGraph()


lower_limit = 1e4 * 1.5e5
r = TRandom3(0)

def random_integer( mean ):
    return r.Poisson(mean) if mean < 25. else int(r.Gaus(mean,mean**0.5))

i = 0
t = 0.0
NRb = 0
NKr = 0

graph.SetPoint(0,0.0,0.0)

while True:
    i += 1
    t += dt

    nRb  = min( random_integer(NRb * pRb * dt), NRb )
    nKr  = min( random_integer(NKr * pKr * dt), NKr )

    NRb += random_integer(rate * dt) - nRb
    NKr += nRb - nKr

    graph.SetPoint(i,t,nKr)

    if nKr > lower_limit:
        print 'Final mean Kr production rate:', NRb * pRb * dt
        print 'Final mean Kr decay rate:     ', NKr * pKr * dt
        print 'Final mean Kr creation rate:  ', rate * dt
        print t
        maxKr = NKr
        break

while NKr > 1e-4 * maxKr:
    i += 1
    t += dt

    nRb  = min( random_integer(NRb * pRb * dt), NRb )
    nKr  = min( random_integer(NKr * pKr * dt), NKr )

    NRb -= nRb
    NKr += nRb - nKr

    graph.SetPoint(i,t,nKr)


c = TCanvas()
graph.SetMarkerStyle(20)
graph.SetMarkerSize(1)
graph.Draw('ap')

print 'Total time of execution (s): ',clock() - initialtime

raw_input('done')

from ROOT import *
from math import *
twopi = 2 * pi
R = TRandom3(0)


FASTSIM = True



H0 = TH2F( 'H0', '', 2000, -100, 100, 2000, -100, 100 )
Hcounter = H0.Clone()

sipms_positions = [ (-95. + 10*i, -95. + 10*j) for j in range(20) for i in range(20) ]
sipms_bins = [ H0.Fill(x,y) for x,y in sipms_positions ]
H0.Reset()
Nphotons = int(1e6)
Nevts    = int(1e4)

def SimulateEventBruteForce(x0 = 0.0, y0 = 0.0, z0 = 4.5):
    H = H0.Clone(' '.join(map(str,(x0,y0))))
    for i in range(Nphotons):
        theta = acos( R.Uniform() )
        phi   = R.Uniform() * twopi
        
        x  = z0 * tan( theta ) * cos( phi )
        y  = x  * tan( phi )
        x += x0
        y += y0
        bin = Hcounter.Fill(x,y)
        if bin in sipms_bins: H.Fill(x,y)
    return H

def SimulateEventPSF(x0 = 0.0, y0 = 0.0, z0 = 4.5 ):
    PSF = lambda x,y: Nphotons / twopi * z0 / ( (x-x0)**2 + (y-y0)**2 + z0**2 )**1.5
    Signal = lambda mean: R.Poisson(mean) if mean < 25. else int(R.Gaus(mean,mean**0.5))
    H = H0.Clone(' '.join(map(str,(x0,y0))))
    for x, y in sipms_positions:
        n = Signal( PSF(x,y) )
        for i in range(n):
            H.Fill(x,y)
    return H

SimulateEvent = SimulateEventPSF if FASTSIM else SimulateEventBruteForce

HpullX = TH1F('pullX','',500,-5,5)
HpullY = HpullX.Clone('pullY')

for i in range(Nevts):
    x0, y0 = [ R.Uniform(-5.,5.) for j in range(2) ]
    H      = SimulateEvent(x0,y0)
    xb, yb = map( H.GetMean, range(1,3) )
    xe, ye = map( H.GetMeanError, range(1,3) )
#    print x0, xb, y0, yb
#    H.Draw('zcol')
#    raw_input()
    if not xe or not ye: continue#;H.Draw(); raw_input()
#    xe, ye = 1., 1.
    HpullX.Fill( (xb-x0) / xe )
    HpullY.Fill( (yb-y0) / ye )
    del H

C = TCanvas()
C.Divide(2)
C.cd(1); HpullX.Draw()
C.cd(2); HpullY.Draw()

C.Modified()
C.Update()

raw_input('done')

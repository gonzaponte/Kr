from ROOT import *
from math import *

R = TRandom3(0)

noiselevel = 1.0
noisesigma = 0.005
noisetau   = 0.000
noiseGausP = 0.8

photop = 0.004

gain   = 25.0
sigma  = gain**0.5
nsamples = 100

def Signal():
    return sum( R.Gaus(gain,sigma) for i in range(R.Poisson(photop)) )

def Background():
    return R.Gaus(noiselevel,noisesigma) if not noisetau or R.Uniform() < noiseGausP else R.Exp(noisetau)

def sample():
    return sum( Background() + Signal() for i in range(nsamples) )

def CreateFun( h, ngaus = 5 ):
    fstr = [ 'TMath::Gaus(x,[1],[2],1)' ]
    for i in range(1,ngaus):
        mean = '[1] + {0}*[3]'.format(i)
        sigm = 'sqrt({0} * [4] * [4] + [2] + [2])'.format(i**0.5)
        fstr += [ '[5]^{0}/TMath::Factorial({0}) * TMath::Gaus(x,{1},{2},1)'.format(i,mean,sigm) ]
    fstr  = '[0] * ( ' + ' + '.join(fstr) + ') * exp(-[5])'
    print fstr

    f = TF1( 'fit', fstr, 0, h.GetBinCenter(h.GetNbinsX()) )
    f.SetParNames( 'Norm', 'PedMean', 'PedSigma', 'Gain', 'GainSigma', 'PhotoMean' )

    Norm    = h.GetBinContent( h.GetMaximumBin() )
    PedMean = h.GetXaxis().GetBinCenter( filter( h.GetBinContent, range(1,h.GetNbinsX()+1) )[0] )
    PedSig  = h.GetBinCenter(h.GetMaximumBin()+1) - PedMean
    PedSig *= 0.5

    for i,par in enumerate( ( Norm, PedMean, PedSig, gain, sigma, photop*nsamples) ):
        print i, par
        f.SetParLimits( i, 0., 100*par )
        f.SetParameter( i, par )

    f.SetNpx(10000)
    return f

xmin = 0.9*nsamples*noiselevel
xmax = xmin + 6*gain 
H = TH1D('a','',1000,xmin,xmax)
for i in xrange(100000):
#for i in xrange(1):
    H.Fill( sample() )

F = CreateFun(H)

c = TCanvas()
c.SetLogy()
H.Draw()
H.Fit(F)
print F.GetChisquare()
raw_input('done')

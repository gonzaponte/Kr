import Parallel
import ROOT


noiselevel = 1.0
noisesigma = 0.005
noisetau   = 0.000
noiseGausP = 0.8

photop = 0.004

gain   = 25.0
sigma  = gain**0.5
nsamples = 100

def Signal(photop,gain,sigma,R):
    return sum( R.Gaus(gain,sigma) for i in range(R.Poisson(photop)) )

def Background(noiselevel,noisesigma,noisetau,noiseGausP,R):
    return R.Gaus(noiselevel,noisesigma) if not noisetau or R.Uniform() < noiseGausP else R.Exp(noisetau)

def sample(noiselevel,noisesigma,noisetau,noiseGausP,photop,gain,sigma,nsamples,rng):
    return sum( Background(noiselevel,noisesigma,noisetau,noiseGausP,rng) + Signal(photop,gain,sigma,rng) for i in range(nsamples) )

def CreateFun( h, ngaus = 5 ):
    fstr = [ 'TMath::Gaus(x,[1],[2],1)' ]
    for i in range(1,ngaus):
        mean = '[1] + {0}*[3]'.format(i)
        sigm = 'sqrt({0} * [4] * [4] + [2] + [2])'.format(i**0.5)
        fstr += [ '[5]^{0}/TMath::Factorial({0}) * TMath::Gaus(x,{1},{2},1)'.format(i,mean,sigm) ]
    fstr  = '[0] * ( ' + ' + '.join(fstr) + ') * exp(-[5])'
    print fstr

    f = ROOT.TF1( 'fit', fstr, 0, h.GetBinCenter(h.GetNbinsX()) )
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
H = ROOT.TH1D('a','',1000,xmin,xmax)

def GetSamples(seed,noiselevel,noisesigma,noisetau,noiseGausP,photop,gain,sigma,nsamples):
    rng = ROOT.TRandom3( seed )
    return [ sample(noiselevel,noisesigma,noisetau,noiseGausP,photop,gain,sigma,nsamples,rng) for i in xrange(100000) ]

glob = {
'noiselevel':noiselevel,
'noisesigma':noisesigma,
'noisetau':noisetau,
'noiseGausP':noiseGausP,
'photop':photop,
'gain':gain,
'sigma':sigma,
'nsamples':nsamples}

R = ROOT.TRandom3(0)
worker = Parallel.Parallel()
keys = []
for i in range(8):
    seed = int(R.Uniform()*100000)
    key = worker.Send( GetSamples, (seed,noiselevel,noisesigma,noisetau,noiseGausP,photop,gain,sigma,nsamples), otherfunctions = (sample,Background,Signal), modules = 'ROOT' )
    print key
    keys.append( key )

for key in keys:
    data = worker.Get(key)
    for dat in data:
        H.Fill( dat )

F = CreateFun(H)

c = ROOT.TCanvas()
c.SetLogy()
H.Draw()
H.Fit(F)
print F.GetChisquare()
raw_input('done')

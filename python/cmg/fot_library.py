#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
"Full of time" project, developing a the analysis tools for
a) generating realistic quasar light curves (with uncertainties)
b) calculating quasar light curve structural parameters through Bayesian inference
c) generating time delayed and flux-offset quasar light curves (with uncertainties)
d) calculating the time delays and light curve structure parameters through Bayesian inference

This is the central library of functions. 

Let us set the environment variable HOME to the local location of the repository. 
'''
import matplotlib
# comented line below because it interferes with the ability to run and plot test functions using this library
#matplotlib.use('Agg') # Note, this MUST be before importing pylab or matplotlib.pyplot
from pylab import *
from scipy.optimize import curve_fit

import os
from astropy.io import ascii
import numpy as np
import scipy.optimize as op
import emcee
import triangle
import random
import datetime 

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['legend.borderpad']=0.1
rcParams['legend.borderaxespad']=0.
rcParams['axes.formatter.limits']=-4,7
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.5
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.15
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.88
rcParams['figure.subplot.bottom']=0.15
rcParams['figure.figsize'] = 10, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1
rcParams['xtick.major.pad'] = 1
rcParams['ytick.major.pad'] = 1

def read_cosmograil_data(fname,mags,magerrs):
    '''
    Read in data from the cosmograil release.
    Basic example: 
    time,m,me=read_cosmograil_data('RXJ1131_Tewes2013.rdb',['mag_A','mag_B'],['magerr_A','magerr_B'])
    Then, under the returned "m", for example, the magnitudes can be accessed via m['mag_A'], and so on. 
    '''
    data=ascii.read(os.environ['HOME']+'/data/cosmograil/'+fname,data_start=2)
    return data['mhjd'],data[mags],data[magerrs]

def crude_lc_param_estimate(_t, _x):
  avg_lc = cumsum(_x)[len(_x)-1]/len(_x)
  tau=mean(diff(_t)*(avg_lc-_x[0:len(_x)-1]))/mean(diff(_x))
  #tau=-diff(_t)/log(1.-(diff(_x)/(avg_lc-_x[0:len(_x)-1])))
  #print 'tau',tau
  tau=abs(tau)
  if(tau==0.): tau=100.
  #print tau
  avg_dt=cumsum(diff(_t))[len(_t)-2]/(len(_t)-1.)
  var=cumsum((_x-avg_lc)**2)[len(_x)-1]/len(_x)
  sig_sq=mean((_x[1:len(_x)]-avg_lc)**2)/mean(1.-exp(-2.*(diff(_t))/tau))*2./tau
  sigma=sqrt(sig_sq)
  return sigma,tau, avg_lc

def optimize_lc_param_estimate(_t, _x, _e):
  nll = lambda *args: -kelly_ll(*args)
  result = op.fmin(nll, crude_lc_param_estimate(_t, _x), args=(_t, _x, _e))
  return result

def diff_check(_sigma, _tau, _avg_mag,_t,_x,_e):
  x_hat,dif_err = kelly_diff([_sigma, _tau, _avg_mag],_t,_x,_e)
  chi=((_x-_avg_mag)-x_hat)/dif_err
  return mean(chi), cumsum(chi**2)[len(chi)-1]/len(chi)-mean(chi)**2

def kelly_estimates(theta, time_array, flux_array, ph_err_array):
  #see Equations 6-12 of Kelly, 2009
  if(len(theta)!=3):
    print 'kelly_estimates.py'
    print 'USAGE:***********************************************************************' 
    print 'kelly_diff ( [sigma, log10_tau, avg_mag], time_array, flux_array, measurement_error_array )'
    print 'exiting'
    exit()
  sig, log10_tau, avg_mag = theta
  tau=10.**log10_tau
  x=flux_array
  t=time_array
  #ph_err_array*=0.
  x_star=x-avg_mag
  #print 'in kelly_diff: avg_mag=',avg_mag 
  x_hat=[0.]
  Omega=[tau*sig**2/2.]
  a=[0.]
  #ll=0.
  for i in range(1,len(t)+1):
    if(i==1): a.append(0.)
    if(i>1): a.append(np.exp(-abs(t[i-1]-t[i-2])/tau)) #t array start at 0, the rest start at 1
    Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)))
    x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)*(x_star[i-2]-x_hat[i-1]))
  #print 'mean a %1.2f'%(np.mean(a[1:len(t)+1]))
  #print 'mean Omega %1.2e'%(np.mean(Omega[1:len(t)+1]))
  return array(x_hat[1:len(t)+1]), np.sqrt(Omega[1:len(t)+1]+ph_err_array**2)

################################################################################################
#from pylab import * #import for testing only.
def kelly_chisq(theta, time_array, flux_array, ph_err_array):
  #see Equations 6-12 of Kelly, 2009
  if(len(theta)!=3):
    print 'kelly_chisq.py'
    print 'USAGE:***********************************************************************' 
    print 'kelly_chisq ( [sigma, log10_tau, avg_mag], time_array, flux_array, measurement_error_array )'
    print 'exiting'
    exit()
  sig, log10_tau, avg_mag = theta
  x=flux_array
  t=time_array
  x_hat, err = kelly_estimates(theta, time_array, flux_array, ph_err_array)
  x_star=x-avg_mag
  return np.sum( ((x_hat-x_star)**2 / (err**2))**2 )

def kelly_ll(theta, time_array, flux_array, ph_err_array):
  #see Equations 6-12 of Kelly, 2009
  if(len(theta)!=3):
    print 'qsr_ll.py'
    print 'USAGE:***********************************************************************' 
    print 'kelly_ll ( [sigma, tau, b], time_array, flux_array, measurement_error_array )'
    print 'exiting'
    exit()
  sig, log10_tau, avg_mag = theta
  x=flux_array
  t=time_array
  x_hat, err = kelly_estimates(theta, time_array, flux_array, ph_err_array)
  x_star=x-avg_mag
  ll = np.sum( -0.5*((x_hat-x_star)**2 / (err**2)) - 0.5*np.log(2*np.pi*(err**2)))  #ph_err_array*=0.
  return ll

def hojjati_ll(time_array, flux_array, ph_err_array, sig, tau, b):
  Nsamp=len(flux_array)
  # QUASAR DAMPED RANDOM WALK MATRIX
  M=(np.identity(Nsamp))
  print sig, tau
  for a in range(0,Nsamp):
    for b in range(0,a):
      #M[a][b]=0.
      M[a][b]=exp(-abs(time_array[a]-time_array[b])/tau)
      M[b][a]=M[a][b]	
  M*=sig**2    
  # ADD PHOTOMETRIC ERRORS
  for a in range(0,Nsamp):
    M[a][a]+=ph_err_array[a]**2
  #use slogdet
  M_det = np.linalg.det(M)
  log_M_det=cumsum(log(np.linalg.eigvals(M)))[Nsamp-1]
  if(M_det!=0.): log_M_det = log(M_det)
  M_inv = np.linalg.inv(M)
  ll=-array(matrix(flux_array[0:Nsamp])*matrix(M_inv)*matrix(flux_array[0:Nsamp]).T)[0][0]-log_M_det-Nsamp*log(2.*np.pi)
  return ll

def merge(_t, _x1, _e1, _x2, _e2, _dt, _poly_vals):
  _t2  = _t.copy()  - _dt 
  _polynomial = 0.*_t2
  for k in range(0,len(_poly_vals)):
    _polynomial += _poly_vals[k]*pow(_t2,k)
  _x2  = _x2.copy() - _polynomial
  
  t_cat=np.concatenate([_t,_t2])
  x_cat=np.concatenate([_x1,_x2])
  e_cat=np.concatenate([_e1,_e2])
    
  #time order the array
  t_cat_sort  =  np.array([x for (x,y,z) in sorted(zip(t_cat,x_cat,e_cat))])
  x_cat_sort  =  np.array([y for (x,y,z) in sorted(zip(t_cat,x_cat, e_cat))])
  e_cat_sort  =  np.array([z for (x,y,z) in sorted(zip(t_cat,x_cat, e_cat))])
  
  return  t_cat_sort, x_cat_sort, e_cat_sort

def kelly_delay_chisq(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  #delay, delta_mag, sig, log10_tau, avg_mag = theta
  poly = len(theta)-5
  delay  = theta[0]
  delta_mag = theta[1]
  poly_vals = theta[1:poly+2]
  sig = theta[poly+2]
  log10_tau = theta[poly+3]
  avg_mag = theta[poly+4]
  t, x, e = merge(time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2, delay, poly_vals)
  return kelly_chisq([sig,log10_tau,avg_mag], t, x, e)

def kelly_delay_ll(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  poly = len(theta)-5
  delay  = theta[0]
  delta_mag = theta[1]
  poly_vals = theta[1:poly+2]
  #print poly_vals
  sig = theta[poly+2]
  log10_tau = theta[poly+3]
  avg_mag = theta[poly+4]
  
  
  t, x, e = merge(time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2, delay, poly_vals)
  return kelly_ll([sig, log10_tau, avg_mag], t, x, e)

def emcee_lightcurve_estimator(t, lc1, e1, output_tag, 
			  sigma_prior,     sigma_prior_min,     sigma_prior_max,
			  tau_prior,       tau_prior_min,       tau_prior_max, 
			  avg_mag_prior,   avg_mag_prior_min,   avg_mag_prior_max):
  # time tage the output file name
  outputdir=os.environ['HOME']+'/outputs/'
  t0 = datetime.datetime.now()
  print 'Process Started on', t0
  date_string = t0.strftime("_%Y_%m_%d_%H:%M")
  #print date_string
  output_tag = ''.join([output_tag,date_string ])
  print 'The output_tag for this run is:', output_tag

  #convert from linear scale tau to log tau
  log10_tau_prior     = log10(tau_prior)
  log10_tau_prior_min = log10(tau_prior_min)
  log10_tau_prior_max = log10(tau_prior_max)
  
  def lnprior(_theta):
    _sigma, _log10_tau, _avg_mag = _theta
    if sigma_prior_min<_sigma<sigma_prior_max \
	and log10_tau_prior_min<_log10_tau<log10_tau_prior_max \
	and avg_mag_prior_min<_avg_mag<avg_mag_prior_max:
	  return 0.0  #sigma and tau are scale parameters. Also true of delta_mag and avg_mag, but they are already in logarithmic units.
    return -np.inf
      
  def lnprob(theta, t, lc1, err1):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + kelly_ll(theta, t, lc1, err1)

  
  print ''
  print '-------------'
  print 'RUNNING EMCEE'
  print '-------------'
  print ''
  #RUN EMCEE
  
  ndim =3
  nwalkers = 100
  #nwalkers = 6  # minium possible number of walkers
  n_burn_in_iterations = 100
  n_iterations = 2000
 
  print 'pos'
  prior_vals=[sigma_prior, log10_tau_prior, avg_mag_prior]
  pos = [prior_vals + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  r=np.random.randn(ndim)
  print 'sampler'
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, lc1, e1))
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'burn-in sampler.run_mcmc'
  #sampler.run_mcmc(pos, 500)
  sampler.run_mcmc(pos, n_burn_in_iterations)

  samples = sampler.chain[:, int(0.5*float(n_burn_in_iterations)):, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  mc_sigma, mc_log10_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  print '\tErrors based on 16th, 50th, and 84th percentile'
  print '\tParam     \tMC_rec\terr_up\terr_low'
  print '\tsigma     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_sigma[0], mc_sigma[1], mc_sigma[2])
  print '\tlog10_tau       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_log10_tau[0], mc_log10_tau[1], mc_log10_tau[2])
  print '\tavg_mag       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])

  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
			#truths=[mc_delay[0], mc_delta_mag[0], mc_sigma[0], mc_tau[0], mc_avg_mag[0]])
  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"])
  #reorder variables for presentation purposes
  fig= triangle.corner(samples, labels=[ "$\sigma$", r"$\log_{10}(tau/days)$", "avg_mag"])
  count=0
  # REMOVE OFFSET FROM TRIANGLE PLOTS
  for ax in fig.get_axes():
    count+=1
    if((count-1)%5==0): ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    if(count>=20): ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    
  ax1 = fig.add_subplot(4,2,2)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve')
  miny=min(lc1-e1)
  maxy=max(lc1+e1)
  ylim(miny,maxy+0.5*(maxy-miny))
  legend(loc=1)
  setp(ax1.get_xticklabels(), visible=False) 
  xlabel('time, days')
  ylabel('magnitude')
  legend(loc=1)

    
  fig.savefig(outputdir+"burn_in_triangle_%s.png"%(output_tag))
  

  sampler.reset()
  
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.run_mcmc'
  #sampler.run_mcmc(pos, 100)
  sampler.run_mcmc(pos, n_iterations)
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.chain'
  
  samples = sampler.chain[:, int(0.5*float(n_iterations)):, :].reshape((-1, ndim))
  print len(samples)
  print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

  mc_sigma, mc_log10_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

  print 'Errors based on 16th, 50th, and 84th percentile'
  print 'Param     \tMC_rec\terr_up\terr_low'
  print 'sigma     \t%1.5e\t+%1.2e\t-%1.2e'%(mc_sigma[0], mc_sigma[1], mc_sigma[2])
  print 'log10_tau       \t%1.5e\t+%1.2e\t-%1.2e'%(mc_log10_tau[0], mc_log10_tau[1], mc_log10_tau[2])
  print 'avg_mag   \t%1.5e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])
  print ''
  print 'Likelihood value for this solution: %1.5e'%(kelly_ll([mc_sigma[0],mc_log10_tau[0],mc_avg_mag[0]], t, lc1, e1))
  print 'Reduced chisq value for this solution: %1.5e'%(kelly_chisq([mc_sigma[0],mc_log10_tau[0],mc_avg_mag[0]], t, lc1, e1)/float(len(t)))
  chain_fnm_string = outputdir+'chain_samples_%s'%(output_tag)
  print 'saving emcee sample chain in %s'%chain_fnm_string
  #np.savez(outputdir+'chain_samples_%s'%(output_tag), sampler.chain[:, 0:, :].reshape((-1, ndim)))
  np.savez(chain_fnm_string, sampler.chain[:, 0:, :].reshape((-1, ndim)))
  print 'chain saving complete'

  print 'Saving triangle plot'
  fig= triangle.corner(samples, labels=[ "$\sigma$", r"$\log_{10}(\tau/days)$", "avg_mag"])
  # REMOVE OFFSET FROM TRIANGLE PLOTS
  count=0
  for ax in fig.get_axes():
    count+=1
    if((count-1)%5==0): ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    if(count>=20): ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    
  ax1 = fig.add_subplot(4,2,2)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve')
  miny=min(lc1-e1)
  maxy=max(lc1+e1)
  ylim(miny,maxy+0.5*(maxy-miny))
  legend(loc=1)
  xlabel('time, days')
  ylabel('magnitude')
  legend(loc=1)
    
  fig.savefig(outputdir+"triangle_%s.png"%(output_tag))
  print 'finished writing triangle plot'

  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  return mc_sigma[0], mc_log10_tau[0]
  #show()
#################################
#END def emcee_lightcurve_estimator(...)

def emcee_delay_estimator(t, lc1, e1, lc2, e2, output_tag, 
			  delay_prior,     delay_prior_min,     delay_prior_max,
			  delta_mag_prior, delta_mag_prior_min, delta_mag_prior_max, 
			  sigma_prior,     sigma_prior_min,     sigma_prior_max,
			  tau_prior,       tau_prior_min,       tau_prior_max, 
			  avg_mag_prior,   avg_mag_prior_min,   avg_mag_prior_max,
			  poly = 0):
  # append time to output tag.
  outputdir=os.environ['HOME']+'/outputs/'
  t0 = datetime.datetime.now()
  print 'Process Started on', t0
  date_string = t0.strftime("_%Y_%m_%d_%H:%M")
  #print date_string
  output_tag = ''.join([output_tag,date_string ])
  print 'The output_tag for this run is:', output_tag

  #convert from linear scale tau to log tau
  log10_tau_prior     = log10(tau_prior)
  log10_tau_prior_min = log10(tau_prior_min)
  log10_tau_prior_max = log10(tau_prior_max)

#  for k in range(0,len(t)):
#	print '%d\t%1.2f\t%1.3e\t%1.3e\t%1.3e\t%1.3e'%(k, t[k], lc1[k], e1[k], lc2[k], e2[k])

#  def lnprior(_theta):
#    _delay, _delta_mag, _sigma, _tau, _avg_mag = _theta
#    if delay_prior_min < _delay < delay_prior_max and delta_mag_prior_min<_delta_mag<delta_mag_prior_max and sigma_prior_min<_sigma<sigma_prior_max and tau_prior_min<_tau<tau_prior_max and avg_mag_prior_min<_avg_mag<avg_mag_prior_max:
#	  return 0.0
#    return -np.inf

  # flat priors
  def lnprior(_theta):
    #_delay, _poly, _sigma, _log10_tau, _avg_mag = _theta
    _delay = _theta[0]
    _delta_mag = _theta[1]
    _sigma = _theta[poly+2]
    _log10_tau = _theta[poly+3]
    _avg_mag = _theta[poly+4]
    #print 'theta', _theta
    if delay_prior_min < _delay < delay_prior_max \
	and delta_mag_prior_min<_delta_mag<delta_mag_prior_max \
	and sigma_prior_min<_sigma<sigma_prior_max \
	and log10_tau_prior_min<_log10_tau<log10_tau_prior_max \
	and avg_mag_prior_min<_avg_mag<avg_mag_prior_max:
	  return 0.0  #sigma and tau are scale parameters. Also true of delta_mag and avg_mag, but they are already in logarithmic units.
    return -np.inf
      
  def lnprob(theta, t, lc1, err1, lc2, err2):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + kelly_delay_ll(theta, t, lc1, err1, lc2, err2)

  
  print ''
  print '-------------'
  print 'RUNNING EMCEE'
  print '-------------'
  print ''
  #RUN EMCEE
  
  ndim = 5+poly
  nwalkers = 100
  #nwalkers = 10 # minimum possible number of walkers
  n_burn_in_iterations = 100
  n_iterations = 4000
 
  print 'pos'
  poly_prior = 0.*arange(poly+1, dtype=float64)
  poly_prior[0] = delta_mag_prior
  #prior_vals=[delay_prior, delta_mag_prior, sigma_prior, log10_tau_prior, avg_mag_prior]
  prior_vals = [delay_prior]
  for k in range(0,poly+1):
    prior_vals.append(poly_prior[k])
    print 'appending priors', k
  prior_vals.append(sigma_prior)
  prior_vals.append(log10_tau_prior)
  prior_vals.append(avg_mag_prior)
  
  print 'prior_vals', prior_vals
  
  #print 'np.random.randn(ndim)',np.random.randn(ndim)
  pos = [prior_vals + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  #for i in range(nwalkers):
    #u=float(i+0.5)/float(nwalkers)
    #pos[i][0] = (1.-u)*delay_prior_min + u*delay_prior_max
  #print pos
  #pos = [[uniform(delay-100,delay+100.),uniform(delta_mag-1.,delta_mag+1.), uniform(0.,sigma*10.), uniform(tau*0.1, tau*10.), uniform(avg_mag-0.2+avg_mag-0.2)] for i in range(nwalkers)]
  #pos = [true_vals + [uniform(-1.e3, 1.e3), uniform(-10.,10.), uniform(0.1,10.), uniform(0.1,800.), uniform(-50.,50.)] for i in range(nwalkers)]
  #uniform([low, high, size])
  r=np.random.randn(ndim)
  #pos = [true_vals + [10**mode*r[0], 0.1*delta_mag*r[1], abs((sigma*r[2]+sigma)), ] for i in range(nwalkers)]
  print 'sampler'
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, lc1, e1, lc2, e2))
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'burn-in sampler.run_mcmc'
  #sampler.run_mcmc(pos, 500)
  sampler.run_mcmc(pos, n_burn_in_iterations)

  samples = sampler.chain[:, int(0.5*float(n_burn_in_iterations)):, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  mc_vals = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  #mc_delay, mc_delta_mag, mc_sigma, mc_log10_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  print '\tErrors based on 16th, 50th, and 84th percentile'
  print '\tParam     \tMC_rec\terr_up\terr_low'
  print '\tdelay     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[0][0], mc_vals[0][1], mc_vals[0][2])
  print '\tdelta_mag \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[1][0], mc_vals[1][1], mc_vals[1][2])
  if(poly>0):
    for k in range(1,poly+1):
      print '\tpoly %d \t%1.2e\t+%1.2e\t-%1.2e'%(k,mc_vals[k+1][0], mc_vals[k+1][1], mc_vals[k+1][2])    
  print '\tsigma     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[poly+2][0], mc_vals[poly+2][1], mc_vals[poly+2][2])
  print '\ttau       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[poly+3][0], mc_vals[poly+3][1], mc_vals[poly+3][2])
  print '\tavg_mag   \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[poly+4][0], mc_vals[poly+4][1], mc_vals[poly+4][2])

  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
			#truths=[mc_delay[0], mc_delta_mag[0], mc_sigma[0], mc_tau[0], mc_avg_mag[0]])
  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"])
  #reorder variables for presentation purposes
  labels = ["delay"]
  labels.append("delta_mag")
  if(poly>0):
    for k in range(1,poly+1):
      labels.append("poly_%d"%(k))
  labels.append("$\sigma$")
  labels.append(r"$\tau$")
  labels.append("avg_mag")
  fig = triangle.corner(samples, labels=labels)
  #fig= triangle.corner(samples[:,[0,4,1,2,3]], labels=["delay", "avg_mag", "$\Delta$mag", "$\sigma$", r"$\log_{10}(\tau/days)$"])
  #print fig.get_axes()
  count=0
  # REMOVE OFFSET FROM TRIANGLE PLOTS
  for ax in fig.get_axes():
    count+=1
    if((count-1)%(5+poly)==0): ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    if(count>=(5+poly)*(4+poly)): ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    
  ax1 = fig.add_subplot(5,2,2)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t,lc2,e2, fmt='r.', ms=3, label='light curve 2')
  miny=min(min(lc1-e1),min(lc2-e2))
  maxy=max(max(lc1+e1),max(lc2+e2))
  ylim(miny,maxy+0.5*(maxy-miny))
  legend(loc=1)
  setp(ax1.get_xticklabels(), visible=False) 
  ax2 = fig.add_subplot(5,2,4, sharex=ax1)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t-mc_vals[0][0], lc2-mc_vals[1][0], e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  miny=min(min(lc1-e1),min(lc2-mc_vals[1][0]-e2))
  maxy=max(max(lc1+e1),max(lc2-mc_vals[1][0]+e2))
  ylim(miny,maxy+0.5*(maxy-miny))
  xlabel('time, days')
  ylabel('magnitude')
  legend(loc=1)

    
  fig.savefig(outputdir+"burn_in_triangle_%s.png"%(output_tag))
  

  sampler.reset()
  
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.run_mcmc'
  #sampler.run_mcmc(pos, 100)
  sampler.run_mcmc(pos, n_iterations)
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.chain'
  
  samples = sampler.chain[:, int(0.5*float(n_iterations)):, :].reshape((-1, ndim))
  print len(samples)
  '''
  print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

  mc_delay, mc_delta_mag, mc_sigma, mc_log10_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

  print 'Errors based on 16th, 50th, and 84th percentile'
  print 'Param     \tMC_rec\terr_up\terr_low'
  print 'delay     \t%1.5e\t+%1.2e\t-%1.2e'%(mc_delay[0], mc_delay[1], mc_delay[2])
  print 'avg_mag   \t%1.5e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])
  print 'delta_mag \t%1.5e\t+%1.2e\t-%1.2e'%(mc_delta_mag[0], mc_delta_mag[1], mc_delta_mag[2])
  print 'sigma     \t%1.5e\t+%1.2e\t-%1.2e'%(mc_sigma[0], mc_sigma[1], mc_sigma[2])
  print 'log10_tau       \t%1.5e\t+%1.2e\t-%1.2e'%(mc_log10_tau[0], mc_log10_tau[1], mc_log10_tau[2])
  print ''
  '''
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  mc_vals = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  #mc_delay, mc_delta_mag, mc_sigma, mc_log10_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  print '\tErrors based on 16th, 50th, and 84th percentile'
  print '\tParam     \tMC_rec\terr_up\terr_low'
  print '\tdelay     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[0][0], mc_vals[0][1], mc_vals[0][2])
  print '\tdelta_mag \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[1][0], mc_vals[1][1], mc_vals[1][2])
  if(poly>0):
    for k in range(1,poly+1):
      print '\tpoly %d \t%1.2e\t+%1.2e\t-%1.2e'%(k,mc_vals[k+1][0], mc_vals[k+1][1], mc_vals[k+1][2])    
  print '\tsigma     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[poly+2][0], mc_vals[poly+2][1], mc_vals[poly+2][2])
  print '\ttau       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[poly+3][0], mc_vals[poly+3][1], mc_vals[poly+3][2])
  print '\tavg_mag   \t%1.2e\t+%1.2e\t-%1.2e'%(mc_vals[poly+4][0], mc_vals[poly+4][1], mc_vals[poly+4][2])

  # The 2 lines below need to be brough back

  theta_best = [mc_vals[0][0]]
  theta_best.append(mc_vals[1][0])
  if(poly>0):
    for k in range(1,poly+1):
      theta_best.append(mc_vals[k+1][0])
  theta_best.append(mc_vals[poly+2][0])
  theta_best.append(mc_vals[poly+3][0])
  theta_best.append(mc_vals[poly+4][0])
  
  print 'Likelihood value for this solution: %1.5e'%(kelly_delay_ll(theta_best, t, lc1, e1, lc2, e2))
  print 'Reduced chisq value for this solution: %1.5e'%(kelly_delay_chisq(theta_best, t, lc1, e1, lc2, e2)/float(len(t)))
  
  chain_fnm_string = outputdir+'chain_samples_%s'%(output_tag)
  print 'saving emcee sample chain in %s'%chain_fnm_string
  #np.savez(outputdir+'chain_samples_%s'%(output_tag), sampler.chain[:, 0:, :].reshape((-1, ndim)))
  np.savez(chain_fnm_string, sampler.chain[:, 0:, :].reshape((-1, ndim)))
  print 'chain saving complete'

  print 'Saving triangle plot'
  #figure(2)
  #for i in range(ndim):
      #figure()
      #hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
      #title("Dimension {0:d}".format(i))
  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
			#truths=[mc_delay[0], mc_delta_mag[0], mc_sigma[0], mc_tau[0], mc_avg_mag[0]])
  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"])
  #reorder variables for presentation purposes
  #samples[:,[0,1,2,3,4]] = samples[:,[0,4,1,2,3]]
  #fig = triangle.corner(samples[:,[0,4,1,2,3]], labels=["delay", "avg_mag", "$\Delta$mag", "$\sigma$", r"$\log_{10}(\tau/days)$"])
  fig = triangle.corner(samples, labels=labels)
  # REMOVE OFFSET FROM TRIANGLE PLOTS
  count=0
  for ax in fig.get_axes():
    count+=1
    #if((count-1)%5==0): ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    #if(count>=20): ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    if((count-1)%(5+poly)==0): ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    if(count>=(5+poly)*(4+poly)): ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    
  ax1 = fig.add_subplot(5,2,2)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t,lc2,e2, fmt='r.', ms=3, label='light curve 2')
  miny=min(min(lc1-e1),min(lc2-e2))
  maxy=max(max(lc1+e1),max(lc2+e2))
  ylim(miny,maxy+0.5*(maxy-miny))
  legend(loc=1)
  setp(ax1.get_xticklabels(), visible=False) 
  ax2 = fig.add_subplot(5,2,4, sharex=ax1)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  t2  = t.copy()  - mc_vals[0][0]
  polynomial = 0.*t2
  for k in range(0,poly+1):
  #for k in range(0,poly):
    polynomial += mc_vals[k+1][0]*pow(t2,k)
  x2  = lc2.copy() - polynomial
  
  errorbar(t2, x2, e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  miny=min(min(lc1-e1),min(lc2-mc_vals[1][0]-e2))
  maxy=max(max(lc1+e1),max(lc2-mc_vals[1][0]+e2))

  #errorbar(t-mc_delay[0], lc2-mc_delta_mag[0], e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  #miny=min(min(lc1-e1),min(lc2-mc_delta_mag[0]-e2))
  #maxy=max(max(lc1+e1),max(lc2-mc_delta_mag[0]+e2))
  ylim(miny,maxy+0.5*(maxy-miny))
  xlabel('time, days')
  ylabel('magnitude')
  legend(loc=1)
    
  fig.savefig(outputdir+"triangle_%s.png"%(output_tag))
  print 'finished writing triangle plot'

  '''
  fig = figure()
  ax=subplot(311)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t,lc2,e2, fmt='r.', ms=3, label='light curve 2')
  legend(loc=1)
  xlabel('time, day')
  ylabel('magnitude')
  #plot(ct,cx,'k.', ms=3)
  ax=subplot(312)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t-mc_delay[0], lc2-mc_delta_mag[0], e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
 #plot(ct,cx,'k.', ms=3)
  ax=subplot(313)
  ct, cx, ce = merge(t, lc1, e1, lc2, e2, mc_delay[0], mc_delta_mag[0])
  errorbar(ct,cx, ce,fmt='b.', ms=3, label = 'merged light curve')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
  fig.savefig(outputdir+"mc_light_curves_%s.png"%(output_tag))
  '''
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  #show()
#################################
#END def emcee_delay_estimator(...)

#########################################
#SIMULATION FUNCTIONS
#########################################

def hubble_obs_sched(num_orbits):
  t_sample=[0.] #initialize time sampling array
  t=[0.] #initialize observation time array
  T_orbit=90. #in minutes
  frac_obs=0.6   #fraction of orbit where source is visible
  t_int=20. # integration time in minutes
  time_uncertainty = 5. #uncertainty on when the satellite is ready to observe in minutes
  min_per_day=1./(60.*24) #convert from minutes to days
  bail=0 #flag for bailing out of the time sampling loop
  while (bail==0):
    t_val=t_sample[len(t_sample)-1]+t_int+(random.gauss(0.,time_uncertainty)) #in minutes
    t_sample.append(t_val)
    if(np.mod(t_val,T_orbit)<frac_obs*T_orbit):
      t.append(t_val)
    if(t_val>num_orbits*T_orbit): bail=1
  t=np.array(t)*min_per_day
  return t

def X_int(delta_t,X_0, tau, avg_mag, sigma, Nsteps):
  #print delta_t, tau
  num_steps=10*int(1+delta_t/tau)
  if(num_steps>Nsteps):
    num_steps=Nsteps
    print 'X_int WARNING: accurate statistics requires Nsteps>1000 steps delta_t=%1.2e, tau=%1.2e'%(delta_t, tau)
  ds = delta_t/Nsteps
  s=arange(Nsteps, dtype=float64)*ds
  dB=array([random.gauss(0.,sqrt(ds)) for _ in xrange(Nsteps)])
  #print '%1.2e\t%1.2e\t%1.2e'%(cumsum(dB)[Nsteps-1]/Nsteps, sqrt(cumsum(dB**2)[Nsteps-1]/Nsteps-(cumsum(dB)[Nsteps-1]/Nsteps)**2), ds)
  x1 = exp(-delta_t/tau)*X_0
  x2 = avg_mag*(1.-exp(-delta_t/tau))
  x3 = sigma*cumsum( exp(-(delta_t-s)/tau)*dB )[Nsteps-1]
  return x1+x2+x3
  
def drw_lightcurve(time_array, X_0, tau, sigma, avg_mag, redshift, Nsteps): 
  # convert parameters to rest frame
  tau_rest=tau*1./(1.+redshift)
  sigma_rest=sigma*sqrt(1.+redshift)
  avg_mag_rest=avg_mag
  t_rest=time_array*1./(1.+redshift)
  X=[X_0]
  for k in range(1,len(time_array)):
    X.append(X_int(t_rest[k]-t_rest[k-1], X[k-1], tau_rest, avg_mag_rest, sigma_rest, Nsteps))
  X=array(X)
  return X 
  
def delayed_lightcurve(time_array, delay, delta_mag, redshift, tau, avg_mag, sigma, Nsteps):
    t1=time_array.copy()
    t2=time_array.copy() - delay
    id1=zeros(len(t1))
    id2=ones(len(t1))
    
    t_cat=concatenate([t1,t2])
    id_cat=concatenate([id1,id2])
    
    #time order the array
    t_cat_sort  =  array([x for (x,y) in sorted(zip(t_cat,id_cat))])
    id_cat_sort =  array([y for (x,y) in sorted(zip(t_cat,id_cat))])
    
    X_0 = random.gauss( avg_mag, sigma*sqrt(tau/2.) )
    lc=drw_lightcurve(t_cat_sort, X_0, tau, sigma, avg_mag, redshift, Nsteps)
    
    lc1 = array([y for (x,y) in zip(id_cat_sort,lc) if x==0.])
    lc2 = array([y for (x,y) in zip(id_cat_sort,lc) if x==1.])
    
    lc2+=delta_mag
    return lc1,lc2

def gauss(x, *p):
  A, mu, sigma = p
  return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def remove_outliers(_t,_m1, _m2, _em1, _em2, display = False):
  t_d = (_t[1:])[:-1]
  d1 = (2.*(_m1[1:])[:-1]-_m1[:-2]-_m1[2:])/sqrt((4.*_em1[1:])[:-1]**2+_em1[:-2]**2+_em1[2:]**2)
  d2 = (2.*(_m2[1:])[:-1]-_m2[:-2]-_m2[2:])/sqrt((4.*_em2[1:])[:-1]**2+_em2[:-2]**2+_em2[2:]**2)
  #t_d = _t[:-1]
  #d1 = (_m1[:-1]-_m1[1:])/sqrt(_em1[1:]**2+_em1[:-1]**2)
  #d2 = (_m2[:-1]-_m2[1:])/sqrt(_em2[1:]**2+_em2[:-1]**2)
  #d1 = (diff(_m1))/sqrt(_em1[1:]**2+_em1[:-1]**2)
  #d2 = (diff(_m2))/sqrt(_em2[1:]**2+_em2[:-1]**2)
  #d1 = (diff(_m1)/diff(_t))/sqrt(_em1[1:]**2+_em1[:-1]**2)
  #d2 = (diff(_m2)/diff(_t))/sqrt(_em2[1:]**2+_em2[:-1]**2)
  hist1, bin_edges1 = np.histogram(d1, bins=int(len(_t)/10.), density=False)
  bin_centres1 = (bin_edges1[:-1] + bin_edges1[1:])/2
  p01 = [10., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
  coeff1, var_matrix1 = curve_fit(gauss, bin_centres1, hist1, p0=p01)
  hist_fit1 = gauss(bin_centres1, *coeff1)
  d_th1 = abs(coeff1[2])
  hist2, bin_edges2 = np.histogram(d2, bins=int(len(_t)/10.), density=False)
  bin_centres2 = (bin_edges2[:-1] + bin_edges2[1:])/2
  p02 = [10., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
  coeff2, var_matrix2 = curve_fit(gauss, bin_centres2, hist2, p0=p02)
  # Get the fitted curve
  hist_fit2 = gauss(bin_centres2, *coeff2)
  d_th2 = abs(coeff2[2])

  d_thresh1 = 3.5*d_th1
  d_thresh2 = 3.5*d_th2

  if(display==True):
    figure(figsize=(10,12))
    ax=subplot(321)
    errorbar(_t,_m1,yerr=_em1, fmt='.')
    subplot(322,sharex=ax)
    errorbar(_t,_m2,yerr=_em2, fmt='.')
    subplot(323, sharex=ax)
    plot(t_d,d1, 'b.')
    subplot(324, sharex=ax)
    plot(t_d,d2, 'b.')
    subplot(325)
  # Get the fitted curve
    plot(bin_centres1, hist1, lw=3, drawstyle='steps-mid', label='Test data')
    plot(bin_centres1, gauss(bin_centres1,*coeff1), 'r-', lw=2)
    print len(_t)
    print d_th1
    #hist(d1, bins=100)
    subplot(326)
    plot(bin_centres2, hist2, lw=3, drawstyle='steps-mid', label='Test data')
    plot(bin_centres2, gauss(bin_centres2,*coeff2), 'r-', lw=2)
    print len(_t)
    print d_th2
    #hist(d1, bins=100)      #hist(d2, bins=100)

    ax=subplot(321)
    plot(t_d[abs(d1)>d_thresh1],((_m1[1:])[:-1])[abs(d1)>d_thresh1], '.', color='orange')
    subplot(322,sharex=ax)
    plot(t_d[abs(d2)>d_thresh2],((_m2[1:])[:-1])[abs(d2)>d_thresh2], '.', color='orange')
    subplot(323, sharex=ax)
    plot(t_d[abs(d1)>d_thresh1],d1[abs(d1)>d_thresh1], '.', color='orange')
    subplot(324, sharex=ax)
    plot(t_d[abs(d2)>d_thresh2],d2[abs(d2)>d_thresh2], '.', color='orange')

  # Outlier identification loop
  outlier_indices=[]
  for k in range(1,len(_t)-1):
    d1_val = (2.*(_m1[k])-_m1[k-1]-_m1[k+1])/sqrt((4.*_em1[k])**2+_em1[k-1]**2+_em1[k+1]**2)
    d2_val = (2.*(_m2[k])-_m2[k-1]-_m2[k+1])/sqrt((4.*_em2[k])**2+_em2[k-1]**2+_em2[k+1]**2)
    d_index = k-1
    #print d1_val, d1[d_index]
    if(abs(d1_val)>d_thresh1 or abs(d2_val)>d_thresh2):
      if(_t[k]-_t[k-1] < 90 and _t[k+1]-_t[k] < 90):
	  outlier_indices.append(k)
      #outlier_indices.append[k]
  print outlier_indices
  for k in outlier_indices:
    d_index = k-1
    if(k-1 in outlier_indices):
      print k
      if(abs(d1[d_index]) < abs(d1[d_index-1])):
	outlier_indices.pop(outlier_indices.index(k))
	
  for k in outlier_indices:
    d_index = k-1
    if(k-1 in outlier_indices):
      print k
      if(abs(d1[d_index]) > abs(d1[d_index-1])):
	outlier_indices.pop(outlier_indices.index(k-1))
	
  print outlier_indices
  for k in outlier_indices:
    d_index = k-1
    if(k-1 in outlier_indices):
      print k
      if(abs(d2[d_index]) < abs(d2[d_index-1])):
	outlier_indices.pop(outlier_indices.index(k))
	
  for k in outlier_indices:
    d_index = k-1
    if(k-1 in outlier_indices):
      print k
      if(abs(d2[d_index]) > abs(d2[d_index-1])):
	outlier_indices.pop(outlier_indices.index(k-1))	    
      
  print outlier_indices
  ax=subplot(321)
  plot(_t[outlier_indices],_m1[outlier_indices], 'rx')
  subplot(322,sharex=ax)
  plot(_t[outlier_indices],_m2[outlier_indices], 'rx')
  #subplot(323, sharex=ax)
  #plot(t_d[abs(d1)>d_thresh1],d1[abs(d1)>d_thresh1], 'r.')
  #subplot(324, sharex=ax)
  #plot(t_d[abs(d2)>d_thresh2],d2[abs(d2)>d_thresh2], 'r.')
  if(len(outlier_indices)!=0):
    _t = np.delete(_t, outlier_indices)
    _m1 = np.delete(_m1, outlier_indices)
    _m2 = np.delete(_m2, outlier_indices)
    _em1 = np.delete(_em1, outlier_indices)
    _em2 = np.delete(_em2, outlier_indices)

  if(display==True):
    subplot(321)
    errorbar(_t,_m1,yerr=_em1, fmt='g+')
    subplot(322,sharex=ax)
    errorbar(_t,_m2,yerr=_em2, fmt='g+')
    show()
  return _t, _m1, _m2, _em1, _em2

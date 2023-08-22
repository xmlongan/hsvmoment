'''
:abbr:`SV(Stochastic Volatility)` model with jumps in the return process

Extension to the Heston SV model, adds a jump component in the log price 
process.
The jump process is a compound Poisson process with rate :math:`\lambda` and
independent of everything else, the jump is a random variable distributed
according to a normal :math:`\mathcal{N}(\mu_J,\sigma_J^2)`.
'''
import math
import mdl_moments
import mdl_cmoments

def moment_J(rate, mu, h=1):
  moment = rate * h * mu
  return(moment)

def var_J(rate, mu, sigma, h=1):
  var = rate * h * (sigma**2 + mu**2)
  return(var)

def cm3_J(rate, mu, sigma, h=1):
  cm3 = rate * h * (sigma**3 + 3*mu*sigma**2)
  return(cm3)

def cm4_J(rate, mu, sigma, h=1):
  cm4  = rate * h * (mu**4 + 6*mu**2*sigma**2 + 3*sigma**4)
  cm4 += 3*(rate*h)**2 * (mu**2 + sigma**2)**2
  return(cm4)

# Model Moments
def moment_y(mu, theta, rate, mu_j, h=1):
  m_y = (mu - theta/2)*h + moment_J(rate, mu_j)
  return(m_y)

def var_y(k, theta, sigma_v, rho, rate, mu_j, sigma_j, h=1):
  h_tld = (-math.expm1(-k*h))/k
  var_y_old = theta*h + (sigma_v**2/(4*k**2) - rho*sigma_v/k)*theta*(h-h_tld)
  var = var_y_old + var_J(rate, mu_j, sigma_j, h)
  return(var)

def cov_yy(k, theta, sigma_v, rho, h=1):
  h_tld = (-math.expm1(-k*h))/k
  cov = theta * h_tld**2 *(sigma_v**2/(8*k) - rho*sigma_v/2)
  return(cov)

def cov_yy_lag2(k, theta, sigma_v, rho, h=1):
  cov = math.exp(-k*h) * cov_yy(k, theta, sigma_v, rho, h)
  return(cov)

def cov_y2y(mu, k, theta, sigma_v, rho, rate, mu_j, sigma_j, h=1):
  h_tld = (-math.expm1(-k*h))/k
  h_f = h * math.exp(-k*h) - h_tld
  #
  t1  = theta*sigma_v**2*mu*h/(4*k) 
  t1 -= (theta*sigma_v)**2*h/(8*k)
  t1 -= theta*sigma_v**2/(4*k)
  #
  t2 = (3*sigma_v**2/(2*k**2) - 2*rho*sigma_v/k) * theta * h_f
  t3 = (2*mu*theta - theta**2) * h * h_tld
  #
  cov_y2y_old  = h_tld*h_f * theta * sigma_v**4/(8*k**3)
  cov_y2y_old += t1 * h_tld**2
  cov_y2y_old -= rho * sigma_v * (h_tld/2) * (t1 + t2)
  #
  cov_y2y = cov_y2y_old + 2*moment_J(rate, mu_j)*cov_yy(k, theta, sigma_v, rho)
  return(cov_y2y)

def cov_yy2(mu, k, theta, sigma_v, rho, rate, mu_j, sigma_j, h=1):
  cov  = mdl_moments.cov_yy(1,2)
  # try
  cov = cov.replace('^', '**')
  cov = eval(cov, {'mu': mu, 'k' : k, 'theta': theta, 'sigma_v': sigma_v, 
                   'rho': rho, 'n': 2, 'h': 1})
  cov += 2*(rate*h)*mu_j * cov_yy(k, theta, sigma_v, rho, h)
  return(cov)

def cm3(mu, k, theta, sigma_v, rho, rate, mu_j, sigma_j, h=1):
  cm = mdl_cmoments.moment_y_central(3)
  # try
  cm = cm.replace('^', '**')
  cm = eval(cm,{'mu': mu, 'k' : k, 'theta': theta, 'sigma_v': sigma_v, 
                'rho': rho, 'n': 2, 'h': 1})
  cm += cm3_J(rate, mu_j, sigma_j, h)
  return(cm)

def cm4(mu, k, theta, sigma_v, rho, rate, mu_j, sigma_j, h=1):
  cm = mdl_cmoments.moment_y_central(4)
  # try
  cm = cm.replace('^', '**')
  cm = eval(cm,{'mu': mu, 'k' : k, 'theta': theta, 'sigma_v': sigma_v, 
                'rho': rho, 'n': 2, 'h': 1})
  h_tld = (-math.expm1(-k*h))/k
  var_y = theta*h + (sigma_v**2/(4*k**2) - rho*sigma_v/k) * theta * (h-h_tld)
  cm += 6 * var_y * var_J(rate, mu_j, sigma_j, h)
  cm += cm4_J(rate, mu_j, sigma_j, h)
  return(cm)



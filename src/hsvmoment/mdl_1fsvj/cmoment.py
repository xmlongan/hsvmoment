'''
Central Moments for One-Factor SV with jump
'''
import math

import sys, os
file_path = os.path.abspath(__file__) 
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(file_path)))
if src_dir not in sys.path: sys.path.append(src_dir)
# so that hsvmoment dir will be treated as a package

from hsvmoment.poly import Poly
from hsvmoment.mdl_1fsv.cmoment import cmoment_y as cm_y
# ['e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)']
from hsvmoment.cpp_mmnt import cmcpp
# ['lambda*h','mu_j','sigma_j^2']

def cmoment_y(l):
  '''Central moment of :math:`y_n` of order :math:`l`
  
  :param l: order of the moment.
  
  :return: poly with attribute ``keyfor`` =  
     ('e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)','lambda',
     'mu_j','sigma_j').
  :rtype: Poly
  '''
  poly = Poly()
  kf = ['e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)','lambda',
    'mu_j','sigma_j']
  poly.set_keyfor(kf)
  # if l == 0:
  #   poln = Poly({(0,0,0,0,0,0,0,0,0,0): 1}); poln.set_keyfor(kf); return(poln)
  # 
  for i in range(l+1):
    coef = math.comb(l, i)
    pol1 = cm_y(i)
    # ('e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)')
    pol2 = cmcpp(l-i)
    # ('lambda*h','mu','sigma')
    keyIndexes = [(0,1,2,3,4,5,6,-1,-1,-1),(-1,0,-1,-1,-1,-1,-1,0,1,2)]
    poln = pol1.mul_poly(pol2, keyIndexes, kf)
    poly.merge(coef * poln)
  return(poly)

##########
# scalar and (partial) derivative
##########

def dpoly(poly, wrt):
  '''Partial derivative of central moment w.r.t. parameter wrt
  
  :param poly: poly with attribute ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)','lambda',
     'mu_j','sigma_j')
  :param wrt: with respect to.
  
  :return: poly with attribute ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)','lambda',
     'mu_j','sigma_j').
  :rtype: Poly
  '''
  pold = Poly()
  kf = ('e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)','lambda',
    'mu_j','sigma_j')
  pold.set_keyfor(kf)
  # 
  # partial derivative w.r.t. mu
  if wrt == 'mu':
    return(pold)
  # partial derivative w.r.t. k
  elif wrt == 'k':
    for k in poly:
      if k[0] != 0:
        knw = (k[0],k[1]+1,k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[9])
        val = (-k[0]) * poly[k]
        pold.add_keyval(knw, val)
      if k[2] != 0:
        knw = (k[0],k[1],k[2]+1,k[3],k[4],k[5],k[6],k[7],k[8],k[9])
        val = (-k[2]) * poly[k]
        pold.add_keyval(knw, val)
  # partial derivative w.r.t. theta
  elif wrt == 'theta':
    for k in poly:
      if k[3] != 0:
        knw = (k[0],k[1],k[2],k[3]-1,k[4],k[5],k[6],k[7],k[8],k[9])
        val = k[3] * poly[k]
        pold.add_keyval(knw, val)
  # partial derivative w.r.t. sigma_v
  elif wrt == 'sigma_v':
    for k in poly:
      if k[4] != 0:
        knw = (k[0],k[1],k[2],k[3],k[4]-1,k[5],k[6],k[7],k[8],k[9])
        val = k[4] * poly[k]
        pold.add_keyval(knw, val)
  # partial derivative w.r.t. rho
  elif wrt == 'rho':
    for k in poly:
      if k[5] != 0:
        knw = (k[0],k[1],k[2],k[3],k[4],k[5]-1,k[6],k[7],k[8],k[9])
        val = k[5] * poly[k]
        pold.add_keyval(knw, val)
      if k[6] != 0:
        knw = (k[0],k[1],k[2],k[3],k[4],k[5]+1,k[6]-2,k[7],k[8],k[9])
        val = (-k[6]) * poly[k]
        pold.add_keyval(knw, val)
  elif wrt == 'lambda':
    for k in poly:
      if k[7] != 0:
        knw = (k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]-1,k[8],k[9])
        val = k[7] * poly[k]
        pold.add_keyval(knw, val)
  elif wrt == 'mu_j':
    for k in poly:
      if k[8] != 0:
        knw = (k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8]-1,k[9])
        val = k[8] * poly[k]
        pold.add_keyval(knw, val)
  elif wrt == 'sigma_j':
    for k in poly:
      if k[9] != 0:
        knw = (k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[9]-1)
        val = k[9] * poly[k]
        pold.add_keyval(knw, val)
  else:
    candidates = "'k','theta','sigma_v','rho','lambda','mu_j','sigma_j'"
    raise ValueError(f"wrt must be one of {candidates}!")
  return(pold)

def poly2num(poly, par):
  '''Decode poly back to scalar
  
  :param poly: poly to be decoded with attribute ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)','lambda',
     'mu_j','sigma_j')
  :param par: parameters in dict.
  
  :return: scalar of the poly.
  :rtype: float
  '''
  k = par['k']
  h = par['h']
  theta = par['theta']
  sigma_v = par['sigma_v']
  rho = par['rho']
  lmbd = par['lambda']
  mu_j = par['mu_j']
  sigma_j = par['sigma_j']
  # 
  value = 0
  for K in poly:
    val = poly[K] * math.exp(-K[0]*k*h) * (h ** K[1]) * (k ** (-K[2]))
    val *= (theta ** K[3]) * (sigma_v ** K[4]) * (rho ** K[5])
    val *= (1-rho**2) ** (K[6]/2)
    val *= (lmbd ** K[7]) * (mu_j ** K[8]) * (sigma_j ** K[9])
    value += val
  return(value)

def cm(l, par):
  '''Central moment in scalar
  
  :param l: order of the central moment.
  :param par: parameters in dict.
  
  :return: scalar of the central moment.
  :rtype: float
  '''
  cmoment = cmoment_y(l)
  value = poly2num(cmoment, par)
  return(value)

def dcm(l, par, wrt):
  '''Partial derivative of central moment w.r.t. parameter wrt
  
  :param l: order of the central moment.
  :param par: parameters in dict.
  :param wrt: with respect to.
  
  :return: scalar of the partial derivative.
  :rtype: float
  '''
  cmoment = cmoment_y(l)
  pold = dpoly(cmoment, wrt)
  value = poly2num(pold, par)
  return(value)


if __name__ == "__main__":
  # test the module
  from pprint import pprint
  # 
  keyfor = ['e^{-kh}','h','k^{-}','theta','sigma_v','rho','sqrt(1-rho^2)',
    'lambda','mu_j','sigma_j^2']
  print(f"cmoment_y() returns poly with keyfor = {keyfor}")
  print("cmoment_y(l=0): "); pprint(cmoment_y(l=0)) # verified
  print("cmoment_y(l=1): "); pprint(cmoment_y(l=1)) # verified
  # 
  print("cmoment_y(l=2): "); pprint(cmoment_y(l=2)) # verified (simply)
  # print("cmcpp(n=2): "); pprint(cmcpp(n=2))
  # print("cm_y(l=2): "); pprint(cm_y(l=2))
  print("cmoment_y(l=3): "); pprint(cmoment_y(l=3))

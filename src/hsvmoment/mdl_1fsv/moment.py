'''
Moments for One-Factor SV
'''
import math
from fractions import Fraction as Frac

import sys, os
file_path = os.path.abspath(__file__) 
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(file_path)))
if src_dir not in sys.path: sys.path.append(src_dir)
# so that hsvmoment dir will be treated as a package

from hsvmoment.poly import Poly
from hsvmoment.ito_mmnt import moment_v, moment_eIII
from hsvmoment.mdl_1fsv.cmoment import c_n, ve_eIII

def b_n(n1, n2, n3, n4, n5):
  '''Represent :math:`b_2(\\boldsymbol{n})` in :eq:`b2-n` in Poly
  
  :param n1: times of constant being selected.
  :param n2: times of :math:`v_{n-1}` being selected.
  :param n3: times of :math:`eI_n` being selected.
  :param n4: times of :math:`I_n` being selected.
  :param n5: times of :math:`I_n^{*}` being selected.
  
  :return: a poly with attibute ``keyfor`` =  
     ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)').
  :rtype: Poly
  '''
  b = Poly()
  kf = ['e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)']
  b.set_keyfor(kf)
  # 
  for i1 in range(n1, -1, -1):           # i1: mu * h
    for i2 in range(n1-i1, -1, -1):      # i2: -theta/2 * h
      for i3 in range(n1-i1-i2, -1, -1): # i3: 1/(2k) * theta
        i4 = n1 - i1 - i2 - i3           # i4: -e^{-kh}/(2k) * theta
        for i in range(n2+1):            # i : -e^{-kh}/(2k)
          for j in range(n4+1):          # j : -sigma_v/(2k)
            num  = math.comb(n1, i1) 
            num *= math.comb(n1-i1, i2)
            num *= math.comb(n1-i1-i2, i3)
            num *= math.comb(n2, i)
            num *= math.comb(n4, j)
            #   e^{-kh},h    ,k^{-}        ,mu,theta,sigma_v,rho,sqrt(1-rho^2)
            key = (i4+i,i1+i2,i3+i4+n2+n3+j,i1,n1-i1,n3+j,  n4-j,n5)
            num *= (-1)**(i2+i4 + n2 + i + j)
            den  = 2**(n1-i1 + n2 + n3 + j)
            b.add_keyval(key, Frac(num, den))
  return(b)

def moment_comb(n, n1, n2, n3, n4, n5):
  '''Moment for this combination in expansion of :math:`y_n^l`
  
  :param n: l in :math:`E[\overline{y}_{n}^l]`.
  :param n1: times of constant being selected.
  :param n2: times of :math:`v_{n-1}` being selected.
  :param n3: times of :math:`eI_n` being selected.
  :param n4: times of :math:`I_n` being selected.
  :param n5: times of :math:`I_n^{*}` being selected.
  
  :return: poly with attribute ``keyfor`` = 
     ('e^{-kh}','h','v_{n-1}','k^{-}','mu','theta','sigma_v','rho',
     'sqrt(1-rho^2)').
  :rtype: Poly
  '''
  poly = ve_eIII(n2, n3, n4, n5)
  # ['e^{-kh}','h','v_{n-1}','k^{-}','theta','sigma_v']
  b = b_n(n1, n2, n3, n4, n5)
  # ['e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)']
  # 
  keyfor = ['e^{-kh}','h','v_{n-1}','k^{-}','mu','theta','sigma_v','rho',
    'sqrt(1-rho^2)']
  # 
  keyIndexes = ([0,1,2,3,-1,4,5,-1,-1],[0,1,-1,2,3,4,5,6,7])
  # -1 for not having this term
  poly = poly.mul_poly(b, keyIndexes, keyfor)
  #
  c = c_n(n, n1, n2, n3, n4, n5)
  return(c * poly)

def sub_v(poly):
  '''Substitute :math:`v_{n-1}` with its moment (poly) in Moment of y_n
  
  :param poly: poly with attibute ``keyfor`` = 
     ('e^{-kh}','h','v_{n-1}','k^{-}','mu','theta','sigma_v','rho',
     'sqrt(1-rho^2)').
  
  :return: poly with attibute ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)').
  :rtype: Poly
  '''
  poly_sum = Poly()
  kf = ['e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)']
  poly_sum.set_keyfor(kf)
  # 
  for k in poly:
    poln = moment_v(k[2]) # ['theta','sigma_v^2/k']
    for key in poln:
      i, j = key
      knw = (k[0], k[1], k[3]+j, k[4], k[5]+i, k[6]+2*j, k[7], k[8])
      val = poln[key] * poly[k]
      poly_sum.add_keyval(knw, val)
  return(poly_sum)

def moment_y(l):
  '''Moment of :math:`y_n` of order :math:`l`
  
  :param l: order of the moment.
  
  :return: poly with attibute ``keyfor`` =  
     ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)').
  :rtype: Poly
  '''
  poly = Poly()
  kf = ['e^{-kh}','h','v_{n-1}','k^{-}','mu','theta','sigma_v','rho',
    'sqrt(1-rho^2)']
  poly.set_keyfor(kf)
  # 
  for n1 in range(l, -1, -1):
    for n2 in range(l-n1, -1, -1):
      for n3 in range(l-n1-n2, -1, -1):
        for n4 in range(l-n1-n2-n3, -1, -1):
          n5 = l - n1 - n2 - n3 - n4
          poly.merge(moment_comb(l, n1, n2, n3, n4, n5))
  poly_sum = sub_v(poly)
  poly_sum.remove_zero()
  return(poly_sum)

##########
# scalar and (partial) derivative
##########

def dpoly(poly, wrt):
  '''Partial derivative of central moment w.r.t. parameter wrt
  
  :param poly: central moment represented by the poly whoes attribute
     ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)').
  :param wrt: with respect to.
  
  :return: poly with attribute ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)').
  :rtype: Poly
  '''
  pold = Poly()
  kf = ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)')
  pold.set_keyfor(kf)
  # 
  # partial derivative w.r.t. k
  if wrt == 'k':
    for k in poly:
      if k[0] != 0:
        knw = list(k); knw[1] += 1
        val = (-k[0]) * poly[k]
        pold.add_keyval(tuple(knw), val)
      if k[2] != 0:
        knw = list(k); knw[2] += 1
        val = (-k[2]) * poly[k]
        pold.add_keyval(tuple(knw), val)
  elif wrt in ['mu','theta','sigma_v']:
    if wrt == 'mu': i = 3
    if wrt == 'theta': i = 4
    if wrt == 'sigma_v': i = 5
    for k in poly:
      if k[i] != 0:
        knw = list(k); knw[i] -= 1
        val = k[i] * poly[k]
        pold.add_keyval(tuple(knw), val)
  # partial derivative w.r.t. rho
  elif wrt == 'rho':
    for k in poly:
      if k[6] != 0:
        knw = list(k); knw[6] -= 1
        val = k[6] * poly[k]
        pold.add_keyval(tuple(knw), val)
      if k[7] != 0:
        knw = list(k); knw[6] += 1; knw[7] -= 2
        val = (-k[7]) * poly[k]
        pold.add_keyval(tuple(knw), val)
  else:
    msg = "wrt must be either 'k','mu','theta','sigma_v' or 'rho'!"
    raise ValueError(msg)
  return(pold)

def poly2num(poly, par):
  '''Decode poly back to scalar
  
  :param poly: poly to be decoded with attribute ``keyfor`` = 
     ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)').
  :param par: parameters in dict.
  
  :return: scalar of the poly.
  :rtype: float
  '''
  # ('e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)')
  k = par['k']
  h = par['h']
  mu = par['mu']
  theta = par['theta']
  sigma_v = par['sigma_v']
  rho = par['rho']
  # 
  value = 0
  for K in poly:
    val = poly[K] * math.exp(-K[0]*k*h) * (h ** K[1]) * (k ** (-K[2]))
    val *= (mu ** K[3])
    val *= (theta ** K[4]) * (sigma_v ** K[5]) * (rho ** K[6])
    val *= (1-rho**2) ** (K[7]/2)
    value += val
  return(value)

def m(l, par):
  '''Moment in scalar
  
  :param l: order of the central moment.
  :param par: parameters in dict.
  
  :return: scalar of the central moment.
  :rtype: float
  '''
  moment = moment_y(l)
  value = poly2num(moment, par)
  return(value)

def dm(l, par, wrt):
  '''Partial derivative of moment w.r.t. parameter wrt
  
  :param l: order of the moment.
  :param par: parameters in dict.
  :param wrt: with respect to.
  
  :return: scalar of the partial derivative.
  :rtype: float
  '''
  moment = moment_y(l)
  pold = dpoly(moment, wrt)
  value = poly2num(pold, par)
  return(value)


if __name__ == "__main__":
  # test the module
  from pprint import pprint
  # 
  # print("b_n(n1=1,n2=1,n3=1,n4=1,n5=1) = ")
  # pprint(b_n(n1=1,n2=1,n3=1,n4=1,n5=1))
  # 
  keyfor = ['e^{-kh}','h','k^{-}','mu','theta','sigma_v','rho','sqrt(1-rho^2)']
  print(f"moment_y() returns poly with keyfor = {keyfor}")
  print("moment_y(l=0): "); pprint(moment_y(l=0))
  print("moment_y(l=1): "); pprint(moment_y(l=1)) # verified
  print("moment_y(l=2): "); pprint(moment_y(l=2))
  print("moment_y(l=3): "); pprint(moment_y(l=3))

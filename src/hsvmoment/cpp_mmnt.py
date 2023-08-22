'''
Compound Poisson Process variable moments

Compound Poisson Process

.. math::
   
   J(t) = \sum_{i=0}^{N(t)}j_i,\quad j_i \sim \mathcal{N}(\mu_j,\sigma_j^2)

where :math:`N(t)` is a Poisson process with rate :math:`\lambda`.
For our interest, I define variable 
:math:`J_n = \sum_{i=N((n-1)h)+1}^{N(nh)}j_i`.

Moment-Generating Function
===========================

For variable :math:`J_n`, its moment-generating function 

.. math::
   
   M_{J_n}(s) = e^{\lambda h (M_{j_i}(s)-1)},
   \qquad
   M_{j_i}(s) = e^{\mu_js+\\frac{1}{2}\sigma_j^2s^2}.

where :math:`M_{j_i}(s)` is the moment-generating function of normal variable
:math:`j_i`.

MGF - :abbr:`CPP(Compound Poisson Process)`
--------------------------------------------

For the first three derivatives,

.. math::
   
   M_{J_n}^{'}(s) 
   &= e^{\lambda h (M_{j_i}(s)-1)} (\lambda h) M_{j_i}^{'}(s),\\\\
   M_{J_n}^{''}(s)
   &= e^{\lambda h (M_{j_i}(s)-1)} \left[ (\lambda h)^2 M_{j_i}^{'2}(s)
   + (\lambda h) M_{j_i}^{''}(s) \\right],\\\\
   M_{J_n}^{(3)}(s)
   &= e^{\lambda h (M_{j_i}(s)-1)} \left[ (\lambda h)^3 M_{j_i}^{'3}(s)
   + 3(\lambda h)^2 M_{j_i}^{'}(s) M_{j_i}^{''}(s) 
   + (\lambda h) M_{j_i}^{(3)}(s) \\right].

I propose to represent derivative of any order as

.. math::
   
   M^{(n)}_{J_n}(s)
   = \sum_{(i,(n1,m1),...,(nl,ml))} b_{key} (\lambda h)^i M_{j_i}^{(n1)m1}(s)
   \cdots M_{j_i}^{(nl)ml}(s)

where :math:`n1 < \cdots < nl` and :math:`key=(i,(n1,m1),...,(nl,ml))`.
Then its derivative goes as

.. math::

   M_{J_n}^{(n+1)}(s)
   &= \sum_{(i,(n1,m1),...,(nl,ml))} b_{key} (\lambda h)^{i+1} M_{j_i}^{(1)}(s)
   M_{j_i}^{(n1)m1}(s)\cdots M_{j_i}^{(nl)ml}(s)\\\\
   &\quad+ \sum_{(i,(n1,m1),...,(nl,ml))} b_{key} (\lambda h)^i 
   (M_{j_i}^{(n1)m1}(s)\cdots M_{j_i}^{(nl)ml}(s))'.

Rearrage the derivative to represent it as that of :math:`M^{(n)}_{J_n}(s)`.


MGF - Normal Distribution
--------------------------

For the first three derivatives,

.. math::
   
   M_{j_i}^{'}(s)
   &= e^{\mu_js+\\frac{1}{2}\sigma_j^2s^2} (\mu_j + \sigma_j^2s),\\\\
   M_{j_i}^{''}(s)
   &= e^{\mu_js+\\frac{1}{2}\sigma_j^2s^2} \left[ (\mu_j + \sigma_j^2s)^2
   +  \sigma_j^2 \\right],\\\\
   M_{j_i}^{(3)}(s)
   &= e^{\mu_js+\\frac{1}{2}\sigma_j^2s^2} \left[ (\mu_j + \sigma_j^2s)^3
   +  (\mu_j + \sigma_j^2s)\sigma_j^2 + 
   2(\mu_j + \sigma_j^2s)\sigma_j^2\\right].

I propose to represent derivative of any order as

.. math::
   
   M_{j_i}^{(n)}(s)
   = \sum_{i,j}b_{ij}(\mu+\sigma^2s)^i \sigma^{2j}.
   
Then its derivative is given as

.. math::

   M_{j_i}^{(n+1)}(s)
   = \sum_{i,j}b_{ij}(\mu+\sigma^2s)^{i+1} \sigma^{2j}
    +\sum_{i>0,j}b_{ij}i(\mu+\sigma^2s)^{i-1} \sigma^{2(j+1)}.

Rearrage the derivative to represent it as that of :math:`M^{(n)}_{j_i}(s)`.

In summary, I defined

1. :py:func:`~hsvmoment.cpp_mmnt.mnorm` to compute moment of normal variable
   which uses 
   :py:func:`~hsvmoment.cpp_mmnt.dmgf`.

2. :py:func:`~hsvmoment.cpp_mmnt.mcpp` to compute moment of CPP variable
   which uses
   :py:func:`~hsvmoment.cpp_mmnt.dmgf_cpp`, 
   :py:func:`~hsvmoment.cpp_mmnt.decode`.

3. :py:func:`~hsvmoment.cpp_mmnt.cmcpp` to compute central moment of CPP
   variable.
'''
from fractions import Fraction as Frac
import math

import sys, os
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if src_dir not in sys.path: sys.path.append(src_dir)

from hsvmoment.poly import Poly

#
# Normal distribution
#

def dmgf(poly):
  '''Derivative of normal Moment-Generating Function
  
  :param poly: poly with attribute ``keyfor`` = ('mu','sigma^2').
  
  :return: poly with attribute ``keyfor`` = ('mu','sigma^2').
  :rtype: Poly
  '''
  poln = Poly()
  poln.set_keyfor(['mu','sigma^2'])
  for k in poly:
    # step 1
    key = (k[0]+1, k[1])
    poln.add_keyval(key, poly[k])
    # step 2
    if k[0] > 0:
      knw = (k[0]-1, k[1]+1)
      poln.add_keyval(knw, poly[k]*k[0])
  return(poln)

def mnorm(n):
  '''Moment of Normal distribution
  
  :param n: order of the moment.
  
  :return: poly with attribute ``keyfor`` = ('mu','sigma^2').
  :rtype: Poly
  '''
  if n < 0:
    msg = f'mnorm(n): order of the moment is {n} which must be > 0!'
    raise ValueError(msg)
  M = []
  kf = ['mu','sigma^2']
  # n = 0
  poly = Poly({(0,0): 1}); poly.set_keyfor(kf)
  M.append(poly)
  if n == 0: return(poly)
  # n = 1
  poly = Poly({(1,0): 1}); poly.set_keyfor(kf)
  M.append(poly)
  if n == 1: return(poly)
  #
  for i in range(2, n+1):
    M.append(dmgf(M[-1]))
  return(M[-1])

#
# Compound Poisson Process
#

def merge(poly1, poly2):
  '''Merge two dict poly
    
    Insert new key-value or add value for the existing key.
    
    :param other: poly.
    
    :return: a new poly.
    :rtype: Poly
  '''
  for key in poly2.keys():
    if key not in poly1:
      poly1[key]  = poly2[key]
    else:
      poly1[key] += poly2[key]
  return(poly1)

def d1_times_key(key):
  '''Update the key after multiply with :math:`(\lambda h)M_{j_i}'(s)`
  
  :param key: a tuple of (i,(n1,m1),...,(nl,ml)) which
     stands for :math:`(\lambda h)^{i} M^{(n_1)m_1}(s) ... M^{(n_l)m_l}(s)`.
  
  :return: updated key.
  :rtype: list
  '''
  # multiply with (lambda h)
  knw = list(key)
  knw[0] += 1
  # multiply with M_{j_i}'(s)
  if key[1][0] == 1: # (n1,m1),(n2,m2),...,(nj,mj) where n1 < n2 < ... < nj
    knw[1] = (1, key[1][1] + 1)
  else:
    knw.insert(1, (1,1))
  return(tuple(knw))

def dterm(key, coef):
  '''Take derivative of each term
  
  :param key: key corresponding to a term, i.e.,
     (i,(n1,m1),(n2,m2),...,(nj,mj)) which
     stands for :math:`(\lambda h)^{i} M^{(n_1)m_1}(s) ... M^{(n_l)m_l}(s)`.
  :param coef: corresponding coefficient.
  
  :return: poly representing the derivative.
  :rtype: Poly
  '''
  poly = Poly()
  # key: tuple(i,(n1,m1),(n2,m2),...,(nj,mj)) where n1 < n2 < ... < nj
  for j in range(1, len(key)):
    knw = list(key)
    n_j, m_j = key[j]
    coef_new = coef * m_j # update coef
    knw[j] = (n_j, m_j-1) # update exponent M^{(n_j)m_j} -> M^{(n_j)(m_j-1)}
    # 
    if j < len(key) - 1:  # having next one
      n_nxt, m_nxt = key[j+1]
      if n_nxt == n_j + 1:
        knw[j+1] = (n_nxt, m_nxt+1)
      else:
        knw.insert(j+1, (n_j+1,1))
    else:                # no follower
      knw.append((n_j+1,1))
    #
    if knw[j][1] == 0: del knw[j] # delete M^{(n_j)0}
    poly.add_keyval(tuple(knw), coef_new)
  return(poly)

def dmgf_cpp(poly):
  '''derivative of Moment Generating Function of CPP
  
  :param poly: poly representing MGF with 
     key :math:`(i,(n_1,m_1),...,(n_l,m_l))` which
     stands for :math:`(\lambda h)^{i}M^{(n_1)m_1}(s) ... M^{(n_l)m_l}(s)`.
  
  :return: poly with the same key-val.
  :rtype: Poly
  '''
  # step 1
  poly_sum = Poly()
  for key in poly:
    knw = d1_times_key(key)
    poly_sum.add_keyval(knw, poly[key])
  # step 2
  for key in poly:
    poly_sum = merge(dterm(key, poly[key]), poly_sum)
  #
  return(poly_sum)

def poly_pow(poly, n):
  '''Raise poly to power n
  
  :param poly: poly with attribute ``keyfor`` = ('mu','sigma^2').
  :param n: power.
  
  :return: poly with attribute ``keyfor`` = ('mu','sigma^2').
  :rtype: Poly
  '''
  poln = Poly(poly)
  poln.set_keyfor(['mu','sigma^2'])
  for i in range(n-1):
    poly = poly * poln
  return(poly)

def decode(poly):
  '''Decode intermediate poly to target poly
  
  :param poly: poly representation of MGF with 
     key :math:`(i,(n_1,m_1),...,(n_l,m_l))` which
     stands for :math:`(\lambda h)^{i} M^{(n_1)m_1}(s) ... M^{(n_l)m_l}(s)`.
  
  :return: poly with attribute ``keyfor`` = ('lambda*h','mu','sigma^2').
  :rtype: Poly
  '''
  poln = Poly()
  poln.set_keyfor(['lambda*h','mu','sigma^2'])
  for k in poly:
    i = k[0]
    pol1 = Poly({(0,0): 1})
    pol1.set_keyfor(['mu','sigma^2'])
    for j in range(1, len(k)):
      n_j, m_j = k[j]
      polj = poly_pow(mnorm(n_j), m_j)
      pol1 = pol1 * polj
    # 
    for key in pol1:
      knw = (i, key[0], key[1])
      val = poly[k] * pol1[key]
      poln.add_keyval(knw, val)
  poln.remove_zero()
  return(poln)

def mcpp_original(n):
  '''Moment of Compound Poisson Process variable of order n
  
  :param n: order of the moment.
  
  :return: poly with attribute ``keyfor`` = ('lambda*h','mu','sigma^2').
  :rtype: Poly
  '''
  if n < 0:
    raise ValueError(f'mcpp(n): order of the moment is {n} which must > 0!')
  M = []
  kf = ['(lambda h)^{i} M^{(n_1)m_1}(s) ... M^{(n_l)m_l}(s)']
  # n = 0
  poly = Poly({(0,): 1}); poly.set_keyfor(kf)
  M.append(poly)
  #
  if n == 0:
    poln = Poly({(0,0,0): 1}); poln.set_keyfor(['lambda*h','mu','sigma^2'])
    return(poln)
  # n = 1
  poly = Poly({(1,(1,1)): 1}); poly.set_keyfor(kf)
  M.append(poly)
  #
  if n == 1:
    poln = Poly({(1,1,0): 1}); poln.set_keyfor(['lambda*h','mu','sigma^2'])
    return(poln)
  #
  for i in range(2, n+1):
    M.append(dmgf_cpp(M[-1]))
  return(decode(M[-1]))

def mcpp(n):
  '''Moment of Compound Poisson Process variable of order n
  
  :param n: order of the moment.
  
  :return: poly with attribute ``keyfor`` = ('lambda*h','mu','sigma').
  :rtype: Poly
  '''
  poly = mcpp_original(n) # ('lambda*h','mu','sigma^2')
  poln = Poly()
  poln.set_keyfor(['lambda*h','mu','sigma'])
  for k in poly:
    key = (k[0],k[1],2*k[2])
    val = poly[k]
    poln.add_keyval(key, val)
  return(poln)
  

def cmcpp(n):
  '''Central Moment of Compound Poisson Process of order n
  
  :param n: order of the moment.
  
  :return: poly with attribute ``keyfor`` = ('lambda*h','mu','sigma')
  :rtype: Poly
  '''
  if n < 0:
    raise ValueError(f'cmcpp(n): Order of the moment is {n} which must > 0!')
  poly = Poly()
  poly.set_keyfor(['lambda*h','mu','sigma'])
  for i in range(n+1):
    coef = math.comb(n, i) * ((-1)**(n-i))
    poln = mcpp(i) # ('lambda*h','mu','sigma')
    # * (lambda * h * mu)^{n-i}
    for k in poln:
      key = (k[0]+n-i, k[1]+n-i, k[2])
      val = poln[k] * coef
      poly.add_keyval(key, val)
  poly.remove_zero()
  return(poly)


if __name__ == '__main__':
  # test module cpp_mmnt
  from pprint import pprint
  # test moment calculation of normal distribution
  print("mnorm(n) in poly with keyfor = ('mu','sigma^2'): ")
  print('mnorm(0): '); pprint(mnorm(n=0))
  print('mnorm(1): '); pprint(mnorm(n=1))
  print('mnorm(2): '); pprint(mnorm(n=2))
  print('mnorm(3): '); pprint(mnorm(n=3))
  print('mnorm(8): '); pprint(mnorm(n=8))
  # test moment calculation of Compound Poisson Process variable
  print("mcpp(n) in poly with keyfor = ('lambda*h','mu','sigma'): ")
  print('mcpp(0): '); pprint(mcpp(n=0))
  print('mcpp(1): '); pprint(mcpp(n=1))
  print('mcpp(2): '); pprint(mcpp(n=2))
  print('mcpp(3): '); pprint(mcpp(n=3))
  print('mcpp(4): '); pprint(mcpp(n=4))
  # test central moment calculation
  print("cmcpp(n) in poly with keyfor = ('lambda*h','mu','sigma'): ")
  print('cmcpp(0): '); pprint(cmcpp(n=0))
  print('cmcpp(1): '); pprint(cmcpp(n=1))
  print('cmcpp(2): '); pprint(cmcpp(n=2))
  print('cmcpp(3): '); pprint(cmcpp(n=3))
  print('cmcpp(4): '); pprint(cmcpp(n=4))

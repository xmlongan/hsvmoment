'''
Itô process moment under One Square-Root Diffusion Processes

The content has also been explained in :doc:`../design` page.

Insights
=========

All :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]` can be represented as
a "Polynomial" of following form

.. _polynomial-representation:

.. math::
   
   &E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]\\\\
   &= \sum_{n_3,i,j,l,o,p,q} b_{n_3ijlopq} e^{n_3k(n-1)h} e^{ik[t-(n-1)h]}
   [t-(n-1)h]^jv_{n-1}^l k^{-o}\\theta^p\sigma_v^q

where :math:`b_{ijlopq}` is the coefficient. 
:math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]` can be 
represented similarly.

To facilitate the representation and corresponding operations, I designed
a new *class* :py:class:`~hsvmoment.poly.Poly` which is derived from 
:class:`~collections.UserDict` in the Python Standard Library 
`collections <https://docs.python.org/3/library/collections.html>`_.

Integrals
==========

The essential computation in recursive equations :eq:`ito-moment-i` 
and :eq:`ito-moment-ii` of :doc:`../theory` is that of following integral

.. math::
   
   \int_{(n-1)h}^t e^{ik[s-(n-1)h]} [s-(n-1)h]^j ds.


For the indefinite integral, we have

.. math::
   
   \int e^{nkt} t^m dt =
   \\begin{cases}
   \sum_{i=0}^m c_{nmi} \\frac{1}{k^{i+1}}e^{nkt} t^{m-i} 
    & \\text{if } n\\neq 0, m \\neq 0,\\\\
   \\frac{1}{nk}e^{nkt}t^0 & \\text{if } n\\neq 0, m = 0,\\\\
   \\frac{1}{m+1}e^{0kt}t^{m+1} & \\text{if } n = 0, m \\neq 0,\\\\
   e^{0kt}t^1 & \\text{if } n =0 , m=0,
   \end{cases}

where :math:`c_{nm0} \\triangleq \\frac{1}{n}` and

.. math::
   
   c_{nmi} \\triangleq \\frac{(-1)^{i}}{n^{i+1}} \prod_{j=m-i+1}^{m} j,
   \quad 1\le i \le m.

The coefficient :math:`c_{nmi}` is implemented in function
:py:func:`~hsvmoment.ito_mmnt.c_nmi` which returns :class:`fractions.Fraction`
instead of decimal (float number).

For the definite integral,

.. math::
   
   \int_{(n-1)h}^t e^{ik[s-(n-1)h]}[s-(n-1)h]^jds = F(t-(n-1)h) - F(0)
    
where :math:`F(t) = \int e^{nkt} t^m dt`. The definite integral is implemented 
in :py:func:`~hsvmoment.ito_mmnt.int_et`.


Polynomial Representation
--------------------------

The result of the integral, returned by :py:func:`~hsvmoment.ito_mmnt.int_et`,
is represented as a "polynomial" of following form

.. math::
   
   \int_{(n-1)h}^t e^{ik[s-(n-1)h]} [s-(n-1)h]^j ds
    = \sum_{i,j^{'},l}c_{ij^{'}l}e^{ik[t-(n-1)h]}[t-(n-1)h]^{j^{'}}k^{-l}

which is encoded in a :py:class:`~hsvmoment.poly.Poly`, derived from
:class:`collections.UserDict`, with
:code:`keyfor = ('e^{k[t-(n-1)h]}','[t-(n-1)h]','k^{-}')`,
``key`` = :math:`(i,j^{'},l)` and ``value`` = :math:`c_{ij^{'}l}`.


Code Design
============

Itô process moment - I
-----------------------

With :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]` represented as 
a "polynomial" of following form

.. math::
   
   &E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]\\\\
   &= \sum_{n_3,i,j,l,o,p,q} b_{n_3ijlopq} e^{n_3k(n-1)h} e^{ik[t-(n-1)h]}
   [t-(n-1)h]^jv_{n-1}^l k^{-o}\\theta^p\sigma_v^q,

consequently, we have

.. math::
   
   &e^{-kt}E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]\\\\
   &= \sum_{n_3,i,j,l,o,p,q} b_{n_3ijlopq} e^{(n_3-1)k(n-1)h}
   e^{(i-1)k[t-(n-1)h]}[t-(n-1)h]^jv_{n-1}^l k^{-o}\\theta^p\sigma_v^q,\\\\
   &e^{kt}E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]\\\\
   &= \sum_{n_3,i,j,l,o,p,q} b_{n_3ijlopq} e^{(n_3+1)k(n-1)h}
   e^{(i+1)k[t-(n-1)h]}[t-(n-1)h]^jv_{n-1}^l k^{-o}\\theta^p\sigma_v^q,\\\\
   &e^{2kt}E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]\\\\
   &= \sum_{n_3,i,j,l,o,p,q} b_{n_3ijlopq} e^{(n_3+2)k(n-1)h}
   e^{(i+2)k[t-(n-1)h]}[t-(n-1)h]^jv_{n-1}^l k^{-o}\\theta^p\sigma_v^q.

Therefore, it's profitable to consider following generic integral

.. math::
   
   &\int_{(n-1)h}^t e^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}|v_{n-1}]ds\\\\ 
   &= \sum_{n_3,i,j,l,o,p,q} b_{n_3ijlopq} e^{(n_3+m)k(n-1)h} 
   \cdot int\_et(i+m,j)\cdot v_{n-1}^l k^{-o}\\theta^p\sigma_v^q\\\\
   &= \sum_{n_3+m,i+m,j^{'},l,o^{'},p,q} b_{(n_3+m)(i+m)j^{'}l o^{'}pq} 
   e^{(n_3+m)k(n-1)h} e^{(i+m)k[t-(n-1)h]} [t-(n-1)h]^{j^{'}}\\\\
   &\qquad \cdot v_{n-1}^{l} k^{-o^{'}}\\theta^{p}\sigma_v^{q}

where 

.. math::
   
   int\_et(i+m,j) 
   =\sum_{i+m,j^{'},l^{'}} c_{(i+m)j^{'}l^{'}}e^{(i+m)k[t-(n-1)h]}
   [t-(n-1)h]^{j^{'}} k^{-l^{'}}.

Implementation:

1. Function :py:func:`~hsvmoment.ito_mmnt.int_meII` in module 
   :py:mod:`~hsvmoment.ito_mmnt` is defined to accomplish the computation in
   equation :eq:`int-meII`.

2. Function :py:func:`~hsvmoment.ito_mmnt.recursive_eII` in module
   :py:mod:`~hsvmoment.ito_mmnt` is defined to realize
   the recursive step in equation :eq:`ito-moment-i` of :doc:`../theory`.

3. Function :py:func:`~hsvmoment.ito_mmnt.moment_eII` in module
   :py:mod:`~hsvmoment.ito_mmnt` is implemented to calculate
   :math:`E[eI_n^{n_3}I_n^{n_4}|v_{n-1}]`.

For demonstration, I re-write the following initial three moments in
:ref:`ito-recursive-i` in :doc:`../theory` according to the "polynomial"
representation

.. math::
   
   E[eI_{n-1,t}^2|v_{n-1}]
   &=& \\frac{1}{2}&e^{2k(n-1)h} e^{2k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0\\\\
   && + &e^{2k(n-1)h}e^{k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^1
   k^{-1}\\theta^0\sigma_v^0\\\\
   && - &e^{2k(n-1)h}e^{k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0\\\\
   && - &e^{2k(n-1)h}e^{0k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^1
   k^{-1}\\theta^0\sigma_v^0\\\\
   && + \\frac{1}{2} &e^{2k(n-1)h}e^{0k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0,\\\\
   %
   E[eI_{n-1,t}I_{n-1,t}|v_{n-1}]
   &=& &e^{k(n-1)h} e^{k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0\\\\
   && +&e^{k(n-1)h} e^{0k[t-(n-1)h]}[t-(n-1)h]^1v_{n-1}^1
   k^{-0}\\theta^0\sigma_v^0\\\\
   && -&e^{k(n-1)h} e^{0k[t-(n-1)h]}[t-(n-1)h]^1v_{n-1}^0
   k^{-0}\\theta^1\sigma_v^0\\\\
   && -&e^{k(n-1)h} e^{0k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0,\\\\
   %
   E[I_{n-1,t}^2|v_{n-1}]
   &=&-&e^{0k(n-1)h} e^{-k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^1
   k^{-1}\\theta^0\sigma_v^0\\\\
   && +&e^{0k(n-1)h} e^{-k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0\\\\
   && +&e^{0k(n-1)h} e^{0k[t-(n-1)h]}[t-(n-1)h]^1v_{n-1}^0
   k^{-0}\\theta^1\sigma_v^0\\\\
   && +&e^{0k(n-1)h} e^{0k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^1
   k^{-1}\\theta^0\sigma_v^0\\\\
   && -&e^{0k(n-1)h} e^{0k[t-(n-1)h]}[t-(n-1)h]^0v_{n-1}^0
   k^{-1}\\theta^1\sigma_v^0.


Itô process moment - II
------------------------

Implementation:

1. Define :py:func:`~hsvmoment.ito_mmnt.int_meIII` similarly.

2. Define :py:func:`~hsvmoment.ito_mmnt.recursive_eIII` to realize the 
   recursive step in equation :eq:`ito-moment-ii` of :doc:`../theory`.

3. Define :py:func:`~hsvmoment.ito_mmnt.moment_eIII` to finish the computation 
   of  :math:`E[eI_n^{n_3}I_n^{n_4}I_n^{*n_5}|v_{n-1}]`.

'''
# from pprint import pprint
from fractions import Fraction as Frac

import sys, os
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if src_dir not in sys.path: sys.path.append(src_dir)

from hsvmoment.poly import Poly

def c_nmi(n,m,i):
  '''Coefficent :math:`c_{nmi}` as in :eq:`c-nmi`.
  
  :param n: n in :math:`e^{nkt}`.
  :param m: m in :math:`t^{m-i}`.
  :param i: i in :math:`t^{m-i}`.
  
  :return: the coefficient :math:`c_{nmi}`.
  :rtype: Fraction
  '''
  prod = 1
  for j in range(m-i+1, m+1):
    prod = prod * j
  den = n ** (i+1) # denumerator
  c = ((-1)**i) * Frac(prod, den)
  return(c)

def int_et(n,m):
  ''':math:`\int_{(n-1)h}^{t} e^{ik[s-(n-1)h]}[s-(n-1)h]^jds`
  
  :param n: i in :math:`e^{ik[s-(n-1)h]}[s-(n-1)h]^j`.
  :type integer: int
  :param m: j in :math:`e^{ik[s-(n-1)h]}[s-(n-1)h]^j`.
  :type integer: int
  
  :return: a poly with attribute ``keyfor`` = 
     ('e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'k^{-}').
  :rtype: Poly
  '''
  if m < 0:
    msg = f"m in int_et(n,m) equals {m}, however it shouldn't be negative!"
    raise ValueError(msg)
  # 
  poly = Poly()
  poly.set_keyfor(['e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'k^{-}'])
  # 
  if n == 0 and m == 0:
    poly[(0,1,0)] = Frac(1,1)
  elif n == 0 and m != 0:
    poly[(0,m+1,0)] = Frac(1,m+1)
  elif n != 0 and m == 0:
    poly[(n,0,1)] = Frac(1,n)
    # - F(0)
    poly[(0,0,1)] = Frac(-1,n)
  else:
    poly[(n,m,1)] = Frac(1,n)
    for i in range(1, m+1):
      c = c_nmi(n, m, i)
      poly[(n,m-i,i+1)] = c
      if i == m: # - F(0): - c_nmi
        poly[(0,0,i+1)] = -c # - c * 1/k^{i+1}
  return(poly)

def int_meII(m, n3, n4, eII):
  '''Integral of :math:`\int_{(n-1)h}^t e^{mks}eIIds` 
  
  :param m: m in :math:`\int_{(n-1)h}^t e^{mks}eIIds` where :math:`eII` is 
     :math:`E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}|v_{n-1}]`.
  :param n3: :math:`n_3` in the integral.
  :param n4: :math:`n_4` in the integral.
  :param eII: a dict with key (n3,n4) and value Poly object with attribute 
     ``keyfor`` = 
     ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  
  :return: poly with attribute ``keyfor`` = 
     ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  :rtype: Poly
  '''
  # poly of E[eI^{n3}I^{n4}]
  b = eII[(n3, n4)]
  # b: poly with ('e^{k(n-1)h}', 'e^{[t-(n-1)h]}', '[t-(n-1)h]','v_{n-1}',
  # 'k^{-}', 'theta', 'sigma_v')
  # 
  poly = Poly()
  kf = ['e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}', 'k^{-}', 
    'theta', 'sigma_v']
  poly.set_keyfor(kf)
  # 
  for k1 in b.keys():
    c = int_et(m+k1[1], k1[2]) # ('e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'k^{-}')
    for k2 in c.keys():
      key = (k1[0]+m, k2[0], k2[1], k1[3], k1[4]+k2[2], k1[5], k1[6])
      # k1[0]+m: compensate e^{mk[s-(n-1)h]} for e^{k(n-1)h}
      val = b[k1] * c[k2]
      poly.add_keyval(key, val)
  return(poly)

def int_meIII(m, n3, n4, n5, eIII):
  '''Integral of :math:`\int_{(n-1)h}^t e^{mks}eIII ds` 
  
  :param m: m in the integral where :math:`eIII` is
     :math:`E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}I_{n-1,s}^{*n_5}|v_{n-1}]`.
  :param n3: :math:`n_3` in the integral.
  :param n4: :math:`n_4` in the integral.
  :param n5: :math:`n_5` in the integral.
  :param eIII: a dict with key (n3,n4,n5) and value Poly object with 
     attribute ``keyfor`` = ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 
     'v_{n-1}','k^{-}', 'theta', 'sigma_v').
  
  :return: poly with attribute ``keyfor`` =
     ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  :rtype: Poly
  '''
  # poly of E[eI^{n3}I^{n4}I^{*n5}]
  b = eIII[(n3, n4, n5)]
  # b: poly with ('e^{k(n-1)h}', 'e^{[t-(n-1)h]}', '[t-(n-1)h]','v_{n-1}',
  # 'k^{-}', 'theta', 'sigma_v')
  # 
  poly = Poly()
  kf = ['e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}', 'k^{-}', 
    'theta', 'sigma_v']
  poly.set_keyfor(kf)
  #
  for k1 in b.keys():
    c = int_et(m+k1[1], k1[2]) # ('e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'k^{-}')
    for k2 in c.keys():
      key = (k1[0]+m, k2[0], k2[1], k1[3], k1[4]+k2[2], k1[5], k1[6])
      # k1[0]+m: compensate e^{mk[s-(n-1)h]} for e^{k(n-1)h}
      val = b[k1] * c[k2]
      poly.add_keyval(key, val)
  return(poly)

def coef_poly(coef, poly, tp):
  '''Multiply poly with different type coefficients
  
  :param coef: Fraction.
  :param poly: poly with attribute ``keyfor`` = 
     ('e^{k(n-1)h}', 'e^{[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  :param tp: type of the multiplication,
  
     +----+----------------------------+
     | tp | multiply with              |
     +====+============================+
     | 1  | :math:`e^{k(n-1)h}v_{n-1}` |
     +----+----------------------------+
     | 2  | :math:`-e^{k(n-1)h}\\theta` |
     +----+----------------------------+
     | 3  | :math:`\\theta`             |
     +----+----------------------------+
     | 4  | :math:`\\sigma_v`           |
     +----+----------------------------+
  
  :return: poly with attribute ``keyfor`` = 
     ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  :rtype: Poly
  '''
  poln = Poly()
  kf = ['e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 
        'v_{n-1}', 'k^{-}', 'theta', 'sigma_v']
  poln.set_keyfor(kf)
  # 
  if tp == 1:
    for k in poly:
      key = (k[0]+1, k[1], k[2], k[3]+1, k[4], k[5], k[6])
      val = coef * poly[k]
      poln.add_keyval(key, val)
  if tp == 2:
    for k in poly:
      key = (k[0]+1, k[1], k[2], k[3],   k[4], k[5]+1,k[6])
      val = (-coef) * poly[k]
      poln.add_keyval(key, val)
  if tp == 3:
    for k in poly:
      key = (k[0],   k[1], k[2], k[3],   k[4], k[5]+1,k[6])
      val = coef * poly[k]
      poln.add_keyval(key, val)
  if tp == 4:
    for k in poly:
      key = (k[0],   k[1], k[2], k[3],   k[4], k[5],  k[6]+1)
      val = coef * poly[k]
      poln.add_keyval(key, val)
  return(poln)

def recursive_eII(n3, n4, eII):
  '''Recursive step in equation :eq:`ito-moment-i`
  
  :param n3: :math:`n_3` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :param n4: :math:`n_4` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :param eII: a dict with key (n3,n4) and value Poly object with attribute 
     ``keyfor`` = 
     ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  
  :return: updated eII.
  :rtype: dict
  '''
  poly = Poly()
  kf = ['e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 
        'v_{n-1}', 'k^{-}', 'theta', 'sigma_v']
  poly.set_keyfor(kf)
  #
  if n3 >= 2 and n4 >= 0:
    c = Frac(n3*(n3-1), 2)
    poly.merge(coef_poly(c, int_meII(1, n3-2, n4, eII), 1))
    poly.merge(coef_poly(c, int_meII(1, n3-2, n4, eII), 2))
    poly.merge(coef_poly(c, int_meII(2, n3-2, n4, eII), 3))
    poly.merge(coef_poly(c, int_meII(1, n3-1, n4, eII), 4))
  if n3 >= 0 and n4 >= 2:
    c = Frac(n4*(n4-1), 2)
    poly.merge(coef_poly(c, int_meII(-1, n3,   n4-2, eII), 1))
    poly.merge(coef_poly(c, int_meII(-1, n3,   n4-2, eII), 2))
    poly.merge(coef_poly(c, int_meII( 0, n3,   n4-2, eII), 3))
    poly.merge(coef_poly(c, int_meII(-1, n3+1, n4-2, eII), 4))
  if n3 >= 1 and n4 >= 1:
    c = Frac(n3*n4, 1)
    poly.merge(coef_poly(c, int_meII(0, n3-1, n4-1, eII), 1))
    poly.merge(coef_poly(c, int_meII(0, n3-1, n4-1, eII), 2))
    poly.merge(coef_poly(c, int_meII(1, n3-1, n4-1, eII), 3))
    poly.merge(coef_poly(c, int_meII(0, n3,   n4-1, eII), 4))
  return(poly)

def recursive_eIII(n3, n4, n5, eIII):
  '''Recursive step in equation :eq:`ito-moment-ii`
  
  :param n3: :math:`n_3` 
     in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :param n4: :math:`n_4` 
     in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :param n5: :math:`n_5` 
     in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :param eIII: a dict with key (n3,n4,n5) and value Poly object with attribute 
     ``keyfor`` = 
     ('e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 'v_{n-1}',
     'k^{-}', 'theta', 'sigma_v').
  
  :return: updated eIII.
  :rtype: dict
  '''
  poly = Poly()
  kf = ['e^{k(n-1)h}', 'e^{k[t-(n-1)h]}', '[t-(n-1)h]', 
        'v_{n-1}', 'k^{-}', 'theta', 'sigma_v']
  poly.set_keyfor(kf)
  #
  if n3 >= 2 and n4 >=0 and n5 >= 0:
    c = Frac(n3*(n3-1), 2)
    poly.merge(coef_poly(c, int_meIII(1, n3-2, n4, n5, eIII), 1))
    poly.merge(coef_poly(c, int_meIII(1, n3-2, n4, n5, eIII), 2))
    poly.merge(coef_poly(c, int_meIII(2, n3-2, n4, n5, eIII), 3))
    poly.merge(coef_poly(c, int_meIII(1, n3-1, n4, n5, eIII), 4))
  if n3 >= 0 and n4 >= 2 and n5 >= 0:
    c = Frac(n4*(n4-1), 2)
    poly.merge(coef_poly(c, int_meIII(-1, n3, n4-2, n5, eIII), 1))
    poly.merge(coef_poly(c, int_meIII(-1, n3, n4-2, n5, eIII), 2))
    poly.merge(coef_poly(c, int_meIII( 0, n3, n4-2, n5, eIII), 3))
    poly.merge(coef_poly(c, int_meIII(-1,n3+1,n4-2, n5, eIII), 4))
  if n3 >= 1 and n4 >= 1 and n5 >= 0:
    c = Frac(n3*n4, 1)
    poly.merge(coef_poly(c, int_meIII(0, n3-1, n4-1, n5, eIII), 1))
    poly.merge(coef_poly(c, int_meIII(0, n3-1, n4-1, n5, eIII), 2))
    poly.merge(coef_poly(c, int_meIII(1, n3-1, n4-1, n5, eIII), 3))
    poly.merge(coef_poly(c, int_meIII(0, n3,   n4-1, n5, eIII), 4))
  if n3 >= 0 and n4 >= 0 and n5 >= 2:
    c = Frac(n5*(n5-1), 2)
    poly.merge(coef_poly(c, int_meIII(-1, n3, n4, n5-2, eIII), 1))
    poly.merge(coef_poly(c, int_meIII(-1, n3, n4, n5-2, eIII), 2))
    poly.merge(coef_poly(c, int_meIII( 0, n3, n4, n5-2, eIII), 3))
    poly.merge(coef_poly(c, int_meIII(-1,n3+1,n4, n5-2, eIII), 4))
  return(poly)

def moment_eII(n3, n4, return_all = False):
  ''':math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`
  
  :param n3: :math:`n_3` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :param n4: :math:`n_4` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :param return_pre: whether or not return lower order moments simultaneously,
     default to ``False``.
  
  :return: poly if return_all=False else eII, where poly with attribute 
     ``keyfor`` = 
     ('e^{k(n-1)h}','e^{k[t-(n-1)h]}','[t-(n-1)h]','v_{n-1}','k^{-}','theta',
     'sigma_v').
  :rtype: Poly or dict of Poly
  '''
  # eII: dict of E[eI^{n3}I^{n4}]
  eII = {}
  # n3 + n4 = 0
  # support for special case
  poly = Poly({(0,0,0,0,0,0,0): Frac(1,1)}) # num/den
  kf = ['e^{k(n-1)h}','e^{k[t-(n-1)h]}','[t-(n-1)h]','v_{n-1}','k^{-}',
    'theta','sigma_v']
  poly.set_keyfor(kf)
  eII[(0,0)] = poly
  # n3 + n4 = 1
  eII[(1,0)] = {}
  eII[(0,1)] = {}
  # n3 + n4 = 2
  poly = Poly({(2,2,0,0,1,1,0): Frac(1,2), 
               (2,1,0,1,1,0,0): Frac(1,1),
               (2,1,0,0,1,1,0): Frac(-1,1),
               (2,0,0,1,1,0,0): Frac(-1,1),
               (2,0,0,0,1,1,0): Frac(1,2)})
  poly.set_keyfor(kf)
  eII[(2,0)] = poly
  poly = Poly({(1,1,0,0,1,1,0): Frac(1,1),
               (1,0,1,1,0,0,0): Frac(1,1),
               (1,0,1,0,0,1,0): Frac(-1,1),
               (1,0,0,0,1,1,0): Frac(-1,1)})
  poly.set_keyfor(kf)
  eII[(1,1)] = poly
  poly = Poly({(0,-1,0,1,1,0,0): Frac(-1,1),
               (0,-1,0,0,1,1,0): Frac(1,1),
               (0, 0,1,0,0,1,0): Frac(1,1),
               (0, 0,0,1,1,0,0): Frac(1,1),
               (0, 0,0,0,1,1,0): Frac(-1,1)})
  poly.set_keyfor(kf)
  eII[(0,2)] = poly
  #
  if n3 + n4 <= 2:
    return( eII if return_all else eII[(n3,n4)] )
  #
  if n3 + n4 > 3:
    # compute all lower-order moments to get ready
    for n in range(3, n3+n4):
      for i in range(n, -1, -1):
        poly = recursive_eII(i, n-i, eII)
        poly.remove_zero()
        eII[(i,n-i)] = poly
  # the last one
  poly = recursive_eII(n3, n4, eII)
  poly.remove_zero()
  eII[(n3,n4)] = poly
  return( eII if return_all else poly )

def moment_eIII(n3, n4, n5, return_all = False):
  ''':math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`
  
  :param n3: :math:`n_3` 
     in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :param n4: :math:`n_4` 
     in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :param n5: :math:`n_5` 
     in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :param return_all: whether or not return lower order moments simultaneously,
     default to ``False``.
  
  :return: poly if return_all=False else eIII, where poly with attribute 
     ``keyfor`` = 
     ('e^{k(n-1)h}','e^{k[t-(n-1)h]}','[t-(n-1)h]','v_{n-1}','k^{-}','theta',
     'sigma_v').
  :rtype: dict or dict of dict
  '''
  # eIII: a dict of moments of E[eI^{n3}I^{n4}I^{*n5}]
  if n3 + n4 + n5 < 0:
    raise ValueError(f"moment_eIII({n3},{n4},{n5}) is called!")
  eIII = {}
  # n3 + n4 + n5 = 0
  # support for special case
  poly = Poly({(0,0,0,0,0,0,0): Frac(1, 1)})
  kf = ['e^{k(n-1)h}','e^{k[t-(n-1)h]}','[t-(n-1)h]','v_{n-1}','k^{-}',
    'theta','sigma_v']
  poly.set_keyfor(kf)
  eIII[(0,0,0)] =  poly
  # n3 + n4 + n5 = 1
  eIII[(1,0,0)] = {}
  eIII[(0,1,0)] = {}
  eIII[(0,0,1)] = {}
  # n3 + n4 + n5 = 2
  eIII[(2,0,0)] = moment_eII(2,0)
  eIII[(1,1,0)] = moment_eII(1,1)
  eIII[(1,0,1)] = {}
  eIII[(0,2,0)] = moment_eII(0,2)
  eIII[(0,1,1)] = {}
  eIII[(0,0,2)] = moment_eII(0,2)
  #
  if n3 + n4 + n5 <= 2:
    return(eIII if return_all else eIII[(n3,n4,n5)])
  #
  if n3 + n4 + n5 > 3:
    # compute all lower-order moments to get ready
    for n in range(3, n3+n4+n5):
      for i in range(n, -1, -1):
        for j in range(n-i, -1, -1):
          poly = recursive_eIII(i, j, n-i-j, eIII)
          poly.remove_zero()
          eIII[(i,j,n-i-j)] = poly
  # the last one
  poly = recursive_eIII(n3, n4, n5, eIII)
  poly.remove_zero()
  eIII[(n3,n4,n5)] = poly
  return(eIII if return_all else poly)

def moment_v(n):
  '''Moment of :math:`v_{n-1}` as in equation :eq:`moment-v`
  
  :param n: order of the moment.
  
  :return: a poly with attribute ``keyfor`` = ('theta','sigma_v^2/k').
  :rtype: Poly
  '''
  v = []
  kf = ['theta','sigma_v^2/k']
  # n = 0
  poly = Poly({(0,0): Frac(1,1)}); poly.set_keyfor(kf)
  v.append(poly)
  if n == 0: return(poly)
  # n = 1
  poly = Poly({(1,0): Frac(1,1)})
  poly.set_keyfor(kf); v.append(poly)
  if n == 1: return(poly)
  # n >= 2
  for i in range(2, n+1):
    poly = v[-1] # recursively computing
    poln = Poly(); poln.set_keyfor(kf)
    for k in poly:
      # times theta and ((i-1)/2) sigma_v^2/k
      poln.add_keyval((k[0]+1, k[1]),   poly[k])
      poln.add_keyval((k[0],   k[1]+1), poly[k] * Frac(i-1, 2))
    v.append(poln)
  return(poln)

if __name__ == "__main__":
  # test the module
  from pprint import pprint
  # 
  print("moment_v(n): poly with keyfor = ('theta','sigma_v^2/k')")
  # print('moment_v(0): '); pprint(moment_v(0))
  # print('moment_v(1): '); pprint(moment_v(1))
  # print('moment_v(2): '); pprint(moment_v(2))
  print('moment_v(3): '); pprint(moment_v(3))
  # print('moment_v(4): '); pprint(moment_v(4))
  # 
  kf = ('e^{k(n-1)h}','e^{k[t-(n-1)h]}','[t-(n-1)h]','v_{n-1}','k^{-}',
    'theta','sigma_v')
  print(f"moment_eII(n3,n4): poly with keyfor = {kf}")
  # print('moment_eII(0,0): '); pprint(moment_eII(0,0)) # verified
  # print('moment_eII(1,0): '); pprint(moment_eII(1,0)) #->
  # print('moment_eII(0,1): '); pprint(moment_eII(0,1))
  # print('moment_eII(2,0): '); pprint(moment_eII(2,0))
  # print('moment_eII(1,1): '); pprint(moment_eII(1,1))
  # print('moment_eII(0,2): '); pprint(moment_eII(0,2))
  print('moment_eII(3,0): '); pprint(moment_eII(3,0))
  print('moment_eII(2,1): '); pprint(moment_eII(2,1))
  print('moment_eII(1,2): '); pprint(moment_eII(1,2)) #<-
  print('moment_eII(0,3): '); pprint(moment_eII(0,3)) # verified
  # print('moment_eII(4,0): '); pprint(moment_eII(4,0)) # verified 1st cmpnt
  # print('moment_eII(3,1): '); pprint(moment_eII(3,1))
  # print('moment_eII(2,2): '); pprint(moment_eII(2,2))
  # print('moment_eII(1,3): '); pprint(moment_eII(1,3))
  # print('moment_eII(0,4): '); pprint(moment_eII(0,4)) # verified 1st cmpnt 
  # print('moment_eII(5,0): '); pprint(moment_eII(5,0))
  # print('moment_eII(1,6): '); pprint(moment_eII(1,6))
  # print('moment_eII(2,7): '); pprint(moment_eII(2,7))
  # print('moment_eII(10,3): '); pprint(moment_eII(10,3))
  #
  print(f"moment_eIII(n3,n4,n5): poly with keyfor = {kf}")
  # 
  print('moment_eIII(0,0,0): '); pprint(moment_eIII(0,0,0))
  # 
  print('moment_eIII(1,0,0): '); pprint(moment_eIII(1,0,0))
  print('moment_eIII(0,1,0): '); pprint(moment_eIII(0,1,0))
  print('moment_eIII(0,0,1): '); pprint(moment_eIII(0,0,1))
  # 
  print('moment_eIII(2,0,0): '); pprint(moment_eIII(2,0,0)) # moment_eII(2,0)
  print('moment_eIII(1,1,0): '); pprint(moment_eIII(1,1,0)) # moment_eII(1,1)
  print('moment_eIII(1,0,1): '); pprint(moment_eIII(1,0,1))
  print('moment_eIII(0,2,0): '); pprint(moment_eIII(0,2,0))
  print('moment_eIII(0,1,1): '); pprint(moment_eIII(0,1,1))
  print('moment_eIII(0,0,2): '); pprint(moment_eIII(0,0,2))
  # 
  print('moment_eIII(3,0,0): '); pprint(moment_eIII(3,0,0))
  print('moment_eIII(2,1,0): '); pprint(moment_eIII(2,1,0))
  print('moment_eIII(2,0,1): '); pprint(moment_eIII(2,0,1))
  print('moment_eIII(1,2,0): '); pprint(moment_eIII(1,2,0))
  print('moment_eIII(1,1,1): '); pprint(moment_eIII(1,1,1))
  print('moment_eIII(1,0,2): '); pprint(moment_eIII(1,0,2))
  print('moment_eIII(0,3,0): '); pprint(moment_eIII(0,3,0))
  print('moment_eIII(0,2,1): '); pprint(moment_eIII(0,2,1))
  print('moment_eIII(0,1,2): '); pprint(moment_eIII(0,1,2))
  print('moment_eIII(0,0,3): '); pprint(moment_eIII(0,0,3))
  # 
  print('moment_eIII(4,0,0): '); pprint(moment_eIII(4,0,0))
  print('moment_eIII(5,1,2): '); pprint(moment_eIII(5,1,2))
  print('moment_eIII(7,2,3): '); pprint(moment_eIII(7,2,3)) # {}
  print('moment_eIII(7,0,2): '); pprint(moment_eIII(7,0,2))
  print('moment_eIII(4,2,5): '); pprint(moment_eIII(4,2,5)) # {}

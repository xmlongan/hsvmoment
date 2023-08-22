'''
Module mdl_2FSVJ

Functions to compute moments of the Two Factor SV model with jumps in return 
process.
'''
import math
from ito_moments import merge, from_nh_to_np1h
from simplify_terms import times_pow, times_exp
from cpp_moments import mcpp

def c_n1n2mi(n1,n2,m,i):
  '''Coefficent :math:`c_{n1n2mi}` as in :eq:`c-n_1n_2mi`.
  
  :param n1: :math:`n_1` in :math:`e^{(n_1k_1+n_2k_2)t}`.
  :type integer: int
  :param n2: :math:`n_2` in :math:`e^{(n_1k_1+n_2k_2)t}`.
  :type integer: int
  :param m: m in :math:`t^{m-i}`.
  :type integer: int
  :param i: i in :math:`t^{m-i}`.
  :type integer: int
  
  :return: the coefficient :math:`c_{n1n2mi}`.
  :rtype: string
  '''
  prod = 1
  for j in range(m-i+1, m+1):
    prod = prod * j
  den = f'({n1}*k_1 + {n2}*k_2)^({i+1})' if i+1 != 1 else f'({n1}*k_1 + {n2}*k_2)'
  c = f'{(-1)**i * prod}/{den}'
  return(c)

def int_ennt(n1,n2,m):
  ''':math:`\int_{(n-1)h}^t e^{(n_1k_1+n_2k_2)s}s^mds=\sum_{ij}c_{ii'j}e^{(ik_1+i'k_2)t}t^j`
  
  :param n1: :math:`n_1` in :math:`e^{(n_1k_1+n_2k_2)s}s^m`.
  :type integer: int
  :param n2: :math:`n_2` in :math:`e^{(n_1k_1+n_2k_2)s}s^m`.
  :type integer: int
  :param m: m in :math:`e^{(n_1k_1+n_2k_2)s}s^m`.
  :type integer: int
  
  :return: a dictionary with key-value pairs (i,i',j)-c_ij'j, 
   representing :math:`\sum_{ii'j}c_{ii'j}e^{(ik_1+i'k_2)t}t^j`.
  :rtype: dict of tuple: str
  '''
  if m < 0:
    msg = f"m in int_et(n,m) equals {m}, however it shouldn't be negative!"
    raise ValueError(msg)
  poly = {}
  if n1 + n2 == 0 and m == 0:
    poly[(0,0,1)] = '1'
    poly[(0,0,0)] = '-(n-1)*h'
  elif n1 + n2 == 0 and m != 0:
    poly[(0,0,m+1)] = f'1/({m+1})'
    poly[(0,0,0)] = f'-((n-1)*h)^({m+1})/({m+1})'
  elif n1 + n2 != 0 and m == 0:
    num = f'exp(({n1}*k_1+{n2}*k_2)*(n-1)*h)'
    den = f'({n1}*k_1+{n2}*k_2)'
    poly[(n1,n2,0)] = f'1/{den}'
    poly[(0,0,0)] = f' - {num}/{den}'
  else:
    den = f'({n1}*k_1+{n2}*k_2)'
    num = f'exp(({n1}*k_1+{n2}*k_2)*(n-1)*h)'
    num1   = num + times_pow('(n-1)*h', m)
    poly[(n1,n2,m)] = f'1/{den}'
    const    = f' - {num1}/{den}'
    for i in range(1, m+1):
      num2 = num + times_pow('(n-1)*h', m-i)
      poly[(n1,n2,m-i)] = c_n1n2mi(n1,n2,m,i)
      const += f' - {num2}*({poly[(n1,n2,m-i)]})' 
    poly[(0,0,0)] = const
  return(poly)

def int_meIIeIII(i, m, n4, n5, n6, n7, n8, eIIeIII):
  ''':math:`\int_{(n-1)h}^te^{mk_is}E[n_4n_5n_6n_7n_8]ds`
  
  :param i: :math:`i` in :math:`e^{mk_is} (i=1,2)`.
  :param m: :math:`m` in :math:`e^{mk_is} (i=1,2)`.
  :param n4: :math:`m_4` in :math:`eI_{1,n-1,t}^{m_4}`.
  :param n5: :math:`m_5` in :math:`I_{1,n-1,t}^{m_5}`.
  :param n6: :math:`m_6` in :math:`eI_{2,n-1,t}^{m_6}`.
  :param n7: :math:`m_7` in :math:`I_{2,n-1,t}^{m_7}`.
  :param n8: :math:`m_8` in :math:`I_{1,n-1,t}^{*m_8}`.
  :param eIIeIII: a dict with key (n4,n5,n6,n7,n8) and value poly with 
     key(i,i',j,l,l') and value :math:`b_{ii'jll'}`.
  
  :return: the integral result in poly with key (i,i',j,l,l') with 
     value :math:`b_{ii'jll'}`.
  :rtype: dict
  '''
  # get the representation of E[n_4n_5n_6n_7n_8]
  b = eIIeIII[(n4,n5,n6,n7,n8)]
  #
  poly_sum = {}
  for k in b.keys():
    c = int_ennt(k[0]+m, k[1], k[2]) if i == 1 else int_ennt(k[0], k[1]+m, k[2])
    # c: poly dict {(i,j,l): e^{(ik_1+jk_2)t} t^l}
    for ijl in c.keys():
      coef = f'({b[k]}) * ({c[ijl]})'
      key = (ijl[0], ijl[1], ijl[2], k[3], k[4])
      if key in poly_sum:
        poly_sum[key] += f' + {coef}'
      else:
        poly_sum[key]  = coef
  return(poly_sum)

def get_coef(prefix,i):
  '''Get the coefficients in Recursive Equation :eq:`ito-moment-m4m5m6m7m8`
  
  :param prefix: the leading term, such as :math:`\\frac{m_4(m_4-1)}{2}`.
  :param i: :math:`i` in :math:`v_i (i=1,2)`.
  
  :return: a tuple of four coefficients.
  :rtype: a tuple of four str
  '''
  c1 = prefix + f'exp(k{i}*(n-1)*h)'
  c2 = '-' + prefix + f'theta{i}*exp(k{i}*(n-1)*h)'
  c3 = prefix + f'theta{i}'
  c4 = prefix + f'sigma_v{i}'
  return((c1,c2,c3,c4))

def coef_poly(coef, poly):
  '''Multiply poly with coefficient
  
  :param coef: coefficient to multiply by poly.
  :type string: str
  :param poly: one polynomial of
     :math:`\sum_{i,i',j,l,l'} b_{ii'jll'} e^{(ik_1 + i'k_2)t} t^j v_{1,n-1}^{l}v_{2,n-1}^{l'}`.
  :type poly: dict with key (i,i',j,l,l') and value :math:`b_{ii'jll'}`
  
  :return: new polynomial.
  :rtype: dict with key (i,i',j,l,l') and value :math:`b_{ii'jll'}`
  '''
  for key in poly.keys():
    poly[key] = f'{coef} * ({poly[key]})'
  return(poly)

def coek_poly(coef, poly, i):
  '''Multiply poly with coefficient and update key
  
  :param coef: coefficient to multiply by poly.
  :type string: str
  :param poly: one polynomial of
     :math:`\sum_{i,i',j,l,l'} b_{ii'jll'} e^{(ik_1 + i'k_2)t} t^j v_{1,n-1}^{l}v_{2,n-1}^{l'}`.
  :type poly: dict with key (i,i',j,l,l') and value :math:`b_{ii'jll'}`
  :param i: :math:`i` in :math:`v_i (i=1,2)`.
  
  :return: new polynomial.
  :rtype: dict with key (i,i',j,l,l') and value :math:`b_{ii'jll'}`
  '''
  poly_new = {}
  for k in poly.keys():
    if i == 1:
      key = (k[0],k[1],k[2],k[3]+1,k[4])
    else:
      key = (k[0],k[1],k[2],k[3],k[4]+1)
    poly_new[key] = f'{coef} * ({poly[k]})'
  return(poly_new)

def recursive_eIIeIII(n4, n5, n6, n7, n8, eIIeIII):
  '''Recursive step of :math:`E[m_4m_5m_6m_7m_8]` as in :eq:`ito-moment-m4m5m6m7m8`
  
  :param n4: :math:`m_4` in :math:`eI_{1,n-1,t}^{m_4}`.
  :param n5: :math:`m_5` in :math:`I_{1,n-1,t}^{m_5}`.
  :param n6: :math:`m_6` in :math:`eI_{2,n-1,t}^{m_6}`.
  :param n7: :math:`m_7` in :math:`I_{2,n-1,t}^{m_7}`.
  :param n8: :math:`m_8` in :math:`I_{1,n-1,t}^{*m_8}`.
  :param eIIeIII: a dict with key (n4,n5,n6,n7,n8) and value poly with 
     key(i,i',j,l,l') and value :math:`b_{ii'jll'}`.
  
  :return: updated eIIeIII.
  :rtype: dict
  '''
  poly_sum = {}
  #
  if n4 >= 2 and n5 >= 0 and n6 >= 0 and n7 >= 0 and n8 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n4*(n4-1)/2)}*' if n4 > 2 else '', 1)
    p1 = coev_poly(c1, int_meIIeIII(1, 1, n4-2, n5, n6, n7, n8, eIIeIII), 1)
    p2 = coef_poly(c2, int_meIIeIII(1, 1, n4-2, n5, n6, n7, n8, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(1, 2, n4-2, n5, n6, n7, n8, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(1, 1, n4-1, n5, n6, n7, n8, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  elif n4 >= 0 and n5 >=2 and n6 >= 0 and n7 >= 0 and n8 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n5*(n5-1)/2)}*' if n5 > 2 else '', 1)
    p1 = coev_poly(c1, int_meIIeIII(1, -1, n4,  n5-2, n6, n7, n8, eIIeIII), 1)
    p2 = coef_poly(c2, int_meIIeIII(1, -1, n4,  n5-2, n6, n7, n8, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(1,  0, n4,  n5-2, n6, n7, n8, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(1, -1, n4+1,n5-2, n6, n7, n8, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  elif n4 >= 1 and n5 >= 1 and n6 >= 0 and n7 >= 0 and n8 >= 0:
    c1,c2,c3,c4 = get_coef(f'{n4*n5}*' if n4*n5 != 1 else '', 1)
    p1 = coev_poly(c1, int_meIIeIII(1, 0, n4-1, n5-1, n6, n7, n8, eIIeIII), 1)
    p2 = coef_poly(c2, int_meIIeIII(1, 0, n4-1, n5-1, n6, n7, n8, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(1, 1, n4-1, n5-1, n6, n7, n8, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(1, 0, n4,   n5-1, n6, n7, n8, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  elif n4 >= 0 and n5 >= 0 and n6 >= 2 and n7 >= 0 and n8 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n6*(n6-1)/2)}*' if n6 > 2 else '', 2)
    p1 = coev_poly(c1, int_meIIeIII(2, 1, n4, n5, n6-2, n7, n8, eIIeIII), 2)
    p2 = coef_poly(c2, int_meIIeIII(2, 1, n4, n5, n6-2, n7, n8, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(2, 2, n4, n5, n6-2, n7, n8, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(2, 1, n4, n5, n6-1, n7, n8, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  elif n4 >= 0 and n5 >= 0 and n6 >= 0 and n7 >= 2 and n8 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n7*(n7-1)/2)}*' if n7 > 2 else '', 2)
    p1 = coev_poly(c1, int_meIIeIII(2, -1, n4, n5, n6, n7-2, n8, eIIeIII), 2)
    p2 = coef_poly(c2, int_meIIeIII(2, -1, n4, n5, n6, n7-2, n8, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(2,  0, n4, n5, n6, n7-2, n8, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(2, -1, n4, n5,n6+1,n7-2, n8, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  elif n4 >= 0 and n5 >= 0 and n6 >= 1 and n7 >= 1 and n8 >= 0:
    c1,c2,c3,c4 = get_coef(f'{n6*n7}*' if n6*n7 != 1 else '', 2)
    p1 = coev_poly(c1, int_meIIeIII(2, 0, n4, n5, n6-1, n7-1, n8, eIIeIII), 2)
    p2 = coef_poly(c2, int_meIIeIII(2, 0, n4, n5, n6-1, n7-1, n8, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(2, 1, n4, n5, n6-1, n7-1, n8, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(2, 0, n4, n5, n6,   n7-1, n8, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  elif n4 >= 0 and n5 >= 0 and n6 >= 0 and n7 >= 0 and n8 >= 2:
    c1,c2,c3,c4 = get_coef(f'{int(n8*(n8-1)/2)}*' if n8 > 2 else '', 1)
    p1 = coev_poly(c1, int_meIIeIII(1, -1, n4, n5, n6, n7, n8-2, eIIeIII), 1)
    p2 = coef_poly(c2, int_meIIeIII(1, -1, n4, n5, n6, n7, n8-2, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(1,  0, n4, n5, n6, n7, n8-2, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(1, -1,n4+1,n5, n6, n7, n8-2, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
    #
    c1,c2,c3,c4 = get_coef(f'{int(n8*(n8-1)/2)}*' if n8 > 2 else '', 2)
    p1 = coev_poly(c1, int_meIIeIII(2, -1, n4, n5, n6, n7, n8-2, eIIeIII), 2)
    p2 = coef_poly(c2, int_meIIeIII(2, -1, n4, n5, n6, n7, n8-2, eIIeIII))
    p3 = coef_poly(c3, int_meIIeIII(2,  0, n4, n5, n6, n7, n8-2, eIIeIII))
    p4 = coef_poly(c4, int_meIIeIII(2, -1, n4, n5,n6+1,n7, n8-2, eIIeIII))
    poly_sum = merge(p4, merge(p3, merge(p2, merge(p1, poly_sum))))
  #
  return(poly_sum)

def moment_eIIeIII(n4, n5, n6, n7, n8, return_all=False):
  ''':math:`E[eI_{1,n-1,t}^{m_4} I_{1,n-1,t}^{m_5} eI_{2,n-1,t}^{m_6} I_{2,n-1,t}^{m_7}I_{n-1,t}^{*m_8}|v_{1,n-1},v_{2,n-1}]`
  
  Moment of the combination of Ito processes, conditioning on the previous 
  volatility components.
  
  :param n4: :math:`m_4` in :math:`eI_{1,n-1,t}^{m_4}`.
  :param n5: :math:`m_5` in :math:`I_{1,n-1,t}^{m_5}`.
  :param n6: :math:`m_6` in :math:`eI_{2,n-1,t}^{m_6}`.
  :param n7: :math:`m_7` in :math:`I_{2,n-1,t}^{m_7}`.
  :param n8: :math:`m_8` in :math:`I_{1,n-1,t}^{*m_8}`.
  :param return_all: whether or not return lower order moments simultaneously,
     defaults to ``False``.
  
  :return: poly if return_all=False else eIIeIII
  :rtype: dict or dict of dict
  '''
  # eIIeIII: a dict of moments of E[eI_{1,n-1,t}^{m_4} I_{1,n-1,t}^{m_5}
  # eI_{2,n-1,t}^{m_6} I_{2,n-1,t}^{m_7}I_{n-1,t}^{*m_8}|v_{1,n-1},v_{2,n-1}]
  if n4+n5+n6+n7+n8 < 0:
    raise ValueError(f'moment_eIIeIII({n_4},{n_5},{n_6},{n_7},{n_8}) is called')
  eIIeIII = {}
  # n4+n5+n6+n7+n8 = 0
  eIIeIII[(0,0,0,0,0)] = {(0,0,0,0,0): '1'} # support for special case
  # n4+n5+n6+n7+n8 = 1
  eIIeIII[(1,0,0,0,0)] = {}
  eIIeIII[(0,1,0,0,0)] = {}
  eIIeIII[(0,0,1,0,0)] = {}
  eIIeIII[(0,0,0,1,0)] = {}
  eIIeIII[(0,0,0,0,1)] = {}
  # n4+n5+n6+n7+n8 = 2
  eIIeIII[(2,0,0,0,0)] = {(2,0,0,0,0): 'theta1/(2*k1)',
                          (1,0,0,1,0): 'exp(k1*(n-1)*h)/k1',
                          (1,0,0,0,0): '-exp(k1*(n-1)*h)*theta1/k1',
                          (0,0,0,1,0): '-exp(2*k1*(n-1)*h)/k1',
                          (0,0,0,0,0): 'exp(2*k1*(n-1)*h)*theta1/(2*k1)'}
  eIIeIII[(1,1,0,0,0)] = {(1,0,0,0,0): 'theta1/k1',
                          (0,0,1,1,0): 'exp(k1*(n-1)*h)',
                          (0,0,1,0,0): '-exp(k1*(n-1)*h)*theta1',
                          (0,0,0,1,0): '-exp(k1*(n-1)*h)*(n-1)*h',
                          (0,0,0,0,0): 'exp(k1*(n-1)*h)*((n-1)*h*theta1-theta1/k1)'}
  eIIeIII[(1,0,1,0,0)] = {}
  eIIeIII[(1,0,0,1,0)] = {}
  eIIeIII[(1,0,0,0,1)] = {}
  eIIeIII[(0,2,0,0,0)] = {(-1,0,0,1,0): '-exp(k1*(n-1)*h)/k1',
                          (-1,0,0,0,0): 'exp(k1*(n-1)*h)*theta1/k1',
                          ( 0,0,1,0,0): 'theta1',
                          ( 0,0,0,1,0): '1/k1',
                          ( 0,0,0,0,0): '-(1/k1+(n-1)*h)*theta1'}
  eIIeIII[(0,1,1,0,0)] = {}
  eIIeIII[(0,1,0,1,0)] = {}
  eIIeIII[(0,1,0,0,1)] = {}
  eIIeIII[(0,0,2,0,0)] = {(0,2,0,0,0): 'theta2/(2*k2)',
                          (0,1,0,0,1): 'exp(k2*(n-1)*h)/k2',
                          (0,1,0,0,0): '-exp(k2*(n-1)*h)*theta2/k2',
                          (0,0,0,0,1): '-exp(2*k2*(n-1)*h)/k2',
                          (0,0,0,0,0): 'exp(2*k2*(n-1)*h)*theta2/(2*k2)'}
  eIIeIII[(0,0,1,1,0)] = {(0,1,0,0,0): 'theta2/k2',
                          (0,0,1,0,1): 'exp(k2*(n-1)*h)',
                          (0,0,1,0,0): '-exp(k2*(n-1)*h)*theta2',
                          (0,0,0,0,1): '-exp(k2*(n-1)*h)*(n-1)*h',
                          (0,0,0,0,0): 'exp(k2*(n-1)*h)*theta2*((n-1)*h-1/k2)'}
  eIIeIII[(0,0,1,0,1)] = {}
  eIIeIII[(0,0,0,2,0)] = {(0,-1,0,0,1): '-exp(k2*(n-1)*h)/k2',
                          (0,-1,0,0,0): 'exp(k2*(n-1)*h)*theta2/k2',
                          (0, 0,1,0,0): 'theta2',
                          (0, 0,0,0,1): 'k2',
                          (0, 0,0,0,0): '-((n-1)*h + 1/k2)*theta2'}
  eIIeIII[(0,0,0,1,1)] = {}
  eIIeIII[(0,0,0,0,2)] = merge(eIIeIII[(0,2,0,0,0)], eIIeIII[(0,0,0,2,0)])
  #
  if n4+n5+n6+n7+n8 <= 2:
    return( eIIeIII if return_all else eIIeIII[(n4,n5,n6,n7,n8)])
  #
  if n4+n5+n6+n7+n8 > 3:
    # compute all lower-order moments to get ready for the last computation
    for n in range(3, n4+n5+n6+n7+n8):
      for i1 in range(n, -1, -1):
        for i2 in range(n-i1, -1, -1):
          for i3 in range(n-i1-i2, -1, -1):
            for i4 in range(n-i1-i2-i3, -1, -1):
              i5 = n - i1 - i2 - i3 - i4
              poly = recursive_eIIeIII(i1,i2,i3,i4,i5, eIIeIII)
              eIIeIII[(i1,i2,i3,i4,i5)] = poly
  #
  poly = recursive_eIIeIII(n4,n5,n6,n7,n8, eIIeIII)
  eIIeIII[(n4,n5,n6,n7,n8)] = poly
  return(eIIeIII if return_all else poly)

###

###

def c_n(l, n1, n2, n3, n4, n5, n6, n7, n8):
  '''Number of occurrences of this combination
  
  :param l: l in :math:`E[\overline{y}_{n}^l]`.
  :param n1: times of constant being selected.
  :param n2: times of :math:`v_{1,n-1}` being selected.
  :param n3: times of :math:`v_{2,n-1}` being selected.
  :param n4: times of :math:`eI_{1,n}` being selected.
  :param n5: times of :math:`I_{1,n}` being selected.
  :param n6: times of :math:`eI_{2,n}` being selected.
  :param n7: times of :math:`I_{2,n}` being selected.
  :param n8: times of :math:`I_{n}^{*}` being selected.
  
  :return: number of occurrences of this combination.
  :rtype: int
  '''
  num  = math.comb(l, n1)
  num *= math.comb(l-n1, n2)
  num *= math.comb(l-n1-n2, n3)
  num *= math.comb(l-n1-n2-n3, n4)
  num *= math.comb(l-n1-n2-n3-n4, n5)
  num *= math.comb(l-n1-n2-n3-n4-n5, n6)
  num *= math.comb(l-n1-n2-n3-n4-n5-n6, n7)
  num *= math.comb(l-n1-n2-n3-n4-n5-n6-n7, n8)
  return(num)

def b_n(n1, n2, n3, n4, n5, n6, n7):
  '''Coefficient corresponding to the combination
  
  :param n1: times of constant being selected.
  :param n2: times of :math:`v_{1,n-1}` being selected.
  :param n3: times of :math:`v_{2,n-1}` being selected.
  :param n4: times of :math:`eI_{1,n}` being selected.
  :param n5: times of :math:`I_{1,n}` being selected.
  :param n6: times of :math:`eI_{2,n}` being selected.
  :param n7: times of :math:`I_{2,n}` being selected.
  
  :return: coefficient of this combination.
  :rtype: str
  '''
  b  = '1' if (n2+n3+n5+n7)%2 == 0 else '(-1)'
  h1 = '(1-exp(-k1*h))/k1'
  h2 = '(1-exp(-k2*h))/k2'
  constant = f'-(lambda*h)*mu_j + (theta1*{h1} + theta2*{h2})/2'
  b += times_pow(f'{constant}', n1)
  b += times_pow(f'{h1}/2', n2)
  b += times_pow(f'{h2}/2', n3)
  b += times_pow('sigma_v1/(2*k1)', n4+n5)
  b += times_pow('sigma_v2/(2*k2)', n6+n7)
  return(b)

def moment_v(n, i):
  '''Moment of :math:`v_{i,n-1}, i=1,2`
  
  :param n: order of the moments.
  
  :return: n order moment of :math:`v_i`.
  :rtype: str
  '''
  if n == 0:
    return('1')
  for j in range(n):
    if j == 0:
      moment = f'theta{i}'
    else:
      moment = f'{moment}*(theta{i}+{j}*sigma_v{i}^2/(2*k{i}))'
  return(moment)

def moment_comb(l, n1, n2, n3, n4, n5, n6, n7, n8):
  '''Moment for this combination in expansion of :math:`\overline{y}_n^l`
  
  :param l: l in :math:`E[\overline{y}_{n}^l]`.
  :param n1: times of constant being selected.
  :param n2: times of :math:`v_{1,n-1}` being selected.
  :param n3: times of :math:`v_{2,n-1}` being selected.
  :param n4: times of :math:`eI_{1,n}` being selected.
  :param n5: times of :math:`I_{1,n}` being selected.
  :param n6: times of :math:`eI_{2,n}` being selected.
  :param n7: times of :math:`I_{2,n}` being selected.
  :param n8: times of :math:`I_{n}^{*}` being selected.
  
  :return: poly with key (i,i',j,l,l') and value :math:`b_{ii'jll'}`.
  :rtype: dict
  '''
  c = c_n(l, n1, n2, n3, n4, n5, n6, n7, n8)
  b = b_n(n1, n2, n3, n4, n5, n6, n7)
  if b.startswith('1'):
    coef = f'{c}{b[1:]}'   # b = '1...'
  else:
    coef = f'(-{c}){b[4:]}'# b = '(-1)...'
  #
  poly = moment_eIIeIII(n4, n5, n6, n7, n8)
  #
  n9 = l-n1-n2-n3-n4-n5-n6-n7-n8
  eJ = mcpp(n9)
  coef = f'({coef}) * ({eJ})'
  # 
  # multiply e^{(ik1+i'k2)t} t^j v_{1,n-1}^l v_{2,n-1}^{l'}
  #       by v_{1,n-1}^{n2} v_{2,n-1}^{n3}
  #          e^{-(n4k1 + n6k2)t} where t = nh
  poly_new = {(k[0]-n4, k[1]-n6, k[2], k[3]+n2, k[4]+n3): f'{coef} * (poly[k])'
              for k in poly.keys()} # k: (i,i',j,l,l')
  return(poly_new)

def decode_to_str(poly):
  expr = ''
  for key in poly.keys():
    i1, i2, j, l1, l2 = key
    t1 = f'exp(({i1}*k1+{i2}*k2)*n*h)'
    t2 = times_pow('n*h', j)
    t3 = moment_v(l1, 1)
    t4 = moment_v(l2, 2)
    expr += f' + ({poly[key]}) * {t1}{t2} * {t3} *{t4}'
  return(expr[3:]) # strip the leading ' + '

def moment_y_central(l):
  poly_sum = {}
  for n1 in range(l, -1, -1):
    cum = n1
    for n2 in range(l-cum, -1, -1):
      cum += n2
      for n3 in range(l-cum, -1, -1):
        cum += n3
        for n4 in range(l-cum, -1, -1):
          cum += n4
          for n5 in range(l-cum, -1, -1):
            cum += n5
            for n6 in range(l-cum, -1, -1):
              cum += n6
              for n7 in range(l-cum, -1, -1):
                cum += n7
                for n8 in range(l-cum, -1, -1):
                  n9 = l-cum-n8
                  poly_sum = merge(poly_sum, 
                                   moment_comb(l,n1,n2,n3,n4,n5,n6,n7,n8,n9))
  expr = decode_to_str(poly_sum) # a str expression of the moment
  return(expr)

# 

def b2_n(n1, n2, n3, n4, n5, n6, n7):
  b = '1' if (n2+n3+n5+n7)%2 == 0 else '(-1)'
  h1 = '(1-exp(-k1*h))/k1 - h'
  h2 = '(1-exp(-k2*h))/k2 - h'
  constant = f'mu*h + (({h1})*theta1 + ({h2})*theta2)/2'
  b += times_pow(f'{constant}', n1)
  b += times_pow(f'{h1}/2', n2)
  b += times_pow(f'{h2}/2', n3)
  b += times_pow('sigma_v1/(2*k1)', n4+n5)
  b += times_pow('sigma_v2/(2*k2)', n6+n7)
  return(b)

def moment_comb2(l, n1, n2, n3, n4, n5, n6, n7, n8):
  c = c_n(l, n1, n2, n3, n4, n5, n6, n7, n8)
  b = b2_n(n1, n2, n3, n4, n5, n6, n7)
  if b.startswith('1'):
    coef = f'{c}{b[1:]}'   # b = '1...'
  else:
    coef = f'(-{c}){b[4:]}'# b = '(-1)...'
  #
  poly = moment_eIIeIII(n4, n5, n6, n7, n8)
  #
  n9 = l-n1-n2-n3-n4-n5-n6-n7-n8
  eJ = mcpp(n9)
  coef = f'({coef}) * ({eJ})'
  # 
  # multiply e^{(ik1+i'k2)t} t^j v_{1,n-1}^l v_{2,n-1}^{l'}
  #       by v_{1,n-1}^{n2} v_{2,n-1}^{n3}
  #          e^{-(n4k1 + n6k2)t} where t = nh
  poly_new = {(k[0]-n4, k[1]-n6, k[2], k[3]+n2, k[4]+n3): f'{coef} * (poly[k])'
              for k in poly.keys()} # k: (i,i',j,l,l')
  return(poly_new)

def moment_y(l):
  poly_sum = {}
  for n1 in range(l, -1, -1):
    cum = n1
    for n2 in range(l-cum, -1, -1):
      cum += n2
      for n3 in range(l-cum, -1, -1):
        cum += n3
        for n4 in range(l-cum, -1, -1):
          cum += n4
          for n5 in range(l-cum, -1, -1):
            cum += n5
            for n6 in range(l-cum, -1, -1):
              cum += n6
              for n7 in range(l-cum, -1, -1):
                cum += n7
                for n8 in range(l-cum, -1, -1):
                  n9 = l-cum-n8
                  poly_sum = merge(poly_sum,
                                   moment_comb2(l,n1,n2,n3,n4,n5,n6,n7,n8,n9))
  expr = decode_to_str(poly_sum) # a str expression of the moment
  return(expr)

# 
def cb_vn(l, n1, n2, i):
  n3 = l - n1 -n2
  cv = math.comb(l, n1) * math.comb(l-n1, n2)
  coef  = f'{cv}'
  coef += times_exp(-n1, f'k{i}*h')
  coef += times_pow(f'(1-exp(-k{i}*h))*theta{i}', n2)
  coef += times_pow(f'sigma_v{i}', n3) + times_exp(-n3, f'k{i}*n*h')
  return(coef)

def g(n2,n3,n4,n5,n6,n7,n8,n9):
  poly = from_nh_to_np1h(moment_eIIeIII(n4, n5, n6, n7, n8))
  # 
  # multiply by outsiders
  # 
  eJ = mcpp(n9)
  plyn = {(k[0]-n4, k[1]-n6, k[2], k[3]+n2, k[4]+n3):
          f'{poly[k]} * ({eJ})' for k in poly.keys()}
  # 
  poly_eIeIvv = {} # eI_{1,n}^{i1} eI_{2,n}^{i2} v_{1,n-1}^{j1} v_{2,n-1}^{j2}
  #
  for key in plyn.keys():
    i1,i2,j,l1,l2 = key
    # decode e^{(i1k1 + i2k2)t} t^j, t=(n+1)*h back to coef
    pre_coef  = f'({plyn[key]}) * ({eJ})'
    pre_coef += f' * exp(({i1}*k1+{i2}*k2)*(n+1)*h)'
    pre_coef += times_pow('(n+1)*h', j)
    # expand v_{1,n}^{l1} v_{2,n}^{l2}
    for m1 in range(l1, -1, -1):
      for m2 in range(l1-m1, -1, -1):
        m3 = l1 - m1 - m2
        coef = cb_vn(l1, m1, m2, 1)
        coef = f'({pre_coef}) * {coef}'
        for m4 in range(l2, -1, -1):
          for m5 in range(l2-m4, -1, -1):
            m6 = l2 - m4 - m5
            coef = f'{coef} * {cn_vn(l2,m4,m5,2)}'
            # poly_eIeIvv
            k = (m3,m6,m1,m4) # To-check
            if k in poly_eIeIvv:
              poly_eIeIvv[k] += f' + {coef}'
            else:
              poly_eIeIvv[k]  = coef
  return(poly_eIeIvv)

def moment_inner_comb2(l1, m1, m2, m3, m4, m5, m6, m7, m8, poly_eIeIvv):
  c = c_n(l1, m1, m2, m3, m4, m5, m6, m7, m8)
  b = b2_n(m1, m2, m3, m4, m5, m6, m7)
  coef = f'{c}{b[1:]}' if b.startswith('1') else f'(-{c}){b[4:]}'
  # 
  poly_sum = {}
  for key in poly_eIeIvv.keys():
    i1,i2,j1,j2 = key
    b_iijj = poly_eIeIvv[key]
    poly = moment_eIIeIII(m4+i1, m5, m6+i2, m7, m8)
    plyn = {(k[0]-m4, k[1]-m6, k[2], k[3]+m2+j1, k[4]+m3+j2):
            f'({b_iijj}) * ({poly[k]})' for k in poly.keys()}
    poly_sum = merge(poly_sum, plyn)
  # 
  poly_sum = coef_poly(coef, poly_sum)
  return(poly_sum)

def moment_outer_comb2(l2, n1, n2, n3, n4, n5, n6, n7, n8, l1):
  n9 = l2 - n1-n2-n3-n4-n5-n6-n7-n8
  c = c_n(l2, n1, n2, n3, n4, n5, n6, n7, n8)
  b = b2_n(n1, n2, n3, n4, n5, n6, n7)
  coef = f'{c}{b[1:]}' if b.startswith('1') else f'(-{c}){b[4:]}'
  # 
  poly_eIeIvv = g(n2,n3,n4,n5,n6,n7,n8,n9)
  #
  poly_sum = {}
  for m1 in range(l1, -1, -1):
    cum = m1
    for m2 in range(l-cum, -1, -1):
      cum += m2
      for m3 in range(l-cum, -1, -1):
        cum += m3
        for m4 in range(l-cum, -1, -1):
          cum += m4
          for m5 in range(l-cum, -1, -1):
            cum += m5
            for m6 in range(l-cum, -1, -1):
              cum += m6
              for m7 in range(l-cum, -1, -1):
                cum += m7
                for m8 in range(l-cum, -1, -1):
                  m9 = l-cum-m8
                  poly = moment_inner_comb2(l1,m1,m2,m3,m4,m5,m6,m7,m8,
                                            poly_eIeIvv)
                  poly_sum = merge(poly_sum, poly)
  # 
  poly_sum = coef_poly(coef, poly_sum)
  return(poly_sum)

def moment_yy(l1, l2):
  poly = {}
  for n1 in range(l, -1, -1):
    cum = n1
    for n2 in range(l-cum, -1, -1):
      cum += n2
      for n3 in range(l-cum, -1, -1):
        cum += n3
        for n4 in range(l-cum, -1, -1):
          cum += n4
          for n5 in range(l-cum, -1, -1):
            cum += n5
            for n6 in range(l-cum, -1, -1):
              cum += n6
              for n7 in range(l-cum, -1, -1):
                cum += n7
                for n8 in range(l-cum, -1, -1):
                  n9 = l-cum-n8
                  poly = merge(poly, 
                               moment_outer_comb2(l2,n1,n2,n3,n4,n5,n6,n7,n8))
  expr = decode_to_str(poly)
  return(expr)

def cov_yy(l1, l2):
  m_yy = moment_yy(l1,l2)
  m_y_l1 = moment_y(l1)
  m_y_l2 = moment_y(l2)
  cov    = f'{m_yy} - ({m_y_l1}) * (m_y_l2)'
  return(cov)

'''
Module ito_moments

Functions to compute Ito process moments.
'''

from simplify_terms import times_pow, times_exp

def c_nmi(n,m,i):
  '''Coefficent :math:`c_{nmi}` as in :eq:`c-nmi`.
  
  :param n: n in :math:`e^{nkt}`.
  :type integer: int
  :param m: m in :math:`t^{m-i}`.
  :type integer: int
  :param i: i in :math:`t^{m-i}`.
  :type integer: int
  
  :return: the coefficient :math:`c_{nmi}`.
  :rtype: string
  '''
  prod = 1
  for j in range(m-i+1, m+1):
    prod = prod * j
  den = f'({n}*k)^({i+1})' if i+1 != 1 else f'({n}*k)'
  c = f'{(-1)**i * prod}/{den}'
  return(c)

def int_et(n,m):
  ''':math:`\int_{(n-1)h}^t e^{nks}s^mds=\sum_{ij}c_{ij}e^{ikt}t^j`
  
  :param n: n in :math:`e^{nks}s^m`.
  :type integer: int
  :param m: m in :math:`e^{nks}s^m`.
  :type integer: int
  
  :return: a dictionary with key-value pairs (i,j)-c_ij, 
   representing :math:`\sum_{ij}c_{ij}e^{ikt}t^j`.
  :rtype: dict of tuple: str
  '''
  if m < 0:
    msg = f"m in int_et(n,m) equals {m}, however it shouldn't be negative!"
    raise ValueError(msg)
  poly = {}
  if n == 0 and m == 0:
    poly[(0,1)] = '1'
    poly[(0,0)] = '-(n-1)*h'
  elif n == 0 and m != 0:
    poly[(0,m+1)] = f'1/({m+1})'
    poly[(0,0)] = f'-((n-1)*h)^({m+1})/({m+1})'
  elif n != 0 and m == 0:
    num = f'exp({n}*k*(n-1)*h)' if n != 1 else 'exp(k*(n-1)*h)'
    den = f'({n}*k)' if n != 1 else 'k'
    poly[(n,0)] = f'1/{den}'
    poly[(0,0)] = f' - {num}/{den}'
  else:
    den = f'({n}*k)' if n != 1 else 'k'
    num = f'exp({n}*k*(n-1)*h)' if n != 1 else 'exp(k*(n-1)*h)'
    num1   = num + times_pow('(n-1)*h', m)
    poly[(n,m)] = f'1/{den}'
    const    = f' - {num1}/{den}'
    for i in range(1, m+1):
      num2 = num + times_pow('(n-1)*h', m-i)
      poly[(n,m-i)] = c_nmi(n,m,i)
      const += f' - {num2}*({poly[(n,m-i)]})' 
    poly[(0,0)] = const
  return(poly)

def merge(poly1, poly2):
  '''Merge two polys by insert new ones and increment existing ones
  
  :param poly1: one polynomial of :math:`\sum_{i,j,l} b_{ijl} e^{ikt}t^j v_{n-1}^l`.
  :type poly: dict with key (i,j,l) and value :math:`b_{ijl}`
  :param poly2: another polynomial.
  :type poly: dict with key (i,j,l) and value :math:`b_{ijl}`
  
  :return: new polynomial.
  :rtype: dict with key (i,j,l) and value :math:`b_{ijl}`
  '''
  for ij in poly2.keys():
    if ij in poly1:
      poly1[ij] = f'{poly1[ij]} + {poly2[ij]}'
    else:
      poly1[ij] = poly2[ij]
  return(poly1)

def int_meII(m, n3, n4, eII):
  ''':math:`\int_{(n-1)h}^t e^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}|v_{n-1}]ds= \sum_{i,j,l} b_{ijl} e^{ikt}t^j v_{n-1}^l`
  
  :param m: m in :math:`e^{mkt}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}|v_{n-1}]`.
  :type integer: int
  :param n3: :math:`n_3` 
   in :math:`e^{mkt}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}|v_{n-1}]`
  :type integer: int
  :param n4: :math:`n_4` 
   in :math:`e^{mkt}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}|v_{n-1}]`
  :type integer: int
  :param eII: a dict with key (n3,n4) and value poly with key (i,j,l) and value :math:`b_{ijl}`.
  :type eII: dict
  
  :return: the integral result in poly with key (i,j,l) and value :math:`b_{ijl}`.
  :rtype: dict
  '''
  # get the representation of E[eI^{n3}I^{n4}]
  b = eII[(n3, n4)] # b: {(i,j,l): coef,...}
  #
  poly_sum = {}
  for ijl in b.keys():
    i, j, l = ijl
    c = int_et(i+m, j) # {(i,j): coef,...}
    for iijj in c.keys():
      coef = f'({b[ijl]}) * ({c[iijj]})'
      key = (iijj[0], iijj[1], l)
      if key in poly_sum:
        poly_sum[key] = f'{poly_sum[key]} + {coef}'
      else:
        poly_sum[key] = coef
  return(poly_sum)

def int_meIII(m, n3, n4, n5, eIII):
  ''':math:`\int_{(n-1)h}^te^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}I_{n-1,s}^{*n_5}|v_{n-1}]ds= \sum_{i,j,l} b_{ijl} e^{ikt}t^j v_{n-1}^l`
  
  :param m: m in :math:`e^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}I_{n-1,s}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param n3: :math:`n_3` 
   in :math:`e^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}I_{n-1,s}^{*n_5}|v_{n-1}]`
  :type integer: int
  :param n4: :math:`n_4` 
   in :math:`e^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}I_{n-1,s}^{*n_5}|v_{n-1}]`
  :type integer: int
  :param n5: :math:`n_5` 
   in :math:`e^{mks}E[eI_{n-1,s}^{n_3}I_{n-1,s}^{n_4}I_{n-1,s}^{*n_5}|v_{n-1}]`
  :type integer: int
  :param eIII: a dict with key (n3,n4,n5) and value poly with key (i,j,l) and value :math:`b_{ijl}`.
  :type eIII: dict
  
  :return: the integral result in poly with key (i,j,l) and value :math:`b_{ijl}`.
  :rtype: dict
  '''
  # get the representation of E[eI^{n3}I^{n4}I^{*n5}]
  b = eIII[(n3, n4, n5)] # b: {(i,j,l): coef,...}
  #
  poly_sum = {}
  for ijl in b.keys():
    i, j, l = ijl
    c = int_et(i+m, j) # {(i,j): coef,...}
    for iijj in c.keys():
      coef = f'({b[ijl]}) * ({c[iijj]})'
      key = (iijj[0], iijj[1], l)
      if key in poly_sum:
        poly_sum[key] = f'{poly_sum[key]} + {coef}'
      else:
        poly_sum[key] = coef
  return(poly_sum)

def coef_poly(coef, poly):
  '''Multiply poly with coefficient
  
  :param coef: coefficient to multiply by poly.
  :type string: str
  :param poly: one polynomial of :math:`\sum_{i,j,l} b_{ijl} e^{ikt}t^j v_{n-1}^l`.
  :type poly: dict with key (i,j,l) and value :math:`b_{ijl}`
  
  :return: new polynomial.
  :rtype: dict with key (i,j,l) and value :math:`b_{ijl}`
  '''
  for key in poly.keys():
    poly[key] = f'{coef} * ({poly[key]})'
  return(poly)

def coek_poly(coef, poly):
  '''Multiply poly with coefficient and update key
  
  :param coef: coefficient to multiply by poly.
  :type string: str
  :param poly: one polynomial of :math:`\sum_{i,j,l} b_{ijl} e^{ikt}t^j v_{n-1}^l`.
  :type poly: dict with key (i,j,l) and value :math:`b_{ijl}`
  
  :return: new polynomial.
  :rtype: dict with key (i,j,l) and value :math:`b_{ijl}`
  '''
  poly_new = {}
  for key in poly.keys():
    i, j, l = key
    poly_new[(i,j,l+1)] = f'{coef} * ({poly[key]})'
  return(poly_new)

def get_coef(prefix):
  '''Get the coefficients in Recursive Equations :eq:`ito-moment-i` and :eq:`ito-moment-ii`
  
  :param prefix: the leading term, such as :math:`\\frac{n_3(n_3-1)}{2}`.
  :type string: str
  
  :return: a tuple of four coefficients.
  :rtype: a tuple of four str
  '''
  c1 = prefix + 'exp(k*(n-1)*h)'
  c2 = '-' + prefix + 'theta*exp(k*(n-1)*h)'
  c3 = prefix + 'theta'
  c4 = prefix + 'sigma_v'
  return((c1,c2,c3,c4))

def recursive_eII(n3, n4, eII):
  '''Recursive step of :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]` as in equation :eq:`ito-moment-i`
  
  :param n3: :math:`n_3` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :type integer: int
  :param n4: :math:`n_4` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :type integer: int
  :param eII: a dict with key (n3,n4) and value poly with key (i,j,l) and value :math:`b_{ijl}`.
  :type eII: dict
  
  :return: updated eII.
  :rtype: dict
  '''
  poly_sum = {}
  #
  if n3 >= 2 and n4 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n3*(n3-1)/2)}*' if n3 > 2 else '')
    # (0,0,1): c1
    poly_sum = merge(coek_poly(c1, int_meII(1, n3-2, n4, eII)), poly_sum)
    poly_sum = merge(coef_poly(c2, int_meII(1, n3-2, n4, eII)), poly_sum)
    poly_sum = merge(coef_poly(c3, int_meII(2, n3-2, n4, eII)), poly_sum)
    poly_sum = merge(coef_poly(c4, int_meII(1, n3-1, n4, eII)), poly_sum)
  #
  if n3 >= 0 and n4 >= 2:
    c1,c2,c3,c4 = get_coef(f'{int(n4*(n4-1)/2)}*' if n4 > 2 else '')
    # (0,0,1): c1
    poly_sum = merge(coek_poly(c1, int_meII(-1, n3,   n4-2, eII)), poly_sum)
    poly_sum = merge(coef_poly(c2, int_meII(-1, n3,   n4-2, eII)), poly_sum)
    poly_sum = merge(coef_poly(c3, int_meII( 0, n3,   n4-2, eII)), poly_sum)
    poly_sum = merge(coef_poly(c4, int_meII(-1, n3+1, n4-2, eII)), poly_sum)
  #
  if n3 >= 1 and n4 >= 1:
    c1,c2,c3,c4 = get_coef(f'{n3*n4}*' if n3*n4 != 1 else '')
    # (0,0,1): coef1
    poly_sum = merge(coek_poly(c1, int_meII(0, n3-1, n4-1, eII)), poly_sum)
    poly_sum = merge(coef_poly(c2, int_meII(0, n3-1, n4-1, eII)), poly_sum)
    poly_sum = merge(coef_poly(c3, int_meII(1, n3-1, n4-1, eII)), poly_sum)
    poly_sum = merge(coef_poly(c4, int_meII(0, n3,   n4-1, eII)), poly_sum)
  return(poly_sum)

def recursive_eIII(n3, n4, n5, eIII):
  '''Recursive step of :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]` as in equation :eq:`ito-moment-ii`
  
  :param n3: :math:`n_3` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param n4: :math:`n_4` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param n5: :math:`n_5` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param eIII: a dict with key (n3,n4,n5) and value poly with key (i,j,l) and value :math:`b_{ijl}`.
  :type eIII: dict
  
  :return: updated eIII.
  :rtype: dict
  '''
  poly_sum = {}
  #
  if n3 >= 2 and n4 >=0 and n5 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n3*(n3-1)/2)}*' if n3 > 2 else '')
    poly_sum = merge(coev_poly(c1, int_meIII(1, n3-2, n4, n5, eIII)), poly_sum)
    poly_sum = merge(coef_poly(c2, int_meIII(1, n3-2, n4, n5, eIII)), poly_sum)
    poly_sum = merge(coef_poly(c3, int_meIII(2, n3-2, n4, n5, eIII)), poly_sum)
    poly_sum = merge(coef_poly(c4, int_meIII(1, n3-1, n4, n5, eIII)), poly_sum)
  #
  if n3 >= 0 and n4 >= 2 and n5 >= 0:
    c1,c2,c3,c4 = get_coef(f'{int(n4*(n4-1)/2)}*' if n4 > 2 else '')
    poly_sum = merge(coev_poly(c1, int_meIII(-1, n3, n4-2, n5, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c2, int_meIII(-1, n3, n4-2, n5, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c3, int_meIII( 0, n3, n4-2, n5, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c4, int_meIII(-1,n3+1,n4-2, n5, eIII)),poly_sum)
  #
  if n3 >= 1 and n4 >= 1 and n5 >= 0:
    c1,c2,c3,c4 = get_coef(f'{n3*n4}*' if n3*n4 != 1 else '')
    poly_sum = merge(coev_poly(c1, int_meIII(0, n3-1, n4-1, n5, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c2, int_meIII(0, n3-1, n4-1, n5, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c3, int_meIII(1, n3-1, n4-1, n5, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c4, int_meIII(0, n3,   n4-1, n5, eIII)),poly_sum)
  #
  if n3 >= 0 and n4 >= 0 and n5 >= 2:
    c1,c2,c3,c4 = get_coef(f'{int(n5*(n5-1)/2)} * ' if n5 > 2 else '')
    poly_sum = merge(coev_poly(c1, int_meIII(-1, n3, n4, n5-2, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c2, int_meIII(-1, n3, n4, n5-2, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c3, int_meIII( 0, n3, n4, n5-2, eIII)),poly_sum)
    poly_sum = merge(coef_poly(c4, int_meIII(-1,n3+1,n4, n5-2, eIII)),poly_sum)
  #
  return(poly_sum)

def moment_eII(n3, n4, return_all = False):
  ''':math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`
  
  :param n3: :math:`n_3` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :type integer: int
  :param n4: :math:`n_4` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}|v_{n-1}]`.
  :type integer: int
  :param return_pre: whether or not return lower order moments simultaneously,
     default to ``False``.
  :type bool: bool
  
  :return: poly if return_all=False else eII
  :rtype: dict or dict of dict
  '''
  # eII: a dict of moments of E[eI^{n3}I^{n4}]
  eII = {}
  # n3 + n4 = 0
  eII[(0,0)] = {(0,0,0): '1'} # support for special case
  # n3 + n4 = 1
  eII[(1,0)] = {}
  eII[(0,1)] = {}
  # n3 + n4 = 2
  #   each key-value: e^{ikt} t^j v_{n-1}^m - b_{ijm}
  eII[(2,0)] = {(2,0,0): 'theta/(2*k)',
                (1,0,1): 'exp(k*(n-1)*h)/k',
                (1,0,0): '-theta*exp(k*(n-1)*h)/k',
                (0,0,1): '-exp(2*k*(n-1)*h)/k',
                (0,0,0): 'theta*exp(2*k*(n-1)*h)/(2*k)'}

  eII[(1,1)] = {(1,0,0): 'theta/k',
                (0,1,1): 'exp(k*(n-1)*h)',
                (0,1,0): '-theta*exp(k*(n-1)*h)',
                (0,0,1): '-(n-1)*h*exp(k*(n-1)*h)',
                (0,0,0): '((n-1)*h - 1/k)*theta*exp(k*(n-1)*h)'}

  eII[(0,2)] = {(0,1,0): 'theta',
                (0,0,1): '1/k',
                (0,0,0): '-((n-1)*h + 1/k)*theta',
                (-1,0,1):'-exp(k*(n-1)*h)/k',
                (-1,0,0):'theta*exp(k*(n-1)*h)/k'}
  #
  if n3 + n4 <= 2:
    return( eII if return_all else eII[(n3,n4)] )
  #
  if n3 + n4 > 3:
    # compute all lower-order moments to get ready for the last computation
    for n in range(3, n3+n4):
      for i in range(n, -1, -1):
        poly = recursive_eII(i, n-i, eII)
        eII[(i,n-i)] = poly
  #
  poly = recursive_eII(n3, n4, eII)
  eII[(n3,n4)] = poly
  return( eII if return_all else poly )

def moment_eIII(n3, n4, n5, return_all = False):
  ''':math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`
  
  :param n3: :math:`n_3` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param n4: :math:`n_4` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param n4: :math:`n_5` in :math:`E[eI_{n-1,t}^{n_3}I_{n-1,t}^{n_4}I_{n-1,t}^{*n_5}|v_{n-1}]`.
  :type integer: int
  :param return_all: whether or not return lower order moments simultaneously,
     default to ``False``.
  :type bool: bool
  
  :return: poly if return_all=False else eIII
  :rtype: dict or dict of dict
  '''
  # eIII: a dict of moments of E[eI^{n3}I^{n4}I^{*n5}]
  if n3 + n4 + n5 < 0:
    raise ValueError(f"moment_eIII({n3},{n4},{n5}) is called!")
  eIII = {}
  # n3 + n4 + n5 = 0
  eIII[(0,0,0)] = {(0,0,0): '1'} # support for special case
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
    return( eIII if return_all else eIII[(n3,n4,n5)] )
  #
  if n3 + n4 + n5 > 3:
    # compute all lower-order moments to get ready for the last computation
    for n in range(3, n3+n4+n5):
      for i in range(n, -1, -1):
        for j in range(n-i, -1, -1):
          poly = recursive_eIII(i, j, n-i-j, eIII)
          eIII[(i,j,n-i-j)] = poly
  #
  poly = recursive_eIII(n3, n4, n5, eIII)
  eIII[(n3,n4,n5)] = poly
  return( eIII if return_all else poly )

def from_nh_to_np1h(poly):
  '''Change :math:`E[eI_{n}^{n_3} I_{n}^{n_4} I_{n}^{*n_5}|v_{n-1}]` to :math:`E[eI_{n+1}^{n_3} I_{n+1}^{n_4} I_{n+1}^{*n_5}|v_n]`
  
  :param poly: a poly of :math:`\sum{i,j,l} b_{ijl} e^{ikt} t^j v_{n-1}^l`, 
   where *t* = 'n*h'.
  :type poly: dict
  
  :return: a poly of :math:`\sum{i,j,l} b_{ijl} e^{ikt} t^j v_n^l`
   where *t* = '(n+1)*h', '(n-1)*h' in :math:`b_{ijl}` changed to 'n*h'.
  :rtype: dict
  '''
  for key in poly.keys():
    poly[key].replace('(n-1)*h', 'n*h')
  return(poly)


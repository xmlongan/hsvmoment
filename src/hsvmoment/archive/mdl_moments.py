'''
Module mdl_moments

Functions to compute Heston :abbr:`SV(Stochastic Volatility)` model moments.
'''
from simplify_terms import times_pow, times_exp
from ito_moments import moment_eIII, coef_poly, from_nh_to_np1h
from mdl_cmoments import c_n, moment_v, cb_vn, g, decode_to_str

def b2_n(n1,n2,n3,n4,n5):
  ''':math:`b_2(\\boldsymbol{n})` in :eq:`b2-n`
  
  :param n1: times of the constant in :math:`y_n` being selected.
  :type integer: int
  :param n2: times of :math:`v_{n-1}` being selected.
  :type integer: int
  :param n3: times of :math:`eI_n` being selected.
  :type integer: int
  :param n4: times of :math:`I_n` being selected.
  :type integer: int
  :param n5: times of :math:`I_n^{*}` being selected.
  :type integer: int
  
  :return: coefficient of :math:`b_2(\\boldsymbol{n})`.
  :rtype: str
  '''
  b  = f'1' if n2%2 == 0 else '(-1)'
  b += times_pow('(mu-theta/2)*h + theta*(1-exp(-k*h))/(2*k)', n1)
  b += times_pow('(1-exp(-k*h))/(2*k)', n2)
  b += times_pow('sigma_v/(2*k)', n3)
  b += times_pow('rho - sigma_v/(2*k)', n4)
  b += times_pow('1-rho^2', n5/2)
  return(b)

def moment_comb2(n, n1, n2, n3, n4, n5):
  '''Moment for this combination in expansion of :math:`y_n^l`
  
  :param n: l in :math:`E[\overline{y}_{n}^l]`.
  :type integer: int
  :param n1: times of the constant in :math:`y_n` being selected.
  :type integer: int
  :param n2: times of :math:`v_{n-1}` being selected.
  :type integer: int
  :param n3: times of :math:`eI_n` being selected.
  :type integer: int
  :param n4: times of :math:`I_n` being selected.
  :type integer: int
  :param n5: times of :math:`I_n^{*}` being selected.
  :type integer: int
  
  :return: poly with key (i,j,l) and value :math:`b_{ijl}`.
  :rtype: dict
  '''
  b = b2_n(n1, n2, n3, n4, n5)
  c = c_n(n, n1, n2, n3, n4, n5)
  if b.startswith('1'):
    coef = f'{c}{b[1:]}'    # b = '1...'
  else:
    coef = f'(-{c}){b[4:]}' # b = '(-1)...'
  #
  poly = moment_eIII(n3, n4, n5) # (0,0,0) already considered therein
  #
  # multiply e^{ikt} t^j v_{n-1}^l by v_{n-1}^{n2} e^{-n3kt} and coef
  #
  poly_new = {(k[0]-n3, k[1], k[2]+n2): f'{coef} * ({poly[k]})' 
              for k in poly.keys()} # (i,j,l) = k
  return(poly_new)

def moment_inner_comb2(l1, m1, m2, m3, m4, m5, poly_eIv):
  '''Moment for this combination in expansion of :math:`y_n^{l_1}` in :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`
  
  :param l1: *l1* in :math:`y_n^{l_1}`.
  :type integer: int
  :param m1: times of the constant in :math:`y_n` being selected.
  :param m2: times of :math:`v_{n-1}` being selected.
  :param m3: times of :math:`eI_{n}` being selected.
  :param m4: times of :math:`I_{n}` being selected.
  :param m5: times of :math:`I_{n}^{*}` being selected.
  :param poly_eIv: a poly dict of :math:`\sum_{i,j}b_{ij}eI_n^iv_{n-1}^j`.
  
  :return: a poly of :math:`\sum_{i,j,l} b_{ijl} e^{ikt} t^j v_{n-1}^l`
   where t = 'n*h'.
  :rtype: dict
  '''
  b = b2_n(m1,m2,m3,m4,m5)
  c = c_n(l1,m1,m2,m3,m4,m5)
  coef = f'{c}{b[1:]}' if b.startswith('1') else f'(-{c}){b[4:]}'
  #
  poly_sum = {}
  for key in poly_eIv.keys():
    i, j = key
    b_ij = poly_eIv[key]
    poly = moment_eIII(m3+i, m4, m5)
    plyn = {(k[0]-m3, k[1], k[2]+m2+j): f'({b_ij}) * ({poly[k]})'
            for k in poly.keys()}
    poly_sum = merge(poly_sum, plyn)
  #
  poly_sum = coef_poly(coef, poly_sum)
  return(poly_sum)

def moment_outer_comb2(l2, n1, n2, n3, n4, n5, l1):
  '''Moment for this combination in expansion of :math:`y_{n+1}^{l_2}` in :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`
  
  :param l2: *l2* in :math:`y_{n+1}^{l_2}`.
  :type integer: int
  :param n1: times of :math:`\\theta` being selected.
  :param n2: times of :math:`v_{n}` being selected.
  :param n3: times of :math:`eI_{n+1}` being selected.
  :param n4: times of :math:`I_{n+1}` being selected.
  :param n5: times of :math:`I_{n+1}^{*}` being selected.
  :param l1: *l1* in :math:`y_{n}^{l_1}`.
  
  :return: a poly of :math:`\sum_{i,j,l} b_{ijl} e^{ikt} t^j v_{n-1}^l`
   where t = 'n*h'.
  :rtype: dict
  '''
  b = b2_n(n1,n2,n3,n4,n5)
  c = c_n(l2,n1,n2,n3,n4,n5)
  coef = f'{c}{b[1:]}' if b.startswith('1') else f'(-{c}){b[4:]}'
  #
  poly_eIv = g(n2,n3,n4,n5)
  #
  poly_sum = {}
  #
  for m1 in range(l1, -1, -1):
    for m2 in range(l1-m1, -1, -1):
      for m3 in range(l1-m1-m2, -1, -1):
        for m4 in range(l1-m1-m2-m3, -1, -1):
          m5 = l1-m1-m2-m3-m4
          poly_sum = merge(poly_sum,
                           moment_inner_comb2(l1, m1, m2, m3, m4, m5, poly_eIv))
  #
  poly_sum = coef_poly(coef, poly_sum)
  return(poly_sum)

def moment_y(l):
  '''Moment of :math:`y_n` of order *l* as in equation :eq:`moment_y`
  
  :param l: *l* in :math:`E[y_{n}^l]`.
  :type integer: int
  
  :return: a string expression of the moment.
  :rtype: str
  '''
  poly_sum = {}
  for n1 in range(l, -1, -1):
    for n2 in range(l-n1, -1, -1):
      for n3 in range(l-n1-n2, -1, -1):
        for n4 in range(l-n1-n2-n3, -1, -1):
          n5 = l - n1 - n2 - n3 - n4
          poly_sum = merge(poly_sum, moment_comb2(l, n1, n2, n3, n4, n5))
  expr = decode_to_str(poly_sum) # a string expression of the moment
  return(expr)

def moment_yy(l1,l2):
  '''Moment :math:`E[y_n^{l_1}y_{n+1}^{l_2}]` as in equation :eq:`moment_yy`
  
  :param l1: *l1* in :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`.
  :param l2: *l2* in :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`.
  
  :return: a string expression of the moment.
  :rtype: str
  '''
  # cov(y_n^{l1}, y_{n+1}^{l2})
  poly = {} # container
  for n1 in range(l2, -1, -1):
    for n2 in range(l2-n1, -1, -1):
      for n3 in range(l2-n1-n2, -1, -1):
        for n4 in range(l2-n1-n2-n3, -1, -1):
          n5 = l2 - n1 - n2 - n3 - n4
          poly = merge(poly, moment_outer_comb2(l2, n1, n2, n3, n4, n5, l1))
  expr = decode_to_str(poly) # a string expression of the moment
  return(expr)

def cov_yy(l1,l2):
  '''Covariance :math:`cov(y_n^{l_1},y_{n+1}^{l_2})` in equation :eq:`cov_yy`
  
  :param l1: *l1* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  :param l2: *l2* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  
  :return: a string expression of the moment.
  :rtype: str
  '''
  m_yy = moment_yy(l1,l2)
  m_y_l1 = moment_y(l1)
  m_y_l2 = moment_y(l2)
  cov    = f'{m_yy} - ({m_y_l1}) * (m_y_l2)'
  return(cov)

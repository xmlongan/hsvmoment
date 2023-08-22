'''
Module mdl_cmoments

Functions to compute Heston :abbr:`SV(Stochastic Volatility)` model 
central moments.
'''

from ito_moments import moment_eIII, coef_poly, from_nh_to_np1h
from simplify_terms import times_pow, times_exp
import math

def c_n(l, n1, n2, n3, n4, n5):
  ''':math:`c(\\boldsymbol{n})` in :eq:`c-n`
  
  :param l: l in :math:`E[\overline{y}_{n}^l]`.
  :type integer: int
  :param n1: times of :math:`\\theta` being selected.
  :type integer: int
  :param n2: times of :math:`v_{n-1}` being selected.
  :type integer: int
  :param n3: times of :math:`eI_n` being selected.
  :type integer: int
  :param n4: times of :math:`I_n` being selected.
  :type integer: int
  :param n5: times of :math:`I_n^{*}` being selected.
  :type integer: int
  
  :return: number of this special combination.
  :rtype: int
  '''
  num  = math.comb(l, n1)
  num *= math.comb(l-n1, n2)
  num *= math.comb(l-n1-n2, n3)
  num *= math.comb(l-n1-n2-n3, n4)
  return(num)

def b_n(n1,n2,n3,n4,n5):
  ''':math:`b(\\boldsymbol{n})` in :eq:`b-n`
  
  :param n1: times of :math:`\\theta` being selected.
  :type integer: int
  :param n2: times of :math:`v_{n-1}` being selected.
  :type integer: int
  :param n3: times of :math:`eI_n` being selected.
  :type integer: int
  :param n4: times of :math:`I_n` being selected.
  :type integer: int
  :param n5: times of :math:`I_n^{*}` being selected.
  :type integer: int
  
  :return: coefficient of :math:`b(\\boldsymbol{n})`.
  :rtype: str
  '''
  b  = f'1' if n2%2 == 0 else '(-1)'
  b += times_pow('theta', n1)
  b += times_pow('(1-exp(-k*h))/(2*k)', n1+n2)
  b += times_pow('sigma_v/(2*k)', n3)
  b += times_pow('rho - sigma_v/(2*k)', n4)
  b += times_pow('1-rho^2', n5/2)
  return(b)

def moment_v(n):
  '''Moment of :math:`v_{n-1}` as in equation :eq:`moment-v`
  
  :param n: order of the moment.
  :type integer: int
  
  :return: n order moment.
  :rtype: str
  '''
  if n == 0:
    return('1')
  for j in range(n):
    if j == 0:
      moment = 'theta'
    else:
      moment = f'{moment} * (theta + {j}*sigma_v^2/(2*k))'
  return(moment)

def moment_comb(n, n1, n2, n3, n4, n5):
  '''Moment for this combination in expansion of :math:`\overline{y}_n^l`
  
  :param n: l in :math:`E[\overline{y}_{n}^l]`.
  :type integer: int
  :param n1: times of :math:`\\theta` being selected.
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
  b = b_n(n1, n2, n3, n4, n5)
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

def cb_vn(l, n1, n2):
  '''Coefficient :math:`c_v(\\boldsymbol{m})\cdot b_v(\\boldsymbol{m})` in equations :eq:`cv-m` and :eq:`bv-m` when expanding :math:`v_n^l`
  
  :param l: *l* in :math:`v_n^l` where :math:`v_n` is expanded as :eq:`vn-expanded`.
  :type integer: int
  :param n1: times of :math:`v_{n-1}` being selected.
  :type integer: int
  :param n2: times of :math:`\\theta` being selected.
  :type integer: int
  
  :return: coefficient of this combination.
  :rtype: str
  '''
  # v_n^l, expand 
  # v_n = exp(-kh)*v_{n-1} + (1-e^{-kh})*theta + sigma_v*e^{-knh}eI_n
  n3 = l - n1 - n2
  cv = math.comb(l, n1) * math.comb(l-n1, n2)
  coef  = f'{cv}'
  coef += times_exp(-n1, 'k*h')
  coef += times_pow('(1-exp(-k*h))*theta', n2)
  coef += times_pow('sigma_v', n3) + times_exp(-n3, 'k*n*h')
  return(coef)

def g(n2,n3,n4,n5):
  ''' Function :math:`g(n_2,n_3,n_4,n_5)` see equation :eq:`g`.
  
  :param n2: times of :math:`v_n` being selected.
  :type integer: int
  :param n3: times of :math:`eI_{n+1}` being selected.
  :type integer: int
  :param n4: times of :math:`I_{n+1}` being selected.
  :type integer: int
  :param n5: times of :math:`I_{n+1}^{*}` being selected.
  :type integer: int
  
  :return: a poly dict of :math:`\sum_{i,j}b_{ij}eI_n^iv_{n-1}^j`.
  :rtype: dict
  '''
  # get E[eI_{n+1}^{n_3} I_{n+1}^{n_4} I_{n+1}^{*n_5}|v_n]
  poly = from_nh_to_np1h(moment_eIII(n3, n4, n5))
  # \sum{i,j,l} b_{ijl} e^{ikt} t^j v_n^l, t = (n+1)h
  #
  # multiply by v_n^{n2} and e^{-n3*kt}, t = (n+1)h
  plyn = {(k[0]-n3, k[1], k[2]+n2): poly[k] for k in poly.keys()}
  # i, j, l = k
  #
  poly_eIv = {} # container for \sum_{i,j} b_{ij} eI_n^iv_{n-1}^j
  #
  for key in plyn.keys():
    i, j, l = key
    # decode the e^{ikt} t^j, t = (n+1)*h back to coef
    pre_coef  = f'({plyn[key]})'
    pre_coef += times_exp(i, 'k*(n+1)*h') 
    pre_coef += times_pow('(n+1)*h', j)
    #
    # expand v_n^l
    for m1 in range(l, -1, -1):
      for m2 in range(l-m1, -1, -1):
        m3 = l-m1-m2
        coef = cb_vn(l, m1, m2) # coefficient of this combination
        coef = f'({pre_coef}) * {coef}'
        # poly_eIv
        # k = (m1,m3) # -> (m3,m1)? !!! To-check
        k = (m3, m1) # eI is placed first, then v
        if k in poly_eIv:
          poly_eIv[k] += f' + {coef}'
        else:
          poly_eIv[k]  = coef
  return(poly_eIv)

def moment_inner_comb(l1, m1, m2, m3, m4, m5, poly_eIv):
  '''Moment for this combination in expansion of :math:`\overline{y}_n^{l_1}` in :math:`E[\overline{y}_n^{l_1}\overline{y}_{n+1}^{l_2}]`
  
  :param l1: *l1* in :math:`\overline{y}_n^{l_1}`.
  :type integer: int
  :param m1: times of :math:`\\theta` being selected.
  :param m2: times of :math:`v_{n-1}` being selected.
  :param m3: times of :math:`eI_{n}` being selected.
  :param m4: times of :math:`I_{n}` being selected.
  :param m5: times of :math:`I_{n}^{*}` being selected.
  :param poly_eIv: a poly dict of :math:`\sum_{i,j}b_{ij}eI_n^iv_{n-1}^j`.
  
  :return: a poly of :math:`\sum_{i,j,l} b_{ijl} e^{ikt} t^j v_{n-1}^l`
   where t = 'n*h'.
  :rtype: dict
  '''
  b = b_n(m1,m2,m3,m4,m5)
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

def moment_outer_comb(l2, n1, n2, n3, n4, n5, l1):
  '''Moment for this combination in expansion of :math:`\overline{y}_{n+1}^{l_2}` in :math:`E[\overline{y}_n^{l_1}\overline{y}_{n+1}^{l_2}]`
  
  :param l2: *l2* in :math:`\overline{y}_{n+1}^{l_2}`.
  :type integer: int
  :param n1: times of :math:`\\theta` being selected.
  :param n2: times of :math:`v_{n}` being selected.
  :param n3: times of :math:`eI_{n+1}` being selected.
  :param n4: times of :math:`I_{n+1}` being selected.
  :param n5: times of :math:`I_{n+1}^{*}` being selected.
  :param l1: *l1* in :math:`\overline{y}_{n}^{l_1}`.
  
  :return: a poly of :math:`\sum_{i,j,l} b_{ijl} e^{ikt} t^j v_{n-1}^l`
   where t = 'n*h'.
  :rtype: dict
  '''
  b = b_n(n1,n2,n3,n4,n5)
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
                           moment_inner_comb(l1, m1, m2, m3, m4, m5, poly_eIv))
  #
  poly_sum = coef_poly(coef, poly_sum)
  return(poly_sum)

def decode_to_str(poly):
  '''Decode the poly dict back to string
  
  :param poly: a dict in {(i,j,l): b_ijl,...} 
   representing :math:`\sum_{i,j,l}b_{ijl} e^{ikt} t^j v_{n-1}^l`.
  :type poly: dict
  
  :return: string expression of the poly, :math:`t` in :math:`e^{ikt} t^j`
   substituted with :math:`n*h`.
  :rtype: str
  '''
  expr = ''
  for key in poly.keys():
    i, j, l = key
    t1 = times_exp(i, 'k*n*h') # e^{ikt} (t=nh)
    t2 = times_pow('n*h', j)   # t^j     (t=nh)
    t3 = moment_v(l)           # E[v_{n-1}^l]
    expr += f' + ({poly[key]}){t1}{t2} * {t3}'
  return(expr[3:]) # exclude the leading ' + '

def moment_y_central(l):
  '''Central moment of :math:`y_n` of order *l* as in equation :eq:`moment_y_central`
  
  :param l: *l* in :math:`E[\overline{y}_{n}^l]`.
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
          poly_sum = merge(poly_sum, moment_comb(l, n1, n2, n3, n4, n5))
  expr = decode_to_str(poly_sum) # a string expression of the moment
  return(expr)

def moment_yy_central(l1,l2):
  '''Moment :math:`E[\overline{y}_n^{l_1}\overline{y}_{n+1}^{l_2}]` as in equation :eq:`moment_yy_central`
  
  :param l1: *l1* in :math:`E[\overline{y}_n^{l_1}\overline{y}_{n+1}^{l_2}]`.
  :param l2: *l2* in :math:`E[\overline{y}_n^{l_1}\overline{y}_{n+1}^{l_2}]`.
  
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
          poly = merge(poly, moment_outer_comb(l2, n1, n2, n3, n4, n5, l1))
  expr = decode_to_str(poly) # a string expression of the moment
  return(expr)

def cov_yy_central(l1,l2):
  '''Covariance :math:`cov(\overline{y}_n^{l_1},\overline{y}_{n+1}^{l_2})` as in equation :eq:`cov-yy-central`
  
  :param l1: *l1* in :math:`cov(\overline{y}_n^{l_1},\overline{y}_{n+1}^{l_2})`.
  :param l2: *l2* in :math:`cov(\overline{y}_n^{l_1},\overline{y}_{n+1}^{l_2})`.
  
  :return: a string expression of the moment.
  :rtype: str
  '''
  m_cycy = moment_yy_central(l1,l2)
  m_y_l1 = moment_y_central(l1)
  m_y_l2 = moment_y_central(l2)
  cov    = f'{m_cycy} - ({m_y_l1}) * (m_y_l2)'
  return(cov)
  

'''
:abbr:`CPP(Compound Poisson Process)` Moments

Implementing functions to compute moments of Compound Poisson Process

.. math::
   
   J_n = \sum_{N((n-1)h)+1}^{N(nh)}j_i, j_i \sim \mathcal{N}(\mu_j,\sigma_j^2)

'''
from simplify_terms import times_pow

def merge(poly1, poly2):
  for key in poly2.keys():
    if key not in poly1:
      poly1[key] = poly2[key]
    else:
      poly1[key] = poly1[key] + poly2[key]
  return(poly1)

def dmgf(poly):
  '''derivative of moment generating function
  
  :param poly: poly representation of MGF 
     as :math:`\sum_{i,j}b_{ij}(\mu+\sigma^2s)^i \sigma^{2j}`
  
  :return: derivative of MGF
  :return: dict
  '''
  # step 1
  poly1 = {(k[0]+1,k[1]): poly[k] for k in poly.keys()}
  #
  # step 2
  poly2 = {(k[0]-1,k[1]+1): poly[k]*k[0] for k in poly.keys() if k[0] > 0}
  #
  # merge the two poly
  return(merge(poly1, poly2))

def decode_to_str(poly):
  expr = ''
  for ij in poly.keys():
    i, j = ij
    if poly[ij] == 1:
      if i == 0 and j == 0:
        expr += f' + 1'
      elif i == 0:
        expr += f' + ' + times_pow('sigma_j',2*(j))[3:]
      else:
        expr += f' + ' + times_pow('mu_j',i)[3:] + times_pow('sigma_j',2*(j))
    else:
      expr += f' + {poly[ij]}' + times_pow('mu_j',i)+times_pow('sigma_j',2*(j))
  return(expr[3:])

def mnorm(n):
  '''Moment of normal distribution with order n
  
  :param n: order of the moment.
  
  :return: expression of the moment.
  :rtype: str
  '''
  if n < 0:
    raise ValueError(f'mnorm(n): Order of the moment is {n} which must > 0!')
  M = {}
  M[0] = {(0,0): 1}
  M[1] = {(1,0): 1}
  for i in range(2,n+1):
    M[i] = dmgf(M[i-1])
  #
  return(decode_to_str(M[n]))

def d1_times_key(key):
  '''Update the key after multiply :math:`(\lambda h)M'(s)`
  
  :param key: a tuple of (i,(n1,m1),...,(nl,ml)).
  
  :return: updated key.
  :rtype: list
  '''
  # multiply (lambda h)
  knw = list(key)
  knw[0] += 1
  # multiply M'(s)
  if key[1][0] == 1: # (n1,m1),(n2,m2),...,(nj,mj) where n1 < n2 < ... < nj
    knw[1] = (1, key[1][1] + 1)
  else:
    knw.insert(1, (1,1))
  return(tuple(knw))

def dterm(key, coef):
  poly = {}
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
    #
    knw = tuple(knw)
    if knw in poly:
      poly[knw] += coef_new
    else:
      poly[knw]  = coef_new
  #
  return(poly)

def dmgf_cpp(poly):
  '''derivative of Moment Generating Function of CPP
  
  :param poly: poly representation of MGF as 
     :math:`\sum_{i,(n_1,m_1),...,(n_l,m_l)} b (\lambda h)^{i}M^{(n_1)m_1}(s) ... M^{(n_l)m_l}(s)`
  
  :return: derivative of MGF
  :rtype: dict
  '''
  # step 1
  poly_sum = {}
  for key in poly.keys():
    knw = d1_times_key(key)
    poly_sum[knw] = poly[key]
  #
  # step 2
  for key in poly.keys():
    poly_sum = merge(dterm(key, poly[key]), poly_sum)
  #
  return(poly_sum)

def decode_to_str_cpp(poly):
  expr = ''
  ks = sorted(poly.keys(), key=lambda key: (key[1][0],-key[1][1]))
  for k in ks:
    i = k[0]
    txt  = f'{poly[k]}*' if poly[k] != 1 else ''
    txt += times_pow('lambda*h', i)[3:] # i > 0
    for j in range(1, len(k)):
      n_j, m_j = k[j]
      # txt += f'*({mnorm(n_j)})^({m_j})'
      txt += times_pow(mnorm(n_j), m_j) # m_j > 0
      # txt += f'*M^(({n_j}){m_j})'
    expr += f' + {txt}'
  return(expr[3:])

def mcpp(n):
  '''Moment of Compound Poisson Process with order n
  
  :param n: order of the moment.
  :rtype: str
  '''
  if n < 0:
    raise ValueError(f'mcpp(n): Order of the moment is {n} which must > 0!')
  M = {}
  M[0] = {(0): 1}
  M[1] = {(1,(1,1)): 1}
  # M[2] = {(2,(1,2)): 1,
  #         (1,(2,1)): 1}
  # M[3] = {(3,(1,3)): 1,
  #         (2,(1,1),(2,1)): 3,
  #         (1,(3,1)): 1}
  for i in range(2,n+1):
    M[i] = dmgf_cpp(M[i-1])
  #
  return(decode_to_str_cpp(M[n]))


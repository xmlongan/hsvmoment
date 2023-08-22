'''
Covariances for Two-Factor SV
'''
import math
from fractions import Fraction as Frac

import sys, os
file_path = os.path.abspath(__file__) 
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(file_path)))
if src_dir not in sys.path: sys.path.append(src_dir)
# so that hsvmoment dir will be treated as a package

from hsvmoment.poly import Poly
from hsvmoment.ito_mmnt import moment_v
from hsvmoment.itos_mmnt import moment_eIIeIII, t_mul_t0
from hsvmoment.mdl_1fsv.cov import exp_vn

from hsvmoment.mdl_2fsv.cmoment import (
  c_n, 
  vvee_eIIeIII
)
from hsvmoment.mdl_2fsv.moment import (
  b_n, 
  sub_v, 
  moment_y,
  dpoly,
  poly2num
)

def vvee_eIIeIII_vnvn(n2, n3, n4, n5, n6, n7, n8):
  '''Multiply :math:`eIIeIII` by :math:`v_{1,n}v_{2,n}e^{-k_1(n+1)h}e^{-k_2(n+1)h}`
  
  Compute
  :math:`E[v_{1,n}^{n_2}v_{2,n}^{n_3}e^{-n_4k_1(n+1)h}e^{-n_6k_2(n+1)h} 
  eI_{1,n+1}^{n_4} I_{1,n+1}^{n_5} eI_{2,n+1}^{n_6} I_{2,n+1}^{n_7} 
  I_{n+1}^{*n_8}|v_{1,n}v_{2,n}]` and expand :math:`v_{1,n}` and 
  :math:`v_{2,n}`.
  
  :param n2: times of :math:`v_{1,n}` being selected.
  :param n3: times of :math:`v_{2,n}` being selected.
  :param n4: times of :math:`eI_{1,n+1}` being selected.
  :param n5: times of :math:`I_{1,n+1}` being selected.
  :param n6: times of :math:`eI_{2,n+1}` being selected.
  :param n7: times of :math:`I_{2,n+1}` being selected.
  :param n8: times of :math:`I_{n+1}^{*}` being selected.
  
  :return: poly with attribute ``keyfor`` = 
     ('e^{-k1*nh}eI_{1,n}','e^{-k2*nh}eI_{2,n}',
     '(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h',
     'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2').
  :rtype: Poly
  '''
  poly = vvee_eIIeIII(n2, n3, n4, n5, n6, n7, n8)
  # ['(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
  #  'e^{-(n1*k1+n2*k2)h}','h',
  # 'v_{1,n}','theta1','sigma_v1', 'v_{2,n}','theta2','sigma_v2']
  # 
  # expand v1_n v2_n
  poly_eIv = Poly()
  kf = ['e^{-k1*nh}eI_{1,n}','e^{-k2*nh}eI_{2,n}',
    '(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
    'e^{-(n1*k1+n2*k2)h}','h',
    'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2']
  poly_eIv.set_keyfor(kf)
  # 
  poln = Poly()
  for k in poly:
    polv = exp_vn(k[3])
    # ['e^{-k1nh}eI_{1,n}', 'e^{-k1h}', 'v_{1,n-1}', 'theta1', 'sigma_v1']
    for k1 in polv:
      t1 = (k[1][0]+k1[1], k[1][1])
      key = (k1[0], k[0], t1, k[2], 
        k1[2], k[4]+k1[3], k[5]+k1[4],
        k[6], k[7], k[8])
      val = poly[k] * polv[k1]
      poln.add_keyval(key, val)
  # 
  for k in poln:
    polv = exp_vn(k[7])
    # ['e^{-k2nh}eI_{2,n}', 'e^{-k2h}', 'v_{2,n-1}', 'theta2', 'sigma_v2']
    for k1 in polv:
      t1 = (k[2][0], k[2][1]+k1[1])
      key = (k[0], k1[0], k[1], t1, k[3], 
        k[4], k[5], k[6],
        k1[2], k[8]+k1[3], k[9]+k1[4])
      val = poln[k] * polv[k1]
      poly_eIv.add_keyval(key, val)
  return(poly_eIv)

def moment_inner_comb(l1, m1, m2, m3, m4, m5, m6, m7, m8, poly_eIv):
  '''Moment for this inner combination in expansion of :math:`y_n^{l_1}`
  
  :param l1: *l1* in :math:`y_n^{l_1}`.
  :param m1: times of constant being selected.
  :param m2: times of :math:`v_{1,n-1}` being selected.
  :param m3: times of :math:`v_{2,n-1}` being selected.
  :param m4: times of :math:`eI_{1,n}` being selected.
  :param m5: times of :math:`I_{1,n}` being selected.
  :param m6: times of :math:`eI_{2,n}` being selected.
  :param m7: times of :math:`I_{2,n}` being selected.
  :param m8: times of :math:`I_{n}^{*}` being selected.
  :param poly_eIv: poly with attribute ``keyfor`` = 
     ('e^{-k1*nh}eI_{1,n}','e^{-k2*nh}eI_{2,n}',
     '(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h',
     'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2').
  
  :return: poly with attribute ``keyfor`` = 
     ('(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h','mu',
     'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2').
  :rtype: Poly
  '''
  poly = Poly()
  kf = ['(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
    'e^{-(n1*k1+n2*k2)h}','h','mu',
    'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2']
  poly.set_keyfor(kf)
  # 
  # combine and compute
  # k1: ['e^{-k1nh}eI_{1,n}','e^{-k2*nh}eI_{2,n}',
  # 
  # '(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
  # 'e^{-(n1*k1+n2*k2)h}','h', 'v_{1,n-1}','theta1','sigma_v1',
  # 'v_{2,n-1}','theta2','sigma_v2']
  # 
  # k2: ['(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
  # 'e^{-(n1*k1+n2*k2)h}','h', 'v_{1,n-1}','theta1','sigma_v1',
  # 'v_{2,n-1}','theta2','sigma_v2']
  for k1 in poly_eIv:
    poln = vvee_eIIeIII(m2+k1[5], m3+k1[8], m4+k1[0], m5, m6+k1[1], m7, m8)
    for k2 in poln:
      # t0 = k2[0] + k1[2]
      t0 = k1[2]
      for t in k2[0]: t0 = t_mul_t0(t, t0)
      # 
      t1 = (k1[3][0] + k2[1][0], k1[3][1] + k2[1][1])
      #              h,           mu,
      key = (t0, t1, k2[2]+k1[4], 0, 
      # v_{1,n-1},  theta1,    sigma_v1,
        k2[3], k2[4]+k1[6], k2[5]+k1[7],
      # v_{2,n-1},  theta2,    sigma_v2
        k2[6], k2[7]+k1[9], k2[8]+k1[10])
      val = poly_eIv[k1] * poln[k2]
      poly.add_keyval(key, val)
  # 
  b = b_n(m1, m2, m3, m4, m5, m6, m7)
  # ['e^{-k1h}','k1^{-}','theta1','sigma_v1',
  #  'e^{-k2h}','k2^{-}','theta2','sigma_v2', 'h','mu']
  poln = Poly()
  poln.set_keyfor(kf)
  # 
  for k1 in poly:
    for k2 in b:
      # t0 = ((1,0,k2[1]), (0,1,k2[5])) + k1[0]
      t0 = t_mul_t0((1,0,k2[1]), k1[0])
      t0 = t_mul_t0((0,1,k2[5]), t0)
      # 
      t1 = (k1[1][0]+k2[0], k1[1][1]+k2[4])
      #              h,           mu,
      key = (t0, t1, k1[2]+k2[8], k1[3]+k2[9],
      # v_{1,n-1},  theta1,    sigma_v1,
        k1[4], k1[5]+k2[2], k1[6]+k2[3],
      # v_{2,n-1},  theta2,    sigma_v2,
        k1[7], k1[8]+k2[6], k1[9]+k2[7])
      val = poly[k1] * b[k2]
      poln.add_keyval(key, val)
  # 
  c = c_n(l1, m1,m2,m3,m4,m5,m6,m7,m8)
  return(c * poln)

def moment_outer_comb(l2, n1, n2, n3, n4, n5, n6, n7, n8, l1):
  '''Moment for this outer combination in expansion of :math:`y_{n+1}^{l_2}`
  
  :param l2: *l2* in :math:`y_{n+1}^{l_2}`.
  :param n1: times of constant being selected.
  :param n2: times of :math:`v_{1,n}` being selected.
  :param n3: times of :math:`v_{2,n}` being selected.
  :param n4: times of :math:`eI_{1,n+1}` being selected.
  :param n5: times of :math:`I_{1,n+1}` being selected.
  :param n6: times of :math:`eI_{2,n+1}` being selected.
  :param n7: times of :math:`I_{2,n+1}` being selected.
  :param n8: times of :math:`I_{n+1}^{*}` being selected.
  :param l1: *l1* in :math:`y_{n}^{l_1}`.
  
  :return: poly with attribute ``keyfor`` = 
     ('(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h','mu',
     'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2').
  :rtype: Poly
  '''
  poly_eIv = vvee_eIIeIII_vnvn(n2, n3, n4, n5, n6, n7, n8)
  # 
  poly = Poly()
  kf = ['(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
    'e^{-(n1*k1+n2*k2)h}','h','mu',
    'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2']
  poly.set_keyfor(kf)
  # 
  for m1 in range(l1, -1, -1):
    for m2 in range(l1-m1, -1, -1):
      for m3 in range(l1-m1-m2, -1, -1):
        for m4 in range(l1-m1-m2-m3, -1, -1):
          for m5 in range(l1-m1-m2-m3-m4, -1, -1):
            for m6 in range(l1-m1-m2-m3-m4-m5, -1, -1):
              for m7 in range(l1-m1-m2-m3-m4-m5-m6, -1, -1):
                m8 = l1-m1-m2-m3-m4-m5-m6-m7
                poln = moment_inner_comb(l1,m1,m2,m3,m4,m5,m6,m7,m8,poly_eIv)
                poly.merge(poln)
  #
  b = b_n(n1, n2, n3, n4, n5, n6, n7)
  # ['e^{-k1h}','k1^{-}','theta1','sigma_v1',
  #  'e^{-k2h}','k2^{-}','theta2','sigma_v2', 'h','mu']
  poln = Poly()
  poln.set_keyfor(kf)
  # 
  for k1 in poly:
    for k2 in b:
      # t0 = ((1,0,k2[1]), (0,1,k2[5])) + k1[0]
      t0 = t_mul_t0((1,0,k2[1]), k1[0])
      t0 = t_mul_t0((0,1,k2[5]), t0)
      # 
      t1 = (k1[1][0]+k2[0], k1[1][1]+k2[4])
      #              h,           mu,
      key = (t0, t1, k1[2]+k2[8], k1[3]+k2[9],
      # v_{1,n-1},  theta1,    sigma_v1,
        k1[4], k1[5]+k2[2], k1[6]+k2[3],
      # v_{2,n-1},  theta2,    sigma_v2,
        k1[7], k1[8]+k2[6], k1[9]+k2[7])
      val = poly[k1] * b[k2]
      poln.add_keyval(key, val)
  # 
  c = c_n(l2, n1,n2,n3,n4,n5,n6,n7,n8)
  return(c * poln)

def moment_yy(l1, l2):
  '''Co-Moment :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`
  
  :param l1: *l1* in :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`.
  :param l2: *l2* in :math:`E[y_n^{l_1}y_{n+1}^{l_2}]`.
  
  :return: poly with attribute ``keyfor`` =  
     ('(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h','mu', 'theta1','sigma_v1','theta2','sigma_v2').
  :rtype: Poly
  '''
  poly = Poly()
  kf = ['(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
    'e^{-(n1*k1+n2*k2)h}','h','mu',
    'v_{1,n-1}','theta1','sigma_v1', 'v_{2,n-1}','theta2','sigma_v2']
  poly.set_keyfor(kf)
  # 
  for n1 in range(l2, -1, -1):
    for n2 in range(l2-n1, -1, -1):
      for n3 in range(l2-n1-n2, -1, -1):
        for n4 in range(l2-n1-n2-n3, -1, -1):
          for n5 in range(l2-n1-n2-n3-n4, -1, -1):
            for n6 in range(l2-n1-n2-n3-n4-n5, -1, -1):
              for n7 in range(l2-n1-n2-n3-n4-n5-n6, -1, -1):
                n8 = l2-n1-n2-n3-n4-n5-n6-n7
                poln = moment_outer_comb(l2,n1,n2,n3,n4,n5,n6,n7,n8,l1)
                poly.merge(poln)
  poly_sum = sub_v(poly)
  poly_sum.remove_zero()
  return(poly_sum)

def cov_yy(l1, l2):
  '''Covariance :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`
  
  :param l1: *l1* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  :param l2: *l2* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  
  :return: poly with attribute ``keyfor`` =
     ('(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h','mu', 'theta1','sigma_v1','theta2','sigma_v2').
  :rtype: Poly
  '''
  # cov = moment_yy(l1, l2) - (moment_y(l1) * moment_y(l2))
  poln = Poly()
  kf = ['(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
    'e^{-(n1*k1+n2*k2)h}','h','mu', 'theta1','sigma_v1','theta2','sigma_v2']
  poln.set_keyfor(kf)
  # 
  pol1 = moment_y(l1)
  pol2 = moment_y(l2)
  for k1 in pol1:
    for k2 in pol2:
      # t0 = k1[0] + k2[0]
      t0 = k1[0]
      for t in k2[0]: t0 = t_mul_t0(t, t0)
      # 
      t1 = (k1[1][0]+k2[1][0], k1[1][1]+k2[1][1])
      key = (t0, t1, k1[2]+k2[2], k1[3]+k2[3], 
        k1[4]+k2[4], k1[5]+k2[5], k1[6]+k2[6], k1[7]+k2[7])
      val = pol1[k1] * pol2[k2]
      poln.add_keyval(key, val)
  # 
  cov = moment_yy(l1, l2) - poln
  cov.remove_zero()
  return(cov)

##########
# scalar and (partial) derivative
##########

def cov(l1, l2, par):
  '''Covariance in scalar
  
  :param l1: *l1* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  :param l2: *l2* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  :param par: parameters in dict.
  
  :return: scalar of the covariance.
  :rtype: float
  '''
  covariance = cov_yy(l1, l2)
  value = poly2num(covariance, par)
  return(value)

def dcov(l1, l2, par, wrt):
  '''Partial derivative of covariance w.r.t. parameter wrt
  
  :param l1: *l1* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  :param l2: *l2* in :math:`cov(y_n^{l_1},y_{n+1}^{l_2})`.
  :param par: parameters in dict.
  :param wrt: with respect to.
  
  :return: scalar of the partial derivative.
  :rtype: float
  '''
  covariance = cov_yy(l1, l2)
  pold = dpoly(covariance, wrt)
  value = poly2num(pold, par)
  return(value)


if __name__ == "__main__":
  # test the module
  from pprint import pprint
  # 
  keyfor = ('(n_1m*k1+n_2m*k2)^{-i_m},...,(n_11*k1+n_21*k2)^{-i_1}',
     'e^{-(n1*k1+n2*k2)h}','h','mu', 'theta1','sigma_v1','theta2','sigma_v2')
  print(f"cov_yy() returns poly with keyfor = {keyfor}")
  # # usually l1 >= 1 and l2 >= 1
  # print("cov_yy(l1=0,l2=0): "); pprint(cov_yy(l1=0,l2=0))
  # print("cov_yy(l1=0,l2=1): "); pprint(cov_yy(l1=0,l2=1))
  # print("cov_yy(l1=1,l2=0): "); pprint(cov_yy(l1=1,l2=0))
  # # 
  print("cov_yy(l1=1,l2=1): "); pprint(cov_yy(l1=1,l2=1))
  print("cov_yy(l1=2,l2=1): "); pprint(cov_yy(l1=2,l2=1))
  print("cov_yy(l1=1,l2=2): "); pprint(cov_yy(l1=1,l2=2))
  # print("cov_yy(l1=3,l2=2): "); pprint(cov_yy(l1=3,l2=2))

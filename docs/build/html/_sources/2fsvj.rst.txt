========================
Two-Factor SV with jump
========================

.. math::
   
    d\log s(t) &= (\mu- v(t)/2) dt + \sqrt{v(t)}dw(t) + jdN(t),\\
    v(t)       &= v_1(t) + v_2(t),\\
    dv_1(t)    &= k_1(\theta_1 - v_1(t))dt + \sigma_{1v} \sqrt{v_1(t)}dw_1(t),\\
    dv_2(t)    &= k_2(\theta_2 - v_2(t))dt + \sigma_{2v} \sqrt{v_2(t)}dw_2(t),

where :math:`w_1(t)`, :math:`w_2(t)` are two independent Wiener processes, 
which are also independent of :math:`w(t)`. And :math:`n(t)` is a Poisson 
process with rate :math:`\lambda` and 
independent of everything else, :math:`j` is a random variable distributed 
according to a normal with mean 
:math:`\mu_j` and variance :math:`\sigma_j^2`.

We have :math:`y_n = y_{o,n} + J_n` where

.. math::
   
   y_{o,n} = \mu h - \frac{1}{2}IV_n + I_n^{*},
   \quad
   J_n \triangleq \sum_{i=N((n-1)h)+1}^{N(nh)}j_i. 

.. note:: If we want to estimate the parameters through *Method of Moments*, 
   we need following quantities,
   
   .. math::
      
      E[y_n], var(y_n), cov(y_n,y_{n+1}), cov(y_n,y_{n+2}), cov(y_n^2,y_{n+1}),\\
      cov(y_n,y_{n+1}^2), cm_3[y_n],cov(y_n^3,y_{n+1}), cov(y_n^2,y_{n+1}^2),
      cov(y_n, y_{n+1}^3).


Central Moments
================

Similarly, I define :math:`\overline{y}_n \triangleq y_n - E[y_n]` and we have

.. math::
   
   E[\overline{y}_n^l]
   = \sum_{i=0}^l C_l^i E[\overline{y}_{o,n}^i] E[\overline{J}_n^{l-i}],

where :math:`E[\overline{y}_{o,n}]= y_{o,n} - E[y_{o,n}]` and
:math:`E[\overline{J}_n]= J_n - E[J_n]`

In summary, I defined

1. :py:func:`~hsvmoment.mdl_2fsvj.cmoment.cmoment_y`.

Moments
========

.. math::
   
   E[y_n^l]
   = \sum_{i=0}^l C_l^i E[y_{o,n}^i] E[J_n^{l-i}].

I defined

1. :py:func:`~hsvmoment.mdl_2fsvj.moment.moment_y`.

Covariances
============

.. math::
   
   cov(y_n^{l_1}, y_{n+1}^{l_2})
   = E[y_n^{l_1}y_{n+1}^{l_2}] - E[y_n^{l_1}]E[y_{n+1}^{l_2}].

.. math::
   
   E[y_n^{l_1}y_{n+1}^{l_2}]
   &= E[(y_{o,n}+J_n)^{l_1}(y_{o,n+1}+J_{n+1})^{l_2}]\\
   &= \sum_{i=0}^{l_1}C_{l_1}^i \sum_{j=0}^{l_2}C_{l_2}^j 
   E[y_{o,n}^i J_n^{l_1-i}y_{o,n+1}^j J_{n+1}^{l_2-j}]\\
   &= \sum_{i=0}^{l_1}\sum_{j=0}^{l_2}C_{l_1}^i C_{l_2}^j
   E[y_{o,n}^iy_{o,n+1}^j]E[J_n^{l_1-i}] E[J_{n+1}^{l_2-j}]

In summary, I defined

1. :py:func:`~hsvmoment.mdl_2fsvj.cov.moment_yy`,

2. :py:func:`~hsvmoment.mdl_2fsvj.cov.cov_yy`.


API
====

.. autosummary::
   :toctree: generated
   
   hsvmoment.mdl_2fsvj.cmoment
   hsvmoment.mdl_2fsvj.moment
   hsvmoment.mdl_2fsvj.cov

.. automodule:: hsvmoment.mdl_2fsvj.moment
   :members:

.. automodule:: hsvmoment.mdl_2fsvj.cmoment
   :members:

.. automodule:: hsvmoment.mdl_2fsvj.cov
   :members:

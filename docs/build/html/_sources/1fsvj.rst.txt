========================
One-Factor SV with jump
========================

The first extension is the following SV model, which adds a jump component in 
the log price process of the Heston model: 

.. math::
   
   d\log s(t) &= (\mu- v(t)/2) dt + \sqrt{v(t)}dw^s(t) + jdN(t),\\
   dv(t)      &= k(\theta - v(t))dt + \sigma_v \sqrt{v(t)}dw^v(t),

where :math:`N(t)` is a Poisson process with rate :math:`\lambda` and 
independent of everything else, :math:`j` is a random variable distributed 
according to :math:`F_j(\cdot,\boldsymbol{\theta}_j)` and 
:math:`\boldsymbol{\theta}_j` is a parameter. We take normal with mean 
:math:`\mu_j` and variance :math:`\sigma_j^2` as an example of 
:math:`F_j(\cdot,\boldsymbol{\theta}_j)`. For this model,

.. math::
   
   y_n = y_{o,n} + J_n,

where

.. math::
   
  y_{o,n} \triangleq \mu h - \frac{1}{2}IV_{n} + \rho I_n + 
  \sqrt{1-\rho^2}I_n^{*}, \qquad
  J_n \triangleq \sum_{i=N((n-1)h)+1}^{N(nh)}j_i.

.. note:: If we want to estimate the parameters using *Method of Moments*, 
   we have three more parameters :math:`\mu_j, \sigma_j^2` 
   and :math:`\lambda` to estimate.
   Therefore, we need following eight equations: 
   
   .. math::
      
      E [y_n], var(y_n), cov(y_n,y_{n+1}), cov(y_n,y_{n+2}),\\
      cov(y_n^2,y_{n+1}), cov(y_n,y_{n+1}^2), cm_3[y_n], cm_4[y_n]
   
   where :math:`cm_3[\cdot], cm_4[\cdot]` denotes the third and fourth central 
   moments respectively. 

Moments
========

Moments and Central Moments

.. math::
   
   E[\overline{y}_{n}^l] 
   &= E[(\overline{y}_{o,n} + \overline{J}_n)^l]
   = \sum_{i=0}^{l} C_l^i E[\overline{y}_{o,n}^i]E[\overline{J}_n^{l-i}],\\
   E[y_n^l]
   &= E[(y_{o,n} + J_n)^l]
   = \sum_{i=0}^{l} C_l^i E[y_{o,n}^i] E[J_n^{l-i}].

Functions :py:func:`~hsvmoment.mdl_1fsv.moment.moment_y` and 
:py:func:`~hsvmoment.mdl_1fsv.cmoment.cmoment_y` can be used to compute
:math:`E[y_{o,n}^i]` and
:math:`E[\overline{y}_{o,n}^i]`, respectively.
Meanwhile, functions :py:func:`~hsvmoment.cpp_mmnt.mcpp` and
:py:func:`~hsvmoment.cpp_mmnt.cmcpp` can be used to compute
:math:`E[J_n^{l-i}]` and
:math:`E[\overline{J}_n^{l-i}]`, respectively.

Covariances
============

.. math::
   
   cov(y_n^{l_1}, y_{n+1}^{l_2})
   = E[y_n^{l_1}y_{n+1}^{l_2}] - E[y_n^{l_1}]E[y_{n+1}^{l_2}]

which reduces to 

.. math::
   
   &E[y_n^{l_1}y_{n+1}^{l_2}]\\
   &= \sum_{i=0}^{l_2}C_{l_2}^i E[y_n^{l_1}y_{o,n+1}^i]E[J_{n+1}^{l_2-i}]\\
   &= \sum_{i=0}^{l_2}C_{l_2}^i \sum_{j=0}^{l_1}C_{l_1}^j 
   E[y_{o,n}^jy_{o,n+1}^i] E[J_n^{l_1-j}]E[J_{n+1}^{l_2-i}].

Function :py:func:`~hsvmoment.mdl_1fsv.cov.moment_yy` in module 
:py:mod:`hsvmoment.mdl_1fsv.cov` can be used to compute
:math:`E[y_{o,n}^jy_{o,n+1}^i]`.

In summary, I defined

1. :py:func:`~hsvmoment.mdl_1fsvj.moment.moment_y` for moment :math:`E[y_n^l]`.

2. :py:func:`~hsvmoment.mdl_1fsvj.cmoment.cmoment_y` for central moment 
   :math:`E[\overline{y}_{n}^l]`.

3. :py:func:`~hsvmoment.mdl_1fsvj.cov.cov_yy` for covariance 
   :math:`cov(y_n^{l_1}, y_{n+1}^{l_2})`.


API
====

.. autosummary::
   :toctree: generated
   
   hsvmoment.mdl_1fsvj.cmoment
   hsvmoment.mdl_1fsvj.moment
   hsvmoment.mdl_1fsvj.cov

.. automodule:: hsvmoment.mdl_1fsvj.moment
   :members:

.. automodule:: hsvmoment.mdl_1fsvj.cmoment
   :members:

.. automodule:: hsvmoment.mdl_1fsvj.cov
   :members:

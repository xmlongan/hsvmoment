===========
Quickstart
===========

Make sure ``hsvmoment`` is installed.

Four :abbr:`SV(Stochastic Volatility)` Models
==============================================

The so-called **Heston Stochastic Volatility Model(baseline model)** [#f1]_ in 
this package refers to following SDEs [#f2]_ ,

.. math::
    ds(t) &= \mu s(t)dt + \sqrt{v(t)}s(t)dw^s(t),\\
    dv(t) &= k(\theta - v(t))dt + \sigma_v\sqrt{v(t)}dw^v(t),

where :math:`s(t)` is the asset price at time :math:`t`. Refers to :doc:`theory`
for details.

We also denote the baseline model as ``mdl_1fsv``. Besides, there are three 
extendsions. In summary, we consider four models:

+--------------+----------------------------------------------------------------------------+
| Model        |    Description                                                             |
+==============+============================================================================+
|mdl_1fsv      | - baseline model                                                           |
|              | - refers to :doc:`theory` or :doc:`1fsv`                                   |
+--------------+----------------------------------------------------------------------------+
|mdl_1fsvj     | - with jump in return process                                              |
|              | - refers to :doc:`1fsvj`                                                   |
+--------------+----------------------------------------------------------------------------+
|mdl_2fsv      | - with volatility consisting of superposition of two Square-Root Diffusions|
|              | - refers to :doc:`2fsv`                                                    |
+--------------+----------------------------------------------------------------------------+
|mdl_2fsvj     | - with jump in return process of model mdl_2fsv                            |
|              | - refers to :doc:`2fsvj`                                                   |
+--------------+----------------------------------------------------------------------------+

The calculation of (central) moments and covariances of the four models are
implemented in the four subpackages of package :code:`hsvmoment` 
respectively as

+--------------+-----------------------------+--------------------------------------+
| Model        | Subpackage                  | Modules                              |
+==============+=============================+======================================+
|mdl_1fsv      |  :code:`hsvmoment.mdl_1fsv` | :py:mod:`hsvmoment.mdl_1fsv.cmoment` |
|              |                             | :py:mod:`hsvmoment.mdl_1fsv.moment`  |
|              |                             | :py:mod:`hsvmoment.mdl_1fsv.cov`     |
+--------------+-----------------------------+--------------------------------------+
|mdl_1fsvj     |  :code:`hsvmoment.mdl_1fsvj`| :py:mod:`hsvmoment.mdl_1fsvj.cmoment`|
|              |                             | :py:mod:`hsvmoment.mdl_1fsvj.moment` |
|              |                             | :py:mod:`hsvmoment.mdl_1fsvj.cov`    |
+--------------+-----------------------------+--------------------------------------+
|mdl_2fsv      |  :code:`hsvmoment.mdl_2fsv` | :py:mod:`hsvmoment.mdl_2fsv.cmoment` |
|              |                             | :py:mod:`hsvmoment.mdl_2fsv.moment`  |
|              |                             | :py:mod:`hsvmoment.mdl_2fsv.cov`     |
+--------------+-----------------------------+--------------------------------------+
|mdl_2fsvj     |  :code:`hsvmoment.mdl_2fsvj`| :py:mod:`hsvmoment.mdl_2fsvj.cmoment`|
|              |                             | :py:mod:`hsvmoment.mdl_2fsvj.moment` |
|              |                             | :py:mod:`hsvmoment.mdl_2fsvj.cov`    |
+--------------+-----------------------------+--------------------------------------+


Moment Computation
===================

Take model ``mdl_1fsv`` (baseline) as an example, I demonstrate how to use
functions provided by the package ``hsvmoment`` to compute
the (central) moments and covariances, and their partial derivatives with
respect to a given parameter.

The *n*\ th central moment and moment of return (of the 
*i*\ th interval of length :math:`h`) :math:`y_i` can be computed as

.. code-block:: python
   
   ## Central Moment
   from hsvmoment.mdl_1fsv.cmoment import cm, dcm
   
   parameters = {'mu':0.125, 'k':0.1, 'theta':0.25, 'sigma_v':0.1, 
     'rho':-0.7, 'h': 1}
   # 3rd central moment
   cmoment = cm(l=3, par=parameters)
   # partial derivative wrt(with respect to) parameter 'k'
   dcmoment = dcm(l=3, par=parameters, wrt='k')
   
   ## Moment
   from hsvmoment.mdl_1fsv.moment import m, dm
   
   # 3rd moment
   moment = m(l=3, par=parameters)
   # partial derivative wrt(with respect to) parameter 'k'
   dmoment = dm(l=3, par=parameters, wrt='k')


Covariance Computation
========================

The covariance between two successive returns of power :math:`l_1` and 
:math:`l_2`, i.e., :math:`cov(y_n^{l_1}, y_{n+1}^{l_2})`, can be computed as

.. code-block:: python
   
   ## Covariance
   from hsvmoments.mdl_1fsv.cov import cov, dcov
   
   parameters = {'mu':0.125, 'k':0.1, 'theta':0.25, 'sigma_v':0.1, 
     'rho':-0.7, 'h': 1}
   # covariance cov(y_n^2, y_{n+1}^2)
   covariance = cov(l1=2, l2=2)
   # partial derivative wrt(with respect to) parameter 'k'
   dcovariance = dcov(l1=2, l2=2, par=parameters, wrt='k')


The corresponding quantities for other models (mdl_1fsvj, mdl_2fsv, mdl_2fsvj)
can be computed by using the counterparts within their subpackages.

----------

.. [#f1] Whose exact equation varies according to different authors. One alternative setting is :math:`dp(t) = \mu dt + \sqrt{v(t)}dw^s(t)` where :math:`p(t) = \log s(t)`. 
.. [#f2] Stochastic Differential Equations. Notations: :math:`v(t)` is the instantaneous return variance at time :math:`t`, and :math:`w^s(t)` and :math:`w^v(t)` are two Wiener processes with correlation :math:`\rho`. 

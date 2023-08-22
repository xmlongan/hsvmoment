==============
One-Factor SV 
==============

The considered SV model is that in :doc:`theory`, i.e., 

.. math::
    ds(t) &= \mu s(t)dt + \sqrt{v(t)}s(t)dw^s(t),\\
    dv(t) &= k(\theta - v(t))dt + \sigma_v\sqrt{v(t)}dw^v(t),

which is also the baseline SV model of the package.

The computation of its moments has been explained in :doc:`design` page.


API
====

.. autosummary::
   :toctree: generated
   
   hsvmoment.mdl_1fsv.cmoment
   hsvmoment.mdl_1fsv.moment
   hsvmoment.mdl_1fsv.cov

.. automodule:: hsvmoment.mdl_1fsv.moment
   :members:

.. automodule:: hsvmoment.mdl_1fsv.cmoment
   :members:

.. automodule:: hsvmoment.mdl_1fsv.cov
   :members:

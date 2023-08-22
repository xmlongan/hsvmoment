'''
Class Poly derived from :code:`UserDict`

I defined a new class :code:`Poly` to extend :class:`~collections.UserDict` 
in Python Standard Library :code:`collections`.

I defined an attribute :py:attr:`~hsvmoment.poly.Poly.keyfor` for class 
:py:class:`~hsvmoment.poly.Poly`.
I also defined following magic methods:

+-------------+---------------------+-----------------------------------------------+
|Magic Method | Example             |Function                                       |
+=============+=====================+===============================================+
|__add__()    |:code:`poly1 + poly2`|Add another poly with the same keyfor          |
+-------------+---------------------+-----------------------------------------------+
|__sub__()    |:code:`poly1 - poly2`|subtract another poly with the same keyfor     |
+-------------+---------------------+-----------------------------------------------+
|__mul__()    |:code:`poly1 * poly2`|Multiply with another poly with the same keyfor|
+-------------+---------------------+-----------------------------------------------+
|__rmul__()   |:code:`c * poly2`    |Reversely multiply a constant                  |
+-------------+---------------------+-----------------------------------------------+

Besides, I defined following methods:

+---------------------------------------------+------------------------------------------------------------------+
|Method                                       |Function                                                          |
+=============================================+==================================================================+
|:py:meth:`~hsvmoment.poly.Poly.merge`        |Merge with another poly with same keyfor                          |
+---------------------------------------------+------------------------------------------------------------------+
|:py:meth:`~hsvmoment.poly.Poly.add_keyval`   |Insert a new key-val or add val to the existing key               |
+---------------------------------------------+------------------------------------------------------------------+
|:py:meth:`~hsvmoment.poly.Poly.remove_zero`  |Remove any item with 0 value                                      |
+---------------------------------------------+------------------------------------------------------------------+
|:py:meth:`~hsvmoment.poly.Poly.set_keyfor`   |Set the ``keyfor`` attribute for the poly                         |
+---------------------------------------------+------------------------------------------------------------------+
|:py:meth:`~hsvmoment.poly.Poly.mul_poly`     |Multiply with other poly and return a new poly                    |
+---------------------------------------------+------------------------------------------------------------------+
|:py:meth:`~hsvmoment.poly.Poly.is_exact_type`|Check whether the supplied argument is a Poly with the same keyfor|
+---------------------------------------------+------------------------------------------------------------------+
'''
from collections import UserDict

def kv(i, key):
  '''get key[i] where key is a tuple
  
  :param i: index of key or -1 for not having the element.
  :param key: tuple as a key for Poly.
  
  :return: key[i] if i in range(len(k)), otherwise 0.
  :rtype: int
  '''
  if i in range(len(key)):
    return(key[i])
  else:
    return(0)

class Poly(UserDict):
  '''Class for different versions of "Polynomial" '''
  
  keyfor = ()
  "attribute Poly.keyfor: purpose for each key component, a tuple of str."
  
  def __add__(self, other):
    '''Add another poly with the same keyfor'''
    if not self.is_exact_type(other):
      msg = "The operands must be Polys with the same 'keyfor' attribute."
      raise NotImplementedError(msg)
    # 
    poly = Poly(self) # UserDict initialization
    poly.set_keyfor(self.keyfor)
    for k in other:
      poly.add_keyval(k, other[k])
    return(poly)
  
  def merge(self, other):
    '''Merge another poly with the same keyfor
    
    Insert new key-value or add value for the existing key.
    
    :param other: poly with the same keyfor.
    
    :return: the updated first poly.
    :rtype: Poly
    '''
    if not self.is_exact_type(other):
      msg = "The operands must be Polys with the same 'keyfor' attribute."
      raise NotImplementedError(msg)
    #
    for k in other:
      self.add_keyval(k, other[k])
  
  def __sub__(self, other):
    '''Subtract another poly with the same keyfor'''
    if not self.is_exact_type(other):
      msg = "The operands must be Polys with the same 'keyfor' attribute."
      raise NotImplementedError(msg)
    # 
    poly = Poly(self) # UserDict initialization
    poly.set_keyfor(self.keyfor)
    for k in other:
      if k in self:
        poly[k] -= other[k]
      else:
        poly[k] = -other[k]
    return(poly)
  
  def __mul__(self, other):
    '''Multiply with another poly with the same keyfor'''
    if not self.is_exact_type(other):
      msg = "The operands must be Poly with the same 'keyfor' attribute."
      raise NotImplementedError(msg)
    # 
    poly = Poly()
    poly.set_keyfor(self.keyfor)
    for k1 in self:
      for k2 in other:
        key = tuple(k1[i]+k2[i] for i in range(len(k1)))
        val = self[k1] * other[k2]
        poly.add_keyval(key, val)
    return(poly)
  
  def __rmul__(self, c):
    '''Reversely multiply a constant'''
    for k in self:
      self[k] *= c
    return(self)
  
  def add_keyval(self, key, val):
    '''Insert a new key-val or add val to the existing key
    
    :param key: key for the poly, a tuple.
    :param val: corresponding value.
    '''
    if key in self:
      self[key] += val
    else:
      self[key]  = val
  
  def remove_zero(self):
    '''Remove any item with 0 value'''
    # val is n/m, fraction number
    ks = [k for k in self if self[k] == 0] # works for fraction number
    for k in ks:
      del self[k]
  
  def set_keyfor(self, names):
    '''Set the ``keyfor`` attribute for the poly
    
    :param names: a sequence of names, 
       each corresponding to the key tuple counterpart,
       i.e., names[i] v.s. key[i].
    '''
    self.keyfor = tuple(names)
  
  def mul_poly(self, other, keyIndexes, keyfor):
    '''Multiply with a different poly and return a new poly
    
    :param other: another poly with different 'keyfor' attribute.
    :param keyIndexes: a tuple with two lists, 
       keyIndexes[0] for self, 
       keyIndexes[1] for other.
    :param keyfor: the 'keyfor' for the returned poly.
    
    :return: a poly with 'keyfor'.
    :rtype: Poly
    '''
    poly = Poly()
    poly.set_keyfor(keyfor)
    # 
    idx1, idx2 = keyIndexes
    kN = len(idx1)
    # 
    for k1 in self:
      for k2 in other:
        key = tuple(kv(idx1[i],k1) + kv(idx2[i],k2) for i in range(kN))
        val = self[k1] * other[k2]
        poly.add_keyval(key, val)
    return(poly)
  
  def is_exact_type(self, other):
    '''Check whether the supplied argument is a Poly with the same keyfor'''
    flag = True
    if not isinstance(other, Poly):
      flag = False
    else:
      try:
        if len(self.keyfor) != len(other.keyfor):
          flag = False
        else:
          for i in range(len(self.keyfor)):
            if self.keyfor[i] != other.keyfor[i]:
              flag = False
              break
      except AttributeError:
        flag = False
    return(flag)

if __name__ == "__main__":
  # test class Poly
  from pprint import pprint
  # 
  pol1 = Poly({(1,0): 1, (0,1): 1})
  pol1.set_keyfor(['e', 'h'])
  pol2 = Poly({(1,0): 1, (0,2): 1})
  pol2.set_keyfor(['e', 'h'])
  # add
  pol3 = pol1 + pol2; print('pol3 = pol1 + pol2: '); pprint(pol3)
  pol4 = pol2 + pol1; print('pol4 = pol2 + pol1: '); pprint(pol4)
  # sub
  pol5 = pol1 - pol2; print('pol5 = pol1 - pol2: '); pprint(pol5)
  # mul
  4 * pol1; print('4 * pol1'); pprint(pol1)
  pol6 = pol1 * pol2; print('pol6 = pol1 * pol2: '); pprint(pol6)
  #
  pol7 = Poly({(11,0): 2, (0,11): 1})
  # pol8 = pol1 * pol7 # will raise error!
  # merge
  pol1.merge(pol2); print('pol1.merge(pol2): '); pprint(pol1)
  # add_keyval
  pol2.add_keyval(key=(2,0), val=1)
  pol2.add_keyval(key=(0,2), val=100); pprint(f'add_keyval: {pol2}')
  # remove_zero
  pol2[(3,1)] = 0; pol2.remove_zero(); pprint(f'remove_zero: {pol2}')
  # set_keyfor
  print('pol2.keyfor:'); pprint(pol2.keyfor)
  # mul_poly
  poln = Poly({(1,):100, (2,):100})
  print('before mul_poly, pol1: '); pprint(pol1)
  poln = pol1.mul_poly(poln, keyIndexes=([0,1],[0,-1]), keyfor=('e','h'))
  print('mul_poly: '); pprint(poln)
  del pol1, pol2, pol3, pol4, pol5, pol6, pol7, poln

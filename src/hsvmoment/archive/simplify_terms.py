'''
Module collecting together term-simplifying functions

Simplify terms like ' :math:`x^0`' -> '', ' :math:`x^1`' -> 'x'. 
Adding parenthese necessarily.
'''

def p(n):
  '''
  Add parenthese or not
  
  :param n: expression.
  :type string: str
  
  :return: an expression surrounded by parenthese when necessary.
  :rtype: string
  '''
  if n >= 0 and n < 10:
    return(f'{n}')
  else:
    return(f'({n})')

def is_expr(term):
  '''
  Check whether a term is a expression containing operators
  
  :param term: a term or expression.
  :type string: str
  
  :return: True or False.
  :rtype: bool
  '''
  operators = ['+','-','*','/','^','(',')']
  flag = False
  for op in operators:
    if op in term:
      flag = True
      break
  return(flag)

def times_pow(base,n):
  '''
  Multiply a power term in proper string
  
  :param base: the base x in :math:`x^n`.
  :type string or number: str or number
  :param n: the exponent n in :math:`x^n`.
  :type real number: int or float
  
  :return: an expression in string.
  :rtype: str
  '''
  if is_expr(base):
    base = f'({base})'
  if n == 0:
    return('')
  elif n == 1:
    return(f' * {base}')
  elif n < 10:
    return(f' * {base}^{n}')
  else:
    return(f' * {base}^({n})')

def times_exp(n,t):
  '''
  multiply an exponetial term in proper string
  
  :param n: the n in :math:`e^{nt}`.
  :type number: int or float
  :param t: the t in :math:`e^{nt}`, t > 0.
  :type string or number: str or int or float
  
  :return: an expression in string.
  :rtype: str
  '''
  if n == 0:
    return('')
  elif n == 1:
    return(f' * exp({t})')
  else:
    return(f' * exp({n}*{t})')

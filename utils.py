'''
Utility functions
'''


def mprint(a):
    '''
    Prints a matrix 

    Input string with variable name
    '''
    print('\n {} = '.format(a))
    print(eval(a, globals()))


class lazy_property(object):
    '''
    meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    '''

    def __init__(self, fget):
        self.fget = fget
        self.func_name = fget.__name__


    def __get__(self, obj, cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj, self.func_name, value)
        return value

class lazy(property):
    pass

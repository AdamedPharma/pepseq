"""
[[4, 3], [-2, 1]]
"""


class Polynomial(object):
    def __init__(self, params):
        # poly = Polynomial([[4, 3], [-2, 1]])
        # exponents = [term[1] for term in poly.terms]
        self.params = params
        # self.terms = [i[0] for i in self.params]
        self.terms = self.params

        return

    def derivative(self):
        return self.__class__(self.params)

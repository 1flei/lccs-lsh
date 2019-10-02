
def getr(dataset):
    rs_key = {'Mnist':300, 'Sift':200}
    return rs_key[dataset.name]

class MNIST:
    def __init__(self):
        self.n = 60000
        self.qn = 1000
        self.d = 50
        # [293.8953705 328.277313  348.174515  360.718002  370.533325  377.599426
        #     384.100876  390.465744  396.036606  401.089386 ]
        self.r = 401
        self.name = 'Mnist'

class MNIST784:
    def __init__(self):
        self.n = 60000
        self.qn = 1000
        self.d = 784
        self.r = 1270
        self.name = 'Mnist784'
        
class Trevi:
    def __init__(self):
        self.n = 100800
        self.qn = 100
        self.d = 4096
        self.name = 'Trevi'

        self.r = 1600
        
class Netflix:
    def __init__(self):
        self.n = 17500
        self.qn = 270
        self.d = 300
        self.name = 'Netflix'
        
class Sift:
    def __init__(self):
        self.n = 1000000
        self.qn = 1000
        self.d = 128
        # [  0.        152.3875655 189.3990095 201.3231125 210.819786  215.4588015
        #      218.8215715 221.551346  223.794548  225.594757 ]
        self.r = 226
        self.name = 'Sift'

class Sift10M:
    def __init__(self):
        self.n = 11164866-100
        self.qn = 100
        self.d = 128
        self.r = 200
        self.name = 'Sift10M'
        
class Yahoo:
    def __init__(self):
        # self.n = 624961
        self.n = 620000
        # self.qn = 4961
        self.qn = 1000
        self.d = 300
        self.name = 'Yahoo'
        
class Gist:
    def __init__(self):
        self.n = 1000000
        self.qn = 1000
        self.d = 960
        self.r = 11294
        self.name = 'Gist'

class MovieLens:
    def __init__(self):
        # self.n = 50000
        # self.qn = 3889
        self.n = 52889
        self.qn = 1000
        self.d  = 150
        self.name = 'MovieLens'

class Glove:
    def __init__(self):
        self.n = 2196017-1000
        self.qn = 1000
        self.d  = 300
        self.name = 'glove'
        self.r =  8

class Synthetic:
    def __init__(self):
        self.n = 1000000
        self.qn = 1000
        self.d  = 256
        self.name = 'synthetic'
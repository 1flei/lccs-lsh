

class LCCS:
    def __init__(self, Ls=[8, 16, 32, 64, 128, 256, 512], step=1):
        self.Ls = Ls
        self.step = step
        self.name = 'lccs'

    def for_param(self, distance):
        for l in self.Ls:
            yield '-L %d --step 1'%l

class LCCS_MP:
    def __init__(self, Ls=[8, 16, 32, 64, 128, 256, 512], ps=[0, 1, 2, 4, 8]):
        self.Ls = Ls
        self.ps = ps
        self.name = 'lccs_mp'

    def for_param(self, distance):
        for l in self.Ls:
            for p in self.ps:
                yield '-L %d -p %d'%(l, p)

class C2LSH:
    def __init__(self, Ls=[8, 16, 32, 64, 128, 256, 512]):
        self.Ls = Ls
        self.name = 'c2lsh'

    def for_param(self, distance):
        for l in self.Ls:
            yield '-L %d'%l

class E2LSH:
    def __init__(self, Ks=[3, 4, 5, 6, 7, 8, 9, 10], Ls=[8, 16, 32, 64, 128, 256, 512], maxHasher=512):
        self.Ks = Ks
        self.Ls = Ls
        self.name = 'e2lsh'
        self.maxHasher = maxHasher

    def for_param(self, distance):
        for l in self.Ls:
            for k in self.Ks:
                if k*l <=self.maxHasher:
                    yield '-L %d -K %d'%(l, k)

class MPLSH:
    def __init__(self, Ks=[5, 6, 7, 8, 9, 10], Ls=[8, 16, 32, 64, 128, 256, 512], maxHasher=512, Ks_angle=[2, 3, 4, 5, 6, 7, 8, 9, 10]):
        self.Ks = Ks
        self.Ls = Ls
        self.name = 'mplsh'
        self.maxHasher = maxHasher
        self.Ks_angle = Ks_angle

    def for_param(self, distance):
        if distance=='angle':
            for r in self.for_param_angle(distance):
                yield r 
        else:
            for r in self.for_param_l2(distance):
                yield r

    def for_param_l2(self, distance):
        for l in self.Ls:
            for k in self.Ks:
                if k*l <=self.maxHasher:
                    yield '-L %d -K %d'%(l, k)

    # def for_param_angle(self, distance):
    #     #bit-packed version
    #     Ks = self.Ks if self.Ks_angle is None else self.Ks_angle
    #     for l in self.Ls:
    #         for k in Ks:
    #             yield '-L %d -K %d --bit_packed'%(l, k)
    def for_param_angle(self, distance):
        #linear-probing version
        Ks = self.Ks if self.Ks_angle is None else self.Ks_angle
        for l in self.Ls:
            for k in Ks:
                if k*l <=self.maxHasher:
                    yield '-L %d -K %d'%(l, k)
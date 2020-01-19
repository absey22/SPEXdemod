class MOR:
    """
    comp=None
    thick=None
    par=None
    rot=None
    temperature=None
    """
    def __init__(self,composition=['al2o3','mgf2'],thickness=[1.2,2.88],\
                 par=None,rot=None,temp=None ):
        self.comp=composition
        self.thick=thickness
        self.par=par
        self.rot=rot
        self.temperature=temp

    def add_layer(self,newlayer):
        self.comp.append(newlayer)

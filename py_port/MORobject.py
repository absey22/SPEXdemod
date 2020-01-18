class MOR:
    """
    comp=None
    thick=None
    par=None
    rot=None
    temperature=None
    """
    def __init__(self,composition,thickness,par=None,rot=None,temp=None):
        self.comp=composition
        self.thick=thickness
        self.par=par
        self.rot=rot
        self.temperature=temp

    def add_layer(self,newlayer):
        self.comp.append(newlayer)

import os
import importlib

__version__ = "0.1"

class McAstroIS:
    """
    McAstro import status
    """
    def __init__(self):
        """
        Try to import all packages contained in self.packages and save status.
        """
        self.packages = []
        self.findpackages()
        self.importStatus = {}
        for m in self.packages:
            try:
                importlib.import_module("." + m, "McAstro")
                self.importStatus[m] = True
            except Exception:
                self.importStatus[m] = False
            
    def findpackages(self):
        for root, dirs, fns in os.walk(__path__[0]):
            for dir_ in dirs:
                if os.path.isfile(os.path.join(__path__[0], dir_,
                                               "__init__.py")):
                    self.packages.append(dir_)
            break;
    
    def showStatus(self):
        print() 
        print("               McAstro import status")
        print("--------------------------------------------------")
        print("{:41s} | {:6s}".format("Package", "Loaded"))
        print("--------------------------------------------------")
        for m, v in sorted(self.importStatus.items()):
            print("{:41s} | {}".format(m, v))
        

def importCheck():
    m = McAstroIS()
    m.showStatus()

import os
import importlib

__version__ = "0.1"

class WrapperIS:
    """
    Wrapper import status
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
                importlib.import_module("." + m, "wrapper")
                self.importStatus[m] = True
            except Exception:
                self.importStatus[m] = False
            
    def findpackages(self):
        for root, dirs, fns in os.walk(__path__[0]):
            for dir in dirs:
                if os.path.isfile(os.path.join(__path__[0], dir,
                                               "__init__.py")):
                    self.packages.append(dir)
            break;
    
    def showStatus(self):
        print() 
        print("               wrapper import status")
        print("--------------------------------------------------")
        print("{:41s} | {:6s}".format("Package", "Status"))
        print("--------------------------------------------------")
        for m, v in sorted(self.importStatus.items()):
            print("{:41s} | {:6}".format(m, v))
        

def importCheck():
    m = WrapperIS()
    m.showStatus()

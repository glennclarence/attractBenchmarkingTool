
from Configuration import Configuration
import os
def ripExtension(filename):
    base=os.path.basename(filename)
    return os.path.splitext(base)[0]



class Configurator:
    # def __init__(self):
    #     self.config = Configuration()

    def __init__(self,  setting, configfunction): #configuration,
        #self.config = configuration
        self.setting  = setting
        self.configFunc = configfunction
    def setConfig(self, config):
        self.config = config
 
    def completeOutputFile(self,filekey, extensionKey):
        rawName = ripExtension(filename)
        self.config.subPath(self)
        filename = rawName + self.config.getExtension(extensionKey)
        return os
        
    def run(self):
        if not self.config.getSetting(self.setting)['dryRun']:
            if self.config.inputFilesExist(self.setting):
                if not self.config.outputFilesExist(self.setting) or self.config.overwriteOutput(self.setting):
                    if not os.path.isdir(self.config.getOutputFolder(self.setting)):
                        os.system("mkdir -p {} ".format(self.config.getOutputFolder(self.setting)))
                    self.configFunc(self.config,self.setting)
            else:
                self.config.getNotExistingInput(self.setting)
        else:
            self.configFunc(self.config,self.setting)
        #else:
        #    raise ValueError('ConfiguratorBase: input files {} of configuration {} does not exist'.format(self.config.getNotExistingInput(self.setting),self.config.getId() ))

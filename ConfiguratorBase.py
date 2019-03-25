
from Configuration import Configuration
import os
import logging

#ISSUE file is not put into the queue if it exists
#  

def ripExtension(filename):
    base=os.path.basename(filename)
    return os.path.splitext(base)[0]

def createLoggingFile(filename):
    logging.basicConfig(filename=filename,level=logging.DEBUG)

class Configurator:
    # def __init__(self):
    #     self.config = Configuration()

    def __init__(self,  setting, configfunction, logger = None): #configuration,
        #self.config = configuration
        self.setting  = setting
        self.configFunc = configfunction
        self.logger = logger
    def setConfig(self, config):
        self.config = config
 
    def completeOutputFile(self,filekey, extensionKey):
        rawName = ripExtension(filename)
        self.config.subPath(self)
        filename = rawName + self.config.getExtension(extensionKey)
        return os
        
    def run(self):
        if not self.config.getSetting(self.setting)['dryRun']:
            if self.config.inputFilesExist(self.setting) or not self.config.getSetting(self.setting)['checkInput']: #check wether all input files are existing
                if not self.config.outputFilesExist(self.setting) or self.config.overwriteOutput(self.setting) or self.config.outputFileIsEmpty(self.setting): # continue if outputfile not existent, overwriting is allowed or the outputfile is empty
                    if not os.path.isdir(self.config.getOutputFolder(self.setting)):
                        status = os.system("mkdir -p {} ".format(self.config.getOutputFolder(self.setting)))
                        if status != 0:  logging.warning( "could not Create {}".format(self.config.getOutputFolder(self.setting)))
                    try:
                        self.configFunc(self.config,self.setting)
                    except Exception as error:
                        errorstr = "An error occured when running CONFIGURATOR   {:13} ERROR".format(self.setting) + error
                        logging.warning(errorstr)
                else:
                    logstring = "Output file for CONFIGURATOR   {:13s} exists already FILE {}".format(self.setting, self.config.getOutputFile(self.setting, 'out'))
                    logging.warning(logstring)
            else:
                notexist = self.config.getNotExistingInput(self.setting)
                logging.warning("Input files for CONFIGURATOR   {:13s} are not existing. FILES:\n\t\t- {}".format(self.setting,str.join(' ', notexist )))
        else:
            self.configFunc(self.config,self.setting)
        #else:
        #    raise ValueError('ConfiguratorBase: input files {} of configuration {} does not exist'.format(self.config.getNotExistingInput(self.setting),self.config.getId() ))


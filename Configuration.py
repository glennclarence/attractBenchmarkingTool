
import json
import os

def findAndChangeNested(inputDict, k, v):
    """Recursiveley find all entries in a JSON object with key 'k' and change them to value 'v' """
    if type(inputDict) is dict:
        for key in inputDict.keys():
            if key == k:
                inputDict[key] = v
            elif type(inputDict[key]) is dict or type(inputDict[key]) is list:
                findAndChangeNested(inputDict[key], k, v)

    if type(inputDict) is list:
        for entry in inputDict:
            if type(entry) is dict or type(entry) is list:
                findAndChangeNested(entry, k, v)

class Configuration:
    
    def __init__(self,config):
        self.files = config['files']
        self.settings = config["settings"]
        self.id = config["id"]
        self.basePath = self.settings["basePath"]

    def getId(self):
        """Return the Id of the configuration."""
        return self.id

    def subPath(self,folder):
        return os.path.join(self.basePath,folder)

    def getFile(self, fileId):
        try:
            entry = self.files[fileId]
        
        #print("get file", fileId, entry, os.path.join( self.basePath,  os.path.join(entry['folder'], entry['name'] + entry['extension'])))
            return os.path.join( self.basePath,  os.path.join(entry['folder'], entry['name'] + entry['extension']))
        except Exception as error:
            print("getFile: fileId {} not in files. ConfigurationId {}".format(fileId, self.id)," Original message", error)

    def fileExists(self,fileKey):
        """returns true if files exist"""
        return True if os.path.isfile(self.getFile(fileKey)) else False

    def inputFilesExist(self,setting):
        """returns true if all input files exist"""
        
        # if len(self.inputFiles) == 0:
        #     raise ValueError ("CONFIGURATION no inputfiles existing for config {} ".format(self.id))
        inFiles = self.settings[setting]['in']
        exists = True
        for fileid in inFiles.values():
            if not self.fileExists(fileid):
                exists = False
        return exists

    def outputFilesExist(self,setting):
        """returns true if all input files exist"""
        
        # if len(self.inputFiles) == 0:
        #     raise ValueError ("CONFIGURATION no inputfiles existing for config {} ".format(self.id))
        outFiles = self.settings[setting]['out']
        exists = True
        for fileid in outFiles.values():
            if not self.fileExists(fileid):
                exists = False
        return exists

    def getNotExistingInput(self,setting):
        """returns all the inputfiles that do not exist"""
        notExisting = []
        inFiles = self.settings[setting]['in']
        for fileId in inFiles.values():
            if not self.fileExists(fileId):
                notExisting.append(self.getFile(fileId))
        return notExisting

    def getOutputFolder(self,setting):
        fileId = self.settings[setting]['out']['out']
        folder = self.files[fileId]['folder']
        return os.path.join( self.basePath, folder)


    def overwriteOutput(self,setting):
        """return weather output maybe overwritte or not"""
        return self.settings[setting]['overwrite']
    
    def changeAttribute(self, key, value):
        """Find a configuration key and change its value. Slower than addOrChangemembers"""
        findAndChangeNested(self.config, key, value)

    def getSetting(self,key):
        return self.settings[key]

    def setId(self, id):
        self.id = id
    
    def save(self, filename):
        """Saves configuration to a json file"""
        with open(filename, "w") as write_file:
            safeDict = {}
            safeDict['files'] = self.files
            safeDict['settings'] = self.settings

            json.dump(safeDict, write_file)

    def getInputFile(self, setting, inFileId):
        setting = self.settings[setting]
        fileId = setting['in'][inFileId]
        return self.getFile(fileId)

    def getOutputFile(self, setting, inFileId):
        setting = self.settings[setting]
        fileId = setting['out'][inFileId]
        return self.getFile(fileId)
    
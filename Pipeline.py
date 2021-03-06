
import queue
from worker import *
import copy
from  Configurators import configuratorFunctions
from ConfiguratorBase import Configurator

class PipeLine:
    def __init__(self,threads, queues):
        self.threads = threads
        self.queues = queues
        self.size = len(self.queues)
    def join(self, timeout = None):
        for thread in self.threads:
            thread.join(timeout)
    def start(self):
        for thread in self.threads:
            thread.start()
    def stop(self):
        for thread in self.threads:
            thread.stop()
    def put(self,item):
        self.queues[0].put(item)
    def get(self):
        return self.queues[-1].get()

def createPipeline( bufferSize, configurators, numItems):  #inputQueue, outputQueue,
    queues = []


    #queues.append(inputQueue)
    queues.append(queue.Queue(bufferSize))
    for i in range(len(configurators)-1):
        queues.append(queue.Queue(bufferSize))
    queues.append(queue.Queue(bufferSize))
    #queues.append(outputQueue)
    threads = []

    for i, configurator in enumerate(configurators):
        try:
            configurator['numThreads']
            configurator['conf']
        except ValueError as error:
            print ("An element of configurators must be a dictionary that contains a key 'numThreads' followed by the number of threads as well as a key 'conf' followed a single configurator or a list of configurators")
        counter = ThreadCounter(0)
        finishCounter = ThreadCounter(0)
        finishThreshold = configurator['numThreads']
        for k in range(configurator['numThreads']):

            if type(configurator['conf']) != list:
                configuratorList = [configurator['conf']]
            else:
                configuratorList = configurator['conf']
            ctrs =  []
            for con in configuratorList:
                print("create cofigurator ", i, "setting", con['setting'], "with configurator function", con['configurator'], "numThread", k)
                ctrs.append(Configurator(con['setting'],configuratorFunctions[con['configurator']]))

            c = ConsumerThread(configurators = copy.deepcopy(ctrs), threshold = numItems,inputQueue=queues[i],resultQueue = queues[i+1],counter = counter,finishCounter =finishCounter, finishThreshold = finishThreshold,  name= "{}_{}".format(configuratorList[0]['configurator'], k))
            threads.append(c)
    return PipeLine(threads,queues)
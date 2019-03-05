import threading
import time
from multiprocessing import RawValue
import queue

class ThreadCounter(object):
    def __init__(self, value=0):
        self.val = RawValue('i', value)
        self.lock = threading.Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value

class ConsumerThread(threading.Thread):
    def __init__(self,configurators, threshold, inputQueue, resultQueue,counter,finishCounter,finishThreshold, group=None, target=None, name=None,
                 args=(), kwargs=None, verbose=None):
        super(ConsumerThread,self).__init__()
        self.target = target
        self.name = name
        self.counter = counter
        self.threshold = threshold
        self.configurators = configurators
        self.finishCounter = finishCounter
        self.finishThreshold  = finishThreshold
        self.lock = threading.Lock()
        self.kill = threading.Event()
        self.inputQueue = inputQueue
        self.resultQueue = resultQueue
        self.dorun = True
        return

    def run(self):
        while not self.kill.is_set():

            if not self.inputQueue.empty():
                ( item, id ) = self.inputQueue.get()
                
                for configurator in self.configurators:
                    #try:
                        configurator.setConfig(item)
                        configurator.run()
                    #except Exception as error:
                     #   print("Run in thread {} and configurator {} failed ".format(self.name,configurator.setting ))
                self.resultQueue.put((item, id))

                self.counter.increment()

            if self.counter.value() >= self.threshold and self.dorun:
                self.dorun = False
                self.finishCounter.increment()

            if self.finishCounter.value() >= self.finishThreshold:
                self.stop()
                #print("killAll",self.name,self.finishCounter.value())

                return 

                
        return

    def stop(self):
        with self.lock:
            self.kill.set()



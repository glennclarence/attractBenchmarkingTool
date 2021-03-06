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
                        #print("process ",configurator.setting, " ", item.files[item.settings[configurator.setting]['out']['out']]['name'], item.files[item.settings[configurator.setting]['out']['out']]['extension'])

                        configurator.run()
                    #except Exception as error:
                     #   print("Run in thread {} and configurator {} failed ".format(self.name,configurator.setting ))
                self.resultQueue.put((item, id))

                self.counter.increment()
                #print("worker qsize- name:",self.name, self.inputQueue.qsize())
            if self.counter.value() >= self.threshold and self.dorun:
                #print("quit this thread ",self.name, "count",self.counter.value(),"threshold",self.threshold,"input size ",self.inputQueue.qsize(),"finish thresh",self.finishThreshold)
                self.dorun = False
                self.finishCounter.increment()

            if self.finishCounter.value() >= self.finishThreshold:
                self.stop()
                print("Kill Thread - name:",self.name, " Finished Threads:" ,self.finishCounter.value(),"Finished Items: ",self.counter.value(),"Number of target Items:",self.threshold,"Remaning objects: ",self.inputQueue.qsize())
                #return 
            time.sleep(0.01)
                
        return

    def stop(self):
        with self.lock:
            self.kill.set()



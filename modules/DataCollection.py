'''
Created on 21 Feb 2017

@author: jkiesele
'''
#from tensorflow.contrib.labeled_tensor import batch
#from builtins import list
from __future__ import print_function

from Weighter import Weighter
from TrainData import TrainData, fileTimeOut
#for convenience
import logging
from pdb import set_trace
import copy

class DataCollection(object):
    '''
    classdocs
    '''


    def __init__(self, infile = None, nprocs = -1):
        '''
        Constructor
        '''
        self.clear()
        self.nprocs = nprocs       
        self.meansnormslimit=500000 
        if infile:
            self.readFromFile(infile)
        
    def clear(self):
        self.samples=[]
        self.sampleentries=[]
        self.originRoots=[]
        self.nsamples=0
        self.dataDir=""
        self.useweights=True
        self.__batchsize=1
        self.filesPreRead=2
        self.isTrain=True
        self.dataclass=TrainData() #for future implementations
        self.weighter=Weighter()
        self.weightsfraction=0.05
        self.maxConvertThreads=2
        self.maxFilesOpen=5
        self.means=None
        self.classweights={}

    def __iadd__(self, other):
        'A += B'
        if not isinstance(other, DataCollection):
            raise ValueError("I don't know how to add DataCollection and %s" % type(other))
        def _extend_(a, b, name):
            getattr(a, name).extend(getattr(b, name))
        _extend_(self, other, 'samples')
        if len(set(self.samples)) != len(self.samples):
            raise ValueError('The two DataCollections being summed contain the same files!')
        _extend_(self, other, 'sampleentries')
        _extend_(self, other, 'originRoots')
        self.nsamples += other.nsamples
        if self.dataDir != other.dataDir:
            raise ValueError('The two DataCollections have different data directories, still to be implemented!')
        self.useweights = self.useweights and self.useweights
        self.filesPreRead = min(self.filesPreRead, other.filesPreRead)
        self.isTrain = self.isTrain and other.isTrain #arbitrary choice, could also raise exception
        if type(self.dataclass) != type(other.dataclass):
            raise ValueError(
                'The two DataCollections were made with a'
                ' different data class type! (%s, and %s)' % (type(self.dataclass), type(other.dataclass))
                )
        if self.weighter != other.weighter:
            raise ValueError(
                'The two DataCollections have different weights'
                )
        if self.weightsfraction != other.weightsfraction:
            raise ValueError('The two DataCollections have different weight fractions')
        self.maxConvertThreads = min(self.maxConvertThreads, other.maxConvertThreads)
        self.maxFilesOpen = min(self.maxFilesOpen, other.maxFilesOpen)
        if not all(self.means == other.means):
            raise ValueError(
                'The two DataCollections head different means'
                )
        self.classweights.update(other.classweights)
        return self

    def __add__(self, other):
        'A+B'
        if not isinstance(other, DataCollection):
            raise ValueError("I don't know how to add DataCollection and %s" % type(other))
        ret = copy.deepcopy(self)
        ret += other
        return ret

    def __radd__(self, other):
        'B+A to work with sum'
        if other == 0:
            return copy.deepcopy(self)
        elif isinstance(other, DataCollection):
            return self + other #we use the __add__ method
        else:
            raise ValueError("I don't know how to add DataCollection and %s" % type(other))
        
    def removeLast(self):
        self.samples.pop()
        self.nsamples-=self.sampleentries[-1]
        self.sampleentries.pop()
        self.originRoots.pop()
        
        
    def getClassWeights(self):
        if not len(self.classweights):
            self.__computeClassWeights()
        return self.classweights
        
    def __computeClassWeights(self,truthclassesarray):
        if not len(self.samples):
            raise Exception("DataCollection:computeClassWeights: no sample files associated")
        import copy
        td=copy.deepcopy(self.dataclass)
        td.readIn(self.getSamplePath(self.samples[0]))
        arr=td.y[0]
        average=0
        allist=[]
        for i in range(arr.shape[1]):
            entries=float((arr[:,i]>0).sum())
            average=average+entries
            allist.append(entries)
        
        outdict={}
        average=average/float((arr.shape[1]))
        for i in range(len(allist)):
            l=average/allist[i] 
            outdict[i]=l
        self.classweights=outdict
        
    def getInputShapes(self):
        '''
        gets the input shapes from the data class description
        '''
        if len(self.samples)<1:
            return []
        self.dataclass.filelock=None
        td=copy.deepcopy(self.dataclass)
        td.readIn(self.getSamplePath(self.samples[0]))
        shapes=td.getInputShapes()
        td.clear()
        return shapes
    
    def getTruthShape(self):
        return self.dataclass.getTruthShapes()
        
    def getNRegressionTargets(self):
        return (self.dataclass.getNRegressionTargets())
    
    def getNClassificationTargets(self):
        return (self.dataclass.getNClassificationTargets())
        
    def getUsedTruth(self):
        return self.dataclass.getUsedTruth()
    
    def setBatchSize(self,bsize):
        if bsize > self.nsamples:
            raise Exception('Batch size must not be bigger than total sample size')
        self.__batchsize=bsize
        
    def getSamplesPerEpoch(self):
        #modify by batch split
        count=self.getNBatchesPerEpoch()
        if count != 1:
            return count*self.__batchsize #final
        else:
            return self.nsamples
        
    
    def getNBatchesPerEpoch(self):
        if self.__batchsize <= 1:
            return 1
        count=0
        while (count+1)*self.__batchsize <= self.nsamples:
            count+=1
        return count
        
    def writeToFile(self,filename):
        import pickle
        fd=open(filename,'wb')
        self.dataclass.clear()
        pickle.dump(self.samples, fd,protocol=0 )
        pickle.dump(self.sampleentries, fd,protocol=0 )
        pickle.dump(self.originRoots, fd,protocol=0 )
        pickle.dump(self.nsamples, fd,protocol=0 )
        pickle.dump(self.useweights, fd,protocol=0 )
        pickle.dump(self.__batchsize, fd,protocol=0 )
        pickle.dump(self.dataclass, fd,protocol=0 )
        pickle.dump(self.weighter, fd,protocol=0 )
        #pickle.dump(self.means, fd,protocol=0 )
        self.means.dump(fd)
        fd.close()
        
    def readFromFile(self,filename):
        import pickle
        fd=open(filename,'rb')
        self.samples=pickle.load(fd)
        self.sampleentries=pickle.load(fd)
        self.originRoots=pickle.load(fd)
        self.nsamples=pickle.load(fd)
        self.useweights=pickle.load(fd)
        self.__batchsize=pickle.load(fd)
        self.dataclass=pickle.load(fd)
        self.weighter=pickle.load(fd)
        self.means=pickle.load(fd)
        fd.close()
        import os
        self.dataDir=os.path.dirname(os.path.abspath(filename))
        self.dataDir+='/'
        #don't check if files exist
        return 
        for f in self.originRoots:
            if not f.endswith(".root"): continue
            if not os.path.isfile(f):
                print('not found: '+f)
                raise Exception('original root file not found')
        for f in self.samples:
            fpath=self.getSamplePath(f)
            if not os.path.isfile(fpath):
                print('not found: '+fpath)
                raise Exception('sample file not found')
        
        
    def readRootListFromFile(self,file):
        import os
        
        self.samples=[]
        self.sampleentries=[]
        self.originRoots=[]
        self.nsamples=0
        self.dataDir=""
        
        fdir=os.path.dirname(file)
        fdir=os.path.abspath(fdir)
        fdir=os.path.realpath(fdir)
        lines = [line.rstrip('\n') for line in open(file)]
        for line in lines:
            if len(line) < 1: continue
            self.originRoots.append(fdir+'/'+line)
    
        if len(self.originRoots)<1:
            raise Exception('root samples list empty')
        
        
    def split(self,ratio):
        '''
        ratio is self/(out+self)
        returns out
        modifies itself
        '''
        import copy
        
        
        out=DataCollection()
        itself=copy.deepcopy(self)
        
        nsamplefiles=len(self.samples)
        
        out.samples=[]
        out.sampleentries=[]
        out.originRoots=[]
        out.nsamples=0
        out.__batchsize=copy.deepcopy(self.__batchsize)
        out.isTrain=copy.deepcopy(self.isTrain)
        out.dataDir=self.dataDir

        out.dataclass=copy.deepcopy(self.dataclass)
        out.weighter=self.weighter #ref oks
        out.means=self.means
        out.useweights=self.useweights
     
        
        itself.samples=[]
        itself.sampleentries=[]
        itself.originRoots=[]
        itself.nsamples=0
        
        
        
        if nsamplefiles < 2:
            out=copy.deepcopy(self)
            print('DataCollection.split: warning: only one file, split will just return a copy of this')
            return out
    
        for i in range(0, nsamplefiles):
            frac=(float(i))/(float(nsamplefiles))
            if frac < ratio and i < nsamplefiles-1:
                itself.samples.append(self.samples[i])
                itself.sampleentries.append(self.sampleentries[i])
                itself.originRoots.append(self.originRoots[i])
                itself.nsamples+=self.sampleentries[i]
            else:
                out.samples.append(self.samples[i])
                out.sampleentries.append(self.sampleentries[i])
                out.originRoots.append(self.originRoots[i])
                out.nsamples+=self.sampleentries[i]
           
        
        
        self.samples=itself.samples
        self.sampleentries=itself.sampleentries
        self.originRoots=itself.originRoots
        self.nsamples=itself.nsamples
        
        return out
    
    
    def createTestDataForDataCollection(
            self, collectionfile, inputfile, 
            outputDir, outname = 'dataCollection.dc',
            batch_mode = False):
        import copy
        self.readFromFile(collectionfile)
        self.dataclass.remove=False
        self.dataclass.weight=False
        self.readRootListFromFile(inputfile)
        self.createDataFromRoot(
            self.dataclass, outputDir, False,
            dir_check = not batch_mode
        )
        self.writeToFile(outputDir+'/'+outname)
        
    
    def recoverCreateDataFromRootFromSnapshot(self, snapshotfile):
        import os
        snapshotfile=os.path.abspath(snapshotfile)
        self.readFromFile(snapshotfile)
        td=self.dataclass
        #For emergency recover  td.reducedtruthclasses=['isB','isC','isUDSG']
        if len(self.originRoots) < 1:
            return
        #if not self.means:
        #    self.means=td.produceMeansFromRootFile(self.originRoots[0])
        outputDir=os.path.dirname(snapshotfile)+'/'
        self.dataDir=outputDir
        finishedsamples=len(self.samples)
        
        self.__writeData_async_andCollect(finishedsamples,outputDir)
        self.writeToFile(outputDir+'/dataCollection.dc')
        
    
    def createDataFromRoot(
                    self, dataclass, outputDir, 
                    redo_meansandweights=True, means_only=False, dir_check=True
                    ):
        '''
        Also creates a file list of the output files
        After the operation, the object will point to the already processed
        files (not root files)
        Writes out a snapshot of itself after every successfully written output file
        to recover the data until a possible error occurred
        '''
        
        if len(self.originRoots) < 1:
            print('createDataFromRoot: no input root file')
            raise Exception('createDataFromRoot: no input root file')
        
        import os
        outputDir+='/'
        if os.path.isdir(outputDir) and dir_check:
            raise Exception('output dir must not exist')
        elif not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        self.dataDir=outputDir
        self.nsamples=0
        self.samples=[]
        self.sampleentries=[]
        import copy
        self.dataclass=copy.deepcopy(dataclass)
        td=self.dataclass
        ##produce weighter from a larger dataset as one file
        
        if redo_meansandweights and (td.remove or td.weight):
            logging.info('producing weights and remove indices')
            self.weighter = td.produceBinWeighter(
                self.originRoots
                )            
            self.weighter.printHistos(outputDir)
        
        if redo_meansandweights:
            logging.info('producing means and norms')
            self.means = td.produceMeansFromRootFile(
                self.originRoots, limit=self.meansnormslimit
                )
        
        if means_only: return
        self.__writeData_async_andCollect(0,outputDir)
        
        
    
    def __writeData(self,sample,means, weighter,outputDir,dataclass):
        import os
        import copy
        from stopwatch import stopwatch
        sw=stopwatch()
        td=copy.deepcopy(dataclass)
        
        fileTimeOut(sample,120) #once available copy to ram
        ramdisksample= '/dev/shm/'+str(os.getpid())+os.path.basename(sample)
        
        def removefile():
            os.system('rm -f '+ramdisksample)
        
        import atexit
        atexit.register(removefile)
        
        os.system('cp '+sample+' '+ramdisksample)
        try:
            td.readFromRootFile(ramdisksample,means, weighter) 
            newname=os.path.basename(sample).rsplit('.', 1)[0]
            newpath=os.path.abspath(outputDir+newname+'.z')
            td.writeOut(newpath)
            print('converted and written '+newname+'.z in ',sw.getAndReset(),' sec')
            self.samples.append(newname+'.z')
            self.nsamples+=td.nsamples
            self.sampleentries.append(td.nsamples)
            td.clear()
            self.writeToFile(outputDir+'/snapshot.dc')
        except Exception as e:
            removefile()
            raise e
        removefile()
        
        
    def __writeData_async_andCollect(self, startindex, outputDir):
        
        from multiprocessing import Process, Queue, cpu_count
        wo_queue = Queue()
        import os
        thispid=str(os.getpid())
        if not os.path.isfile(outputDir+'/snapshot.dc'):
            self.writeToFile(outputDir+'/snapshot.dc')
        
        tempstoragepath='/dev/shm/'+thispid
        
        print('creating dir '+tempstoragepath)
        os.system('mkdir -p '+tempstoragepath)
        
        def writeData_async(index,woq):
            
            import copy
            from stopwatch import stopwatch
            sw=stopwatch()
            td=copy.deepcopy(self.dataclass)
            sample=self.originRoots[index]
            fileTimeOut(sample,120) #once available copy to ram
            ramdisksample= tempstoragepath+'/'+str(os.getpid())+os.path.basename(sample)
            
            def removefile():
                os.system('rm -f '+ramdisksample)
            
            import atexit
            atexit.register(removefile)
            success=False
            out_samplename=''
            out_sampleentries=0
            newname=os.path.basename(sample).rsplit('.', 1)[0]
            newpath=os.path.abspath(outputDir+newname+'.z')
            
            
            
            try:
                os.system('cp '+sample+' '+ramdisksample)
                td.readFromRootFile(ramdisksample,self.means, self.weighter) 
                td.writeOut(newpath)
                print('converted and written '+newname+'.z in ',sw.getAndReset(),' sec -', index)
                
                out_samplename=newname+'.z'
                out_sampleentries=td.nsamples
                success=True
                td.clear()
                removefile()
                woq.put((index,[success,out_samplename,out_sampleentries]))
                
                
            except:
                print('problem in '+newname)
                removefile()
                woq.put((index,[False,out_samplename,out_sampleentries]))
                raise 
            
            
        
        
        def __collectWriteInfo(successful,samplename,sampleentries,outputDir):
            if not successful:
                raise Exception("write not successful, stopping")
            
            self.samples.append(samplename)
            self.nsamples+=sampleentries
            self.sampleentries.append(sampleentries)
            self.writeToFile(outputDir+'/snapshot.dc')
            
        processes=[]
        for i in range(startindex,len(self.originRoots)):
            processes.append(Process(target=writeData_async, args=(i,wo_queue) ) )
        
        nchilds = int(cpu_count()/2)-2 if self.nprocs <= 0 else self.nprocs
        #import os
        #if 'nvidiagtx1080' in os.getenv('HOSTNAME'):
        #    nchilds=cpu_count()-5
        if nchilds<1: 
            nchilds=1
        
        #nchilds=10
        
        index=0
        alldone=False
        try:
            while not alldone:
                if index+nchilds >= len(processes):
                    nchilds=len(processes)-index
                    alldone=True
                
                
                logging.info('starting %d child processes' % nchilds)
                for i in range(nchilds):
                    logging.info('starting %s...' % self.originRoots[startindex+i+index])
                    processes[i+index].start()
                        
                results=[]
                import time
                time.sleep(1)

                while 1:
                    running = len(results)<nchilds #  any(p.is_alive() for p in processes)
                    while not wo_queue.empty():
                        res=wo_queue.get()
                        results.append(res)
                        logging.info('collected result %d, %d left' % (res[0], (nchilds-len(results))))
                    if not running:
                        break
                    time.sleep(0.1)
                
                logging.info('joining')
                for i in range(nchilds):
                    processes[i+index].join(5)
                    
                results.sort()
                results = [r[1] for r in results]
                for i in range(nchilds):
                    logging.info(results[i])
                    __collectWriteInfo(results[i][0],results[i][1],results[i][2],outputDir)
                
                index+=nchilds
        
        except:
            os.system('rm -rf '+tempstoragepath)
            raise 
        os.system('rm -rf '+tempstoragepath)
        
    def convertListOfRootFiles(
                    self, inputfile, dataclass, outputDir, 
                    takemeansfrom='', means_only = False,
                    output_name = 'dataCollection.dc', batch_mode = False):
        newmeans=True
        if takemeansfrom:
            self.readFromFile(takemeansfrom)
            newmeans=False
        self.readRootListFromFile(inputfile)
        self.createDataFromRoot(
                    dataclass, outputDir, 
                    newmeans, means_only = means_only, 
                    dir_check= not batch_mode
                    )
        self.writeToFile(outputDir+'/'+output_name)
        
    def getAllLabels(self):
        return self.__stackData(self.dataclass,'y')
    
    def getAllFeatures(self):
        return self.__stackData(self.dataclass,'x')
        
    def getAllWeights(self):
        return self.__stackData(self.dataclass,'w')
    
    
    def getSamplePath(self,samplefile):
        #for backward compatibility
        if samplefile[0] == '/':
            return samplefile
        return self.dataDir+'/'+samplefile
    
    def __stackData(self, dataclass, selector):
        import numpy
        td=dataclass
        out=[]
        firstcall=True
        for sample in self.samples:
            td.readIn(self.getSamplePath(sample))
            #make this generic
            thislist=[]
            if selector == 'x':
                thislist=td.x
            if selector == 'y':
                thislist=td.y
            if selector == 'w':
                thislist=td.w
               
            if firstcall:
                out=thislist
                firstcall=False
            else:
                for i in range(0,len(thislist)):
                    if selector == 'w':
                        out[i] = numpy.append(out[i],thislist[i])
                    else:
                        out[i] = numpy.vstack((out[i],thislist[i]))
                
        return out
    
        
    
        
    def generator(self):
        import numpy
        import copy
        from sklearn.utils import shuffle
        
        #helper class
        class tdreader(object):
            def __init__(self,filelist,maxopen,tdclass):
                
                #print('init reader for '+str(len(filelist))+' files:')
                #print(filelist)
                
                self.filelist=filelist
                self.max=maxopen
                self.nfiles=len(filelist)
                self.tdlist=[]
                self.tdopen=[]
                self.tdclass=copy.deepcopy(tdclass)
                self.tdclass.clear()#only use the format, no data
                for i in range(maxopen):
                    self.tdlist.append(copy.deepcopy(tdclass))
                    self.tdopen.append(False)
                    
                self.closeAll() #reset state
                
            def start(self):
                for i in range(self.max):
                    self.__readNext()
                
            def __readNext(self):
                import copy
                readfilename=self.filelist[self.filecounter]
                self.tdlist[self.nextcounter]=copy.deepcopy(self.tdclass)
                self.tdlist[self.nextcounter].readIn_async(readfilename)
                
                #print('reading file '+readfilename)#DEBUG
                
                self.tdopen[self.nextcounter]=True
                self.filecounter=self.__increment(self.filecounter,self.nfiles)
                
                #if self.filecounter==0:
                #    print('file counter reset to 0 after file '+readfilename)
                
                self.nextcounter=self.__increment(self.nextcounter,self.max)
                
            def __getLast(self):
                td=self.tdlist[self.lastcounter]
                td.readIn_join()
                
                #print('got '+td.samplename+' with '+str(len(td.x[0]))+' samples')
                
                self.tdopen[self.lastcounter]=False
                self.lastcounter=self.__increment(self.lastcounter,self.max)
                return td
                
            def __increment(self,counter,maxval):
                counter+=1
                if counter>=maxval:
                    counter=0   
                return counter 
            def NOT__del__(self):
                #print('del')
                self.closeAll()
                
            def closeAll(self):
                #close all
                for i in range(self.max):
                    if self.tdopen[i]:
                        self.tdlist[i].readIn_abort()
                        self.tdlist[i].clear()
                        self.tdopen[i]=False
                
                self.nextcounter=0
                self.lastcounter=0
                self.filecounter=0
                #print('closed')
                
            def get(self):
                #returns last reads next circular
                td=self.__getLast()
                self.__readNext()
                return td
                
        
        td=(self.dataclass)
        totalbatches=self.getNBatchesPerEpoch()
        processedbatches=0
        
        #print(totalbatches,self.__batchsize,self.nsamples)
        
        xstored=[numpy.array([])]
        dimx=0
        ystored=[]
        dimy=0
        wstored=[]
        dimw=0
        nextfiletoread=0
        
        xout=[]
        yout=[]
        wout=[]
        samplefilecounter=0
        
        #prepare file list
        filelist=[]
        for s in self.samples:
            filelist.append(self.getSamplePath(s))
        
        TDReader=tdreader(filelist, self.maxFilesOpen, self.dataclass)
        
        #print('generator: total batches '+str(totalbatches))
        
        TDReader.start()
        #### 
        #
        # make block class for file read with get function that starts the next read automatically
        # and closes all files in destructor?
        #
        #  check if really the right ones are read....
        #
        psamples=0 #for random shuffling
        while 1:
            if processedbatches == totalbatches:
                processedbatches=0
            
            lastbatchrest=0
            if processedbatches == 0: #reset buffer and start new
                #print('DataCollection: new turnaround')
                xstored=[numpy.array([])]
                dimx=0
                ystored=[]
                dimy=0
                wstored=[]
                dimw=0
                lastbatchrest=0
                
                
            else:
                lastbatchrest=xstored[0].shape[0]
            
            batchcomplete=False
            
            
            
            if lastbatchrest >= self.__batchsize:
                batchcomplete = True
                
            # if(xstored[1].ndim==1):
                
            while not batchcomplete:
                import sys, traceback
                try:
                    td=TDReader.get()
                except:
                    traceback.print_exc(file=sys.stdout)
                
                if xstored[0].shape[0] ==0:
                    #print('dc:read direct') #DEBUG
                    xstored=td.x
                    dimx=len(xstored)
                    ystored=td.y
                    dimy=len(ystored)
                    wstored=td.w
                    dimw=len(wstored)
                    if not self.useweights:
                        dimw=0
                    xout=[]
                    yout=[]
                    wout=[]
                    for i in range(0,dimx):
                        xout.append([])
                    for i in range(0,dimy):
                        yout.append([])
                    for i in range(0,dimw):
                        wout.append([])
                        
                else:
                    if td.x[0].shape == 0:
                        print('Found empty (corrupted?) file, skipping')
                        continue
                    #print('dc:append read sample') #DEBUG
                    for i in range(0,dimx):
                        if(xstored[i].ndim==1):
                            xstored[i] = numpy.append(xstored[i],td.x[i])
                        else:
                            xstored[i] = numpy.vstack((xstored[i],td.x[i]))
                    
                    for i in range(0,dimy):
                        if(ystored[i].ndim==1):
                            ystored[i] = numpy.append(ystored[i],td.y[i])
                        else:
                            ystored[i] = numpy.vstack((ystored[i],td.y[i]))
                    
                    for i in range(0,dimw):
                        if(wstored[i].ndim==1):
                            wstored[i] = numpy.append(wstored[i],td.w[i])
                        else:
                            wstored[i] = numpy.vstack((wstored[i],td.w[i]))
                    
                if xstored[0].shape[0] >= self.__batchsize:
                    batchcomplete = True
                    
                    #random shuffle each time
                    for i in range(0,dimx):
                        xstored[i]=shuffle(xstored[i], random_state=psamples)
                    for i in range(0,dimy):
                        ystored[i]=shuffle(ystored[i], random_state=psamples)
                    for i in range(0,dimw):
                        wstored[i]=shuffle(wstored[i], random_state=psamples)
                    
                    
                    #randomize elements
                     
                #limit of the random generator number 
                psamples+=  td.x[0].shape[0]   
                if psamples > 4e8:
                    psamples/=1e6
                    psamples=int(psamples)
                
                td.clear()


                
            if batchcomplete:
                
                #print('batch complete, split')#DEBUG
                
                for i in range(0,dimx):
                    splitted=numpy.split(xstored[i],[self.__batchsize])
                    xstored[i] = splitted[1]
                    xout[i] = splitted[0]
                for i in range(0,dimy):
                    splitted=numpy.split(ystored[i],[self.__batchsize])
                    ystored[i] = splitted[1]
                    yout[i] = splitted[0]
                for i in range(0,dimw):
                    splitted=numpy.split(wstored[i],[self.__batchsize])
                    wstored[i] = splitted[1]
                    wout[i] = splitted[0]
            
            for i in range(0,dimx):
                if(xout[i].ndim==1):
                    xout[i]=(xout[i].reshape(xout[i].shape[0],1)) 
                if not xout[i].shape[1] >0:
                    raise Exception('serious problem with the output shapes!!')
                            
            for i in range(0,dimy):
                if(yout[i].ndim==1):
                    yout[i]=(yout[i].reshape(yout[i].shape[0],1))
                if not yout[i].shape[1] >0:
                    raise Exception('serious problem with the output shapes!!')
                    
            for i in range(0,dimw):
                if(wout[i].ndim==1):
                    wout[i]=(wout[i].reshape(wout[i].shape[0],1))
                if not xout[i].shape[1] >0:
                    raise Exception('serious problem with the output shapes!!')
            
            processedbatches+=1
            #print('generator: processed batch '+str(processedbatches)+' of '+str(totalbatches)+' '+str(psamples))
            
            #safety check
            
            #xout[1]=numpy.zeros(self.__batchsize) #DEBUG
            #xout[1]=(xout[1].reshape(xout[1].shape[0],1))
            
            #yout[1]=numpy.zeros(self.__batchsize) #DEBUG
            #yout[1]=(yout[1].reshape(yout[1].shape[0],1))
            
            if self.useweights:
                yield (xout,yout,wout)
            else:
                yield (xout,yout)
            
            

    
    
    
    

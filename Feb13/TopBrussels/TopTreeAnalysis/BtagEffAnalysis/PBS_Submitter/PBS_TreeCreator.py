# python script to run the KinFits (TopTreeAnalysis) on PBS.

# for command line options
from optparse import OptionParser

# regular expressions
import re

# interacting with the os
from subprocess import Popen, PIPE, STDOUT
import sys
import os, os.path

# working with time
import time
from time import strftime, gmtime
from datetime import datetime

# import packages for multi-threading
import Queue
import threading

##################
### LogHandler ###
##################

class logHandler:

    def __init__ (self,fileName):

	self.logFile = fileName

        if not self.logFile == "" and os.path.exists(self.logFile):

            os.remove(self.logFile)

    def output(self,string):
        
        if not self.logFile == "":

            f = open(self.logFile,"a")

            f.write("\n["+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"] "+string)

            f.write("\n")

	    f.close()

	else:

      			print "\n["+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"] "+string


###################
### MailHandler ###
###################

# importing smtp lib
import smtplib

class MailHandler:

    def __init__(self,recepient):

        self.smtpServer = "mach.vub.ac.be"
        #self.smtpServer = "localhost"

        self.senderAddress = "PBSTreeCreator@mtop.iihe.ac.be"

        #+Popen('hostname', shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().strip()

        #self.toAnnounce = [ "top-brussels-datasets@cern.ch" ]

        self.toAnnounce = recepient.split(',')
            
    def sendMail(self,subject,msg):

        toAddrs = ""
        
        for to in range(0,len(self.toAnnounce)):
            toAddrs = toAddrs+self.toAnnounce[to]+", "

        m = "From: %s\r\nTo: %s\r\nSubject: %s\r\nX-Mailer: My-Mail\r\n\r\n" % (self.senderAddress, toAddrs, subject)
        
        server = smtplib.SMTP(self.smtpServer)
        server.sendmail(self.senderAddress, toAddrs.split(), m+msg)
        server.quit()


###################
## JobHandler ##
###################

class JobHandler:

    def __init__(self, nJob, name, systematic, crossSection, intLumi, pathPNFS,time,jobNum,start,end):

        self.nJob = nJob
        self.pbsFile = ""
        self.pbsLog = ""
        self.pbsID = ""
				
        self.inputName = name
        self.systematicOption = systematic
        self.inputXS = crossSection
        self.inputLumi = intLumi
        self.inputPNFSDir = pathPNFS
        
        self.xmlCFG = ""
        self.workingDir=""
        self.resultsDir=""
        self.log = ""

        self.walltime=time

        self.taskName = "" 

        self.srmList = "" # container for files needed to be srmcp'd on the WN when using srmcp

        self.start = start

        self.end = end

        self.jobNum = jobNum

    def setlog(self,log):

        self.log = log

    def setupWorkingDir (self):

        global options
        global timestamp
        global userName

        self.resultsDir = "./Results/RESULTS_"+options.TaskName+"_"+timestamp+"/"

        if not os.path.exists(self.resultsDir) and self.nJob == 0:
            os.mkdir(self.resultsDir)

        if not options.local:
            self.workingDir = "/localgrid/"+userName+"/"+options.TaskName+"_"+timestamp+"_job_"+str(self.nJob)
        else:
            if not os.path.exists("./WorkArea"):
                os.mkdir("./WorkArea")
                
            self.workingDir = "./WorkArea/"+options.TaskName+"_"+timestamp+"_job_"+str(self.nJob)

        os.mkdir(self.workingDir)

        self.log.output(Popen("cp -vfr "+options.WorkingDir+"/BTagTReeCreator "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr "+options.WorkingDir+"/*.xml "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr "+options.WorkingDir+"/BtagMassPlots.root "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr "+options.WorkingDir+"/ReweightHistos "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr "+options.WorkingDir+"/PileUpReweighting "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr "+options.WorkingDir+"/JECFiles "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr /user/mmaes/lib/* "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        self.log.output(Popen("cp -vfr ./DownloadSample.sh "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        
        if not options.local and self.inputPNFSDir.find("/user/") == 0 and not self.inputPNFSDir.find("/pnfs/") == 0:
            split = self.srmList.split(",")
            for f in split:
                self.log.output(Popen("cp -vfr "+self.inputPNFSDir+"/"+f+" "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

        #self.log.output(Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
        #sys.exit(1)
        Popen("cp $X509_USER_PROXY "+self.workingDir+"/gridProxy", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

        if os.path.exists(self.xmlCFG):
            os.remove(self.xmlCFG)
            
    def createXMLCFG(self):

        global options
        
        self.log.output("Checking command line option:  "+self.systematicOption)

        self.log.output("Creating xml config files for TopTreeAnalysis")
        
        self.xmlCFG = options.WorkingDir+"/myBTAGconfig_"+str(self.nJob)+".xml"

        xmlFile = open(self.xmlCFG,"w")
        xmlFile.write("<?xml version=\"1.0\"?>\n")
        xmlFile.write("<datasets>\n")

        dataSetXML = "<d name=\""+self.inputName+"\" title=\""+self.inputName+"\" add=\"1\" color=\"1\" ls=\"1\" lw=\"2\" normf=\"1\" xsection=\""+self.inputXS+"\" EqLumi=\""+self.inputLumi+"\" filenames=\""

        fname = "job_"+str(self.nJob)+"_"+str(strftime("%d%m%Y_%H%M%S"))
        file=open(fname,"w")
        cmd = "ls "+self.inputPNFSDir+"/*.root"
        InputFiles = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        file.write(InputFiles)
        file.write("EOF")
        file.close()

        #print self.start
        #print self.end
        iFile = 0
        for file in open(fname):

            #print iFile
            if not file.rfind("EOF") == -1 or (iFile >= self.start and iFile < self.end):

                if options.usesrmcp or ((self.inputPNFSDir.find("/home/mmaes/samples/") == 0 or self.inputPNFSDir.find("/user/") == 0) and not self.inputPNFSDir.find("/pnfs/") == 0):
                    if not options.local or (options.local and options.usesrmcp and self.inputPNFSDir.find("/pnfs/") == 0) or self.inputPNFSDir.find("/home/mmaes/samples/") == 0:
                        split = file.split("/")
                        dataSetXML+=split[len(split)-1].strip()+","                        
                        self.srmList+=split[len(split)-1].strip()+","
                    else:
                        options.usesrmcp=bool(False)
                        split = file.split("/")
                        dataSetXML+=self.inputPNFSDir+"/"+split[len(split)-1].strip()+","
                        #break;
                else:
                    dataSetXML+="dcap://maite.iihe.ac.be"+file.strip()+","

            iFile = iFile+1
            
        os.remove(fname)


        if options.usesrmcp or ((self.inputPNFSDir.find("/home/mmaes/samples/") == 0 or self.inputPNFSDir.find("/user/") == 0) and not self.inputPNFSDir.find("/pnfs/") == 0):            
            dataSetXML = dataSetXML.split(",EOF,")[0]
            self.srmList = self.srmList.split(",EOF,")[0]
        else:
            dataSetXML = dataSetXML.split(",dcap://maite.iihe.ac.beEOF,")[0]
       
        dataSetXML += "\"/>"

        #print dataSetXML

        #sys.exit(1)


        self.log.output(dataSetXML)
            
        xmlFile.write(dataSetXML+"\n")

        xmlFile.write("</datasets>\n\n")
        xmlFile.write("<analysis>\n")
        if len(self.inputName.split("Data")) > 1 or len(self.inputName.split("data")) > 1 or len(self.inputName.split("DATA")) > 1:
            xmlFile.write("<a type=\"Collections\" PVCollection=\"PrimaryVertex\" JetType=\"2\" JetCollection=\"PFJets_selectedPatJetsPF2PAT\" METType=\"2\" METCollection=\"PFMET_patType1CorrectedPFMetPF2PAT\" MuonCollection=\"Muons_selectedPatMuonsPF2PAT\" ElectronCollection=\"Electrons_selectedPatElectronsPF2PAT\" loadGenJetCollection=\"1\" GenJetCollection=\"GenJets_ak5GenJetsNoNu\" loadGenEventCollection=\"1\" GenEventCollection=\"GenEvent\" loadNPGenEventCollection=\"0\" NPGenEventCollection=\"NPGenEvent\" loadMCParticles=\"1\" MCParticlesCollection=\"MCParticles\" TrackMETCollection=\"\" loadTrackMET=\"0\"/>")
        else:
            xmlFile.write("<a type=\"Collections\" PVCollection=\"PrimaryVertex\" JetType=\"2\" JetCollection=\"PFJets_selectedPatJetsPF2PAT\" METType=\"2\" METCollection=\"PFMET_patType1CorrectedPFMetPF2PAT\" MuonCollection=\"Muons_selectedPatMuonsPF2PAT\" ElectronCollection=\"Electrons_selectedPatElectronsPF2PAT\" loadGenJetCollection=\"1\" GenJetCollection=\"GenJets_ak5GenJetsNoNu\" loadGenEventCollection=\"1\" GenEventCollection=\"GenEvent\" loadNPGenEventCollection=\"0\" NPGenEventCollection=\"NPGenEvent\" loadMCParticles=\"1\" MCParticlesCollection=\"MCParticles\" TrackMETCollection=\"\" loadTrackMET=\"0\"/>")
            
        xmlFile.write("<a type=\"Selection\" PVertexNdofCut=\"4\" PVertexZCut=\"24.\" PVertexRhoCut=\"2.\" MuonPtCutSR=\"20.\" MuonEtaCutSR=\"2.1\" MuonRelIsoCutSR=\"0.05\" MuonNHitsCutSR=\"10\" MuonD0CutSR=\"0.02\" MuonDRJetsCut=\"0.3\" MuonPtCutVetoSR=\"10.\" MuonEtaCutVetoSR=\"2.5\" MuonRelIsoCutVetoSR=\"0.2\" ElectronPtCut=\"15.\" ElectronEtaCut=\"2.5\" ElectronRelIsoCut=\"0.2\" JetsPtCutSR=\"30.\" JetsEtaCutSR=\"2.4\" applyJetID=\"1\" JetEMFCut=\"0.01\" n90HitsCut=\"1\" fHPDCut=\"0.98\" NofJets=\"4\" NofJetBins=\"2\"/>\n")
        xmlFile.write("<a type=\"Conditions\" isMC=\"1\" MCRound=\"0\" Vars_ByFile=\"0\" VarsFile=\"m0_100_m12_100\" IntToCut=\"4\" Verbose=\"2\" Luminosity=\"9999999\" JES=\"1.\" nPseudoExp=\"0\" nPseudoSession=\"0\" runonTTrees=\"0\" doABCD=\"1\" doVJEstim=\"1\" doVJEstPE=\"1\" doTtJEstim=\"1\" doTemplComp=\"0\" doSystematics=\"0\"/>\n")
        xmlFile.write("<a type=\"CRForTtbarEstimation\" BtagAlgo_ttjEst=\"0\" BtagDiscriCut_ttjEst=\"4.38\" MuonPtCutCR=\"30.\" MuonEtaCutCR=\"2.1\" MuonRelIsoCutCR=\"0.1\" JetsPtCutCR=\"30.\" JetsEtaCutCR=\"2.4\" MblCut=\"160.\" DRBBCut=\"2.3\" HTBBCut=\"500.\" NREvtFraction=\"0.75\"/>\n")
        xmlFile.write("<a type=\"CRForABCDEstimation\" NXbinsABCD=\"200\" NYbinsABCD=\"200\" XbinMinABCD=\"0\" XbinMaxABCD=\"20\" YbinMinABCD=\"0\" YbinMaxABCD=\"20\" cutXmin=\"0\" cutX0=\"0.1\" cutX1=\"0.2\" cutXmax=\"20.\" cutYmin=\"0.\" cutY0=\"3.\" cutY1=\"4.\" cutYmax=\"20.\" region=\"1\"/>\n")
        xmlFile.write("<a type=\"ParamForVJetEstimation\" BtagAlgo_vjEst=\"0\" NofBtagWorkingPoint_vjEst=\"1\" BtagWorkingPoint_vjEst=\"2.03,3.20\" MinMethod=\"Minuit2\" MinOption=\"Combined\" useMJLE=\"0\" useUnBinMLE=\"1\" NVJetPE=\"500\" TagEffInit=\"0.794,0.128,0.097-0.70,0.043,0.02-0.63,0.05,0.010/0.807,0.134,0.124-0.70,0.043,0.02-0.63,0.05,0.010\" NVlikeInit=\"14./4.\" NTTlikeInit=\"6./8.\" EffEbsel=\"0.0515,0.4170,0.5281/0.0187,0.2604,0.7049\" VJEstFixParam=\"0,1,2\" NofIterationsVJestShapeEstim=\"40\"/>\n")
        xmlFile.write("<a type=\"Observables\" runOnObsByString=\"0\" listOfObsInts=\"2\" listOfObsStrings=\"ET1oET2,ET1oET3\" binning=\"../config/Binning.root\" bins=\"20\"/>\n")
        xmlFile.write("<a type=\"CrossSection\" MCExpFilename=\"../config/MCFile.root\" LuminosityError=\"0.1\" TriggerEff=\"1\" TriggerEffError=\"0.05\" SkimEff=\"1.\" SkimEffError=\"0\" MuonSelEff=\"0.43\" MuonSelEffError=\"0.003\" SecondLeptonVetoEff=\"0.4833\" SecondLeptonVetoEffError=\"0.01\" JetSelEff=\"0.8206\" JetSelEffError=\"0.04\" NofSingleTopEvts=\"0\" NofSingleTopEvtsError=\"0\"/>\n")
        xmlFile.write("<a type=\"Search\" doBkgEstim=\"1\" doDumpPseudoExpInfoInTTree=\"1\" DumpTreeName=\"dumpTreeFile.root\" FractionHWEvts=\"0.1,0.2,0.3\"/>\n")
        xmlFile.write("</analysis>")
        xmlFile.close()
        
    def createPBSCFG (self):

        global timestamp
        global userName
        global options

        self.pbsFile = "/localgrid/"+userName+"/"+options.TaskName+"_"+self.inputName+"_"+timestamp+"_job_"+str(self.nJob)+".pbs"
        self.pbsLog = "/localgrid/"+userName+"/"+options.TaskName+"_"+self.inputName+"_"+timestamp+"_job_"+str(self.nJob)+".log"

        if options.local:
            self.pbsLog=""
            self.pbsFile=self.workingDir+"/runme.sh"

        self.taskName = self.inputName+"_job"+str(self.nJob)+"_"+timestamp

        pbs = open(self.pbsFile,"w")

        if not options.local:
            pbs.write("#! /bin/bash\n")
            pbs.write("#PBS -l walltime="+self.walltime+"\n")
            pbs.write("#PBS -r n\n")
            pbs.write("#PBS -N "+self.taskName+"\n")
            pbs.write("#PBS -j oe\n")
            pbs.write("#PBS -k oe\n")
        
            pbs.write("echo dumping some info on the worker node\n")
            pbs.write("hostname\n")
            pbs.write("df -h\n")
            pbs.write("uptime\n")
            pbs.write("free\n")
            pbs.write("ls -l /scratch/\n")
            pbs.write("export X509_USER_PROXY=/scratch/$PBS_JOBID/PBSJOB/gridProxy\n")
            pbs.write("export DCACHE_RA_BUFFER=\"250000000\"\n")
            
        pbs.write("export ROOTSYS=/localgrid/"+userName+"/root\n")
        pbs.write("export PATH=$ROOTSYS/bin:$PATH\n")
        pbs.write("export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH\n")
        pbs.write("export LD_LIBRARY_PATH=/scratch/$PBS_JOBID/PBSJOB:$LD_LIBRARY_PATH\n")
        pbs.write("echo \"Root Version: $(root-config --version)\"\n")

        if not options.local:
            pbs.write("echo moving workingdir to WN scratch\n")
            pbs.write("\nmv -vf "+self.workingDir+" /scratch/$PBS_JOBID/PBSJOB\n")
            pbs.write("cd /scratch/$PBS_JOBID/PBSJOB\n")    
        
        if not options.usesrmcp:
            if not options.local:
                pbs.write("source /jefmount_mnt/jefmount/cmss/slc5_amd64_gcc434/external/dcap/2*/etc/profile.d/init.sh\n") # to fix the fact that the WN's dont have 64bit libdcap

        elif self.inputPNFSDir.find("/home/mmaes/samples/") == 0:

            # download files via curl

            InputFiles = self.srmList.split(',')

            for i in InputFiles:
                pbs.write("echo \"Downloading file:  "+str(i)+"\"\n")

                url = "https://mtop.wn.iihe.ac.be/samples/"+self.inputPNFSDir.split('/')[4]+"/"+str(i)

                if not options.local:
                    #pbs.write("wget --no-check-certificate -O /scratch/$PBS_JOBID/PBSJOB/"+str(i)+" "+url+"\n")
                    pbs.write("curl -k -o /scratch/$PBS_JOBID/PBSJOB/"+str(i)+" "+url+" 2>&1\n")
                else:
                    #pbs.write("wget --no-check-certificate -O ./"+str(i)+" "+url+" 2> output\n")
                    pbs.write("curl -k -o ./"+str(i)+" "+url+" 2>&1\n")
            
        elif not self.inputPNFSDir.find("/user/") == 0:

            pbs.write("sh DownloadSample.sh "+self.inputPNFSDir+" "+self.srmList+"\n")
            
        #sys.exit(1)
        
        pbs.write("ls -ltr\n\n")

        #pbs.write("exit\n\n")
        
        pbs.write("./BTagTReeCreator "+str(self.systematicOption)+" 0 myBTAGconfig_"+str(self.nJob)+".xml\n")

        pbs.write("ls -l BtagTrees/\n")

        if not options.local:
            if options.usesrmcp or (self.inputPNFSDir.find("/user/") == 0 and not self.inputPNFSDir.find("/pnfs/") == 0):
                pbs.write("echo removing locally copied TopTree files before stageout\n")
                split = self.srmList.split(",")
                for f in split:
                    pbs.write("rm -v "+f+"\n")

        if not options.local:
            pbs.write("echo Moving working directory back to localgrid\n")
            pbs.write("\nmv -f /scratch/$PBS_JOBID/PBSJOB "+self.workingDir+"\n")
        pbs.write("echo \"THIS IS THE END\"\n")
        
        pbs.close()
        
    def submitPBSJob (self):

        global timestamp
        global userName

        cmd = "cd /localgrid/"+userName+"/; qsub -q localgrid@cream02.iihe.ac.be "+self.pbsFile

        #print cmd

        #sys.exit(1)

        self.pbsID=Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        
        self.log.output("PBS Job ID: "+self.pbsID)

    def checkPBSJob (self):
        
        global timestamp
        global userName
        global SleepTime

        status = Popen("qstat "+self.pbsID, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        
        while status.find("Unknown Job Id") == -1:

            self.log.output("     -----> It seems that the job is still running, sleeping "+str(SleepTime)+"s")
            
            time.sleep(SleepTime)

            status = Popen("qstat "+self.pbsID, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

        self.log.output("-> JOB Finished, dumping output")
        self.pbsLog = (Popen("echo /localgrid/"+userName+"/"+self.taskName+".o*", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()).strip()
        if os.path.exists(self.pbsLog):
            for line in open(self.pbsLog):
                self.log.output(line.strip())

    def runLocalJob (self):

        global options

        cmd = "cd "+self.workingDir+"; source runme.sh >> output"
        cmdtail = "tail -f "+self.workingDir+"/output"

        self.log.output("Starting job with command: "+cmd)

        pExe = Popen("cd "+self.workingDir+"; cat runme.sh", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

        self.log.output(pExe.stdout.read())
        
        self.log.output("Progress can be tracked with : "+cmdtail)

        pExe = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

        exit=pExe.poll()
        while exit == None:
            exit = pExe.poll()
            time.sleep(30)

        f=open(self.workingDir+"/done","w")
        f.write(pExe.stdout.read())
        f.close()
        
        self.log.output("Job has Finished")

        for line in open(self.workingDir+"/output","r"):
             self.log.output(line.strip())

        #print exit
        
        if not exit == 0:

            self.log.output("  --> Seems that the job has failed!!!!!!!")

    def process(self):

        if not options.local:
						
            self.createXMLCFG()
            
            self.setupWorkingDir()

            #sys.exit(1)

            self.createPBSCFG()

            #sys.exit(1)
            
            self.submitPBSJob()

            #sys.exit(1)

            self.checkPBSJob()

            if self.jobNum == -1:

                self.log.output(Popen("cp -f "+self.workingDir+"/BtagTrees/*.root "+self.resultsDir+"; cp -f "+self.workingDir+"/BtagTrees/*.txt "+self.resultsDir+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

            else:

                if not os.path.exists(self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum)):
                    os.mkdir(self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum))

                self.log.output(Popen("cp -f "+self.workingDir+"/BtagTrees/*.root "+self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum)+"; cp -f "+self.workingDir+"/BtagTrees/*.txt "+self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum)+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

            # clean up

            if os.path.exists(self.pbsFile):
                os.remove(self.pbsFile)
            if os.path.exists(self.pbsLog):
                os.remove(self.pbsLog)
                
        else:

            self.createXMLCFG()
            
            self.setupWorkingDir()

            self.createPBSCFG()

            self.runLocalJob()

            if self.jobNum == -1:

                self.log.output(Popen("cp -f "+self.workingDir+"/BtagTrees/*.root "+self.resultsDir+"; cp -f "+self.workingDir+"/BtagTrees/*.txt "+self.resultsDir+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

            else:

                if not os.path.exists(self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum)):
                    os.mkdir(self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum))

                self.log.output(Popen("cp -f "+self.workingDir+"/BtagTrees/*.root "+self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum)+"; cp -f "+self.workingDir+"/BtagTrees/*.txt "+self.resultsDir+"/"+self.inputName+"_job"+str(self.jobNum)+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

            # clean up

            if os.path.exists(self.pbsFile):
                os.remove(self.pbsFile)
            if os.path.exists(self.pbsLog):
                os.remove(self.pbsLog)
            
        return True

   
##############
## WorkFlow ##
##############

class WorkFlow (threading.Thread ):

    def __init__(self, nThread, *args, **kwds):

        global options
        global log
        global timestamp
        global nJobs

        self.nThread = nThread        

        threading.Thread.__init__(self, *args, **kwds)
    
        self.keepAlive = bool(True)

        if int(nJobs) == 1:
            self.log = log
        else:
            self.log = logHandler("logs/"+options.TaskName+"_"+timestamp+"_thread"+str(self.nThread)+".txt")

    def stop (self):
        
        self.keepAlive = bool(False)

    def run (self):

        global nJobs
        
        # our thread runs forever
        
        while not jobsPool.empty():

            #print self.nThread

            job = jobsPool.get()

            job.setlog(self.log)

            if int(nJobs) == 1:
                log.output("-> Thread "+str(self.nThread)+" Processing Job: "+str(job.nJob+1)+"/"+str(nJobs))
            else:
                log.output("-> Thread "+str(self.nThread)+" Processing Job: "+str(job.nJob+1)+"/"+str(nJobs)+" (LogFile: "+self.log.logFile+")")

            if (job.process()):
                log.output("-> Thread "+str(self.nThread)+" Finished Job: "+str(job.nJob+1)+"/"+str(nJobs))


###############
### METHODS ###
###############

def checkGridProxy ():

    global log

    log.output("* Checking GRID proxy")
        
    cmd = 'voms-proxy-info'
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    
    lines = output.split('\n')
        
    #print lines

    if lines[1] == "Couldn't find a valid proxy.":
        
        log.output(" ---> No Valid GRID proxy was retrieved, please use voms-proxy-init --voms cms:/cms/becms --valid 190:00 before running the skimmer.")
        
        sys.exit(1)
        
    else:

        for i in xrange(len(lines)):
            
            if not lines[i].rfind("timeleft") == -1:
                
                timeleft = lines[i].split("timeleft  : ")
                
                splittedtimeLeft = timeleft[1].split(":");
                
                if int(splittedtimeLeft[0]) < 50:
                    
                    log.output(" ---> GRID proxy is valid < 50h, please renew it using voms-proxy-init --voms cms:/cms/becms --valid 190:00 before running this script!")
                    
                    sys.exit(1)

                else:
                        
                    log.output(" ---> Valid GRID proxy was found (Valid for "+str(timeleft[1])+")")


###############
### OPTIONS ###
###############

optParser = OptionParser()

optParser.add_option("-t","--taskName", dest="TaskName",default="BTagTreeCreator",
                     help="TaskName that will be used as prefix for workingdir and logs", metavar="")
optParser.add_option("-d","--workingDir", dest="WorkingDir",default="",
                     help="Directory containing the needed TopTreeAnalysis setup and the file containing input datasets (inputSamples.txt)", metavar="")
optParser.add_option("-m","--mail", dest="Mail", default="",
                     help="E-mail adress to inform when the script finished", metavar="")
optParser.add_option("-o","--log-stdout", action="store_true", dest="stdout",default=bool(False),
                     help="Write the main log file to the stdout", metavar="")
optParser.add_option("-l","--run-local", action="store_true", dest="local",default=bool(False),
                     help="Use local CPUs", metavar="")
optParser.add_option("","--srmcp", action="store_true", dest="usesrmcp",default=bool(False),
                     help="Use SRMCP on the WNs to download rootfiles first from /pnfs. If false, libDcap is used. THIS IS ONLY FOR SAMPLES ON PNFS", metavar="")
optParser.add_option("-w","--walltime", dest="walltime",default="07:00:00",
                     help="Define the PBS job walltime.", metavar="")
   
(options, args) = optParser.parse_args()

if options.WorkingDir == "":

    print "No working dir provided! For help use python PBS_TreeCreator.py -h"

    sys.exit(1)

############
### MAIN ###
############

## SETTINGS

#RootInstallation = "/home/LocalSoft/root"
#RootInstallation = "/home/LocalSoft/specialofficial_root/root/"
#RootInstallation = "/user/pvmulder/special_rootpatch/v5-30-00-patches"
#RootInstallation = "/user/cmssoft/root_old/"
RootInstallation = "/Software/LocalSoft/root_5.30.02/root"

SleepTime = int(60) # time to sleep between checking of job status
#SleepTime = int(30) # time to sleep between checking of job status

## END SETTINGS

# special dirs

if not os.path.exists("logs"):
    os.mkdir("logs")
if not os.path.exists("Results"):
    os.mkdir("Results")
    
# timestamp

timestamp = strftime("%d%m%Y_%H%M%S") # need a timestamp for dirs and logfiles

# logging

if not options.stdout:
    #log = logHandler("logs/log_"+timestamp+".txt")
    log = logHandler("logs/"+options.TaskName+"_"+timestamp+".txt")
else:
    log = logHandler("")

# start the loop

log.output("*** TopTreeAnalysis BtagTreeCreator Batch Job system ***")

log.output("Making a dump of the configuration used:")
log.output("  --> TaskName:\t\t"+options.TaskName)
log.output("  --> WorkingDir:\t\t"+options.WorkingDir)
log.output("  --> log-stdout:\t\t"+str(options.stdout))
log.output("  --> run-local:\t\t"+str(options.local))
log.output("  --> Use SRMCP instead of libDcap:\t\t"+str(options.usesrmcp))

# get username

userName = Popen('echo $USER', shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().strip()

# check grid proxy and copy it to the WorkingDir if it exists

if not options.local:
    
    checkGridProxy()
    
# copy root to localgrid

if not options.local:
    if not os.path.exists("/localgrid/"+userName+"/root"):

        log.output("-> Copying root to localgrid") 
        
        cmd = "cp -vfr "+RootInstallation+" /localgrid/"+userName+"/root"
        
        Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        
        ver = open("/localgrid/"+userName+"/root/ver","w")
        
        ver.write(RootInstallation)
        
        ver.close()
        
    else:
        
        rootver = ""
        
        for ver in open("/localgrid/"+userName+"/root/ver","r"):
            
            rootver = ver
            
            if not rootver == RootInstallation:
                
                log.output("-> RE-Copying root to localgrid") 

                cmd = "mv /localgrid/"+userName+"/root /localgrid/"+userName+"/root_old_"+timestamp
                cmd += " ;cp -vfr "+RootInstallation+" /localgrid/"+userName+"/root"
                
                Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                
                ver = open("/localgrid/"+userName+"/root/ver","w")
        
                ver.write(RootInstallation)
                
                ver.close()
        

# Read in the inputSamples.txt file and create our Queue to store jobs:

jobsPool = Queue.Queue ( 0 )

nJobs = int(0)

systematics = [0]

for line in open(options.WorkingDir+"/inputSamples.txt"):

    if len(line.split("#")) < 2 and not line.rfind("@SYSTEMATICS:") == -1:
        systematics = ((line.split("@SYSTEMATICS:")[1]).strip()).split(",")

for syst in systematics:
    for line in open(options.WorkingDir+"/inputSamples.txt"):

        if len(line.split("#")) < 2 and line.rfind("@SYSTEMATICS:") == -1:
    
            splitted = line.split(":")
    
            if len(splitted) > 1:

                wall = options.walltime

                if len(splitted) > 6:
                    wall = (splitted[6].strip()).replace("-",":")

                # split samples in different jobs
                nGroupFiles = -1
                
                if len(splitted) > 5:
                    nGroupFiles = int(splitted[5].strip())

                if splitted[1] != "-1": #data -> no syst

                    dir = splitted[4].split("\n")[0]
                    cmd = "ls "+dir+"/*.root | wc -l"
                    nInputFiles = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

                    nInputFiles = float(nInputFiles.strip("\n"))

                    if nGroupFiles == -1: nGroupFiles=int(nInputFiles)

                    nCycles = int(nInputFiles/nGroupFiles)

                    if nCycles*nGroupFiles < nInputFiles:
                        nCycles=nCycles+1

                    #nCycles=1

                    for i in xrange(nCycles):

                        start=i*nGroupFiles
                    
                        end=(i+1)*nGroupFiles

                        jobNum = -1

                        if not nGroupFiles == nInputFiles:
                            jobNum = i

                        if end > nInputFiles:
                            end=int(nInputFiles)
                
                        if i == 0:
                            
                            if nCycles > 1:
                                log.output("-> Adding sample "+str(splitted[0])+" <-> Wall Time: " +str(wall)+" <-> Systematics switch: "+str(splitted[1])+" <-> Systematic: "+str(syst)+" (nJobs: "+str(nCycles)+")")
                            else:
                                log.output("-> Adding sample "+str(splitted[0])+" <-> Wall Time: " +str(wall)+" <-> Systematics switch: "+str(splitted[1])+" <-> Systematic: "+str(syst))

                        job = JobHandler(nJobs, splitted[0], str(syst), splitted[2], splitted[3], splitted[4].split("\n")[0],wall,jobNum,start,end)
                        jobsPool.put(job)
                    
                        nJobs += 1

# extra loop for the samples where no systematics should be run

for line in open(options.WorkingDir+"/inputSamples.txt"):

    if len(line.split("#")) < 2:
    
        splitted = line.split(":")
    
        if len(splitted) > 1:

            wall = options.walltime

            if len(splitted) > 6:
                wall = (splitted[6].strip()).replace("-",":")

            # split samples in different jobs
            nGroupFiles = -1
                
            if len(splitted) > 5:
                nGroupFiles = int(splitted[5].strip())

            if splitted[1] == "-1": #data -> no syst

                dir = splitted[4].split("\n")[0]
                cmd = "ls "+dir+"/*.root | wc -l"
                nInputFiles = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

                nInputFiles = float(nInputFiles.strip("\n"))

                if nGroupFiles == -1: nGroupFiles=int(nInputFiles)

                nCycles = int(nInputFiles/nGroupFiles)

                if nCycles*nGroupFiles < nInputFiles:
                    nCycles=nCycles+1
                        
                for i in xrange(nCycles):

                    start=i*nGroupFiles
                    
                    end=(i+1)*nGroupFiles

                    jobNum = -1

                    if not nGroupFiles == nInputFiles:
                        jobNum = i

                    if end > nInputFiles:
                        end=int(nInputFiles)
                
                    if i == 0:
                            
                        if nCycles > 1:
                            log.output("-> Adding sample "+str(splitted[0])+" <-> Wall Time: " +str(wall)+" <-> No Systematics (nJobs: "+str(nCycles)+")")
                        else:
                            log.output("-> Adding sample "+str(splitted[0])+" <-> Wall Time: " +str(wall)+" <-> No Systematics")

                    job = JobHandler(nJobs, splitted[0], str(0), splitted[2], splitted[3], splitted[4].split("\n")[0],wall,jobNum,start,end)
                    jobsPool.put(job)
                    
                    nJobs += 1

#print nJobs
#sys.exit(1)

# start our threads

workers = []
#for x in xrange ( int(nJobs) ):

nWorkers = int(nJobs)

if nWorkers > 30: nWorkers=30

for x in xrange ( nWorkers ):
   workers.append(WorkFlow(x))
   workers[x].start()

   time.sleep(5)

# check if there are workers still working

notDone=bool(True)

while notDone:

    notDone = False

    for worker in workers:
    
        if worker.isAlive(): # If there is one worker alive, we are still not finished
            
            notDone=bool(True)

    if not notDone:

        log.output("-> All jobs are DONE")

    #time.sleep(10)
		
# send mail when finished
if not options.Mail == "":
    log.output("Sending announcement to:  "+options.Mail)
    mailHandler = MailHandler(options.Mail)
    message = "The script finished ;-)\n\n"+Popen("cat "+options.WorkingDir+"/inputSamples.txt", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    mailHandler.sendMail("[PBS_TreeCreator.py], script finished", message)

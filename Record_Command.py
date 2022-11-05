import multiprocessing
from time import sleep
import os
import sys
import sounddevice as sd
import soundfile as sf



#Define SSH record command for each Raspberry Pi 
def pi1():
    os.system("ssh pi1@192.168.137.240 sudo nice -n -20 arecord  -f S16_LE -r 48000 -d 20 -c 1 -D plughw:1 pi1.wav")
    
def pi2():
    os.system("ssh pi2@192.168.137.71 sudo nice -n -20 arecord  -f S16_LE -r 48000 -d 20 -c 1 -D plughw:1 pi2.wav")

def pi3():
    os.system("ssh pi3@192.168.137.198 sudo nice -n -20 arecord  -f S16_LE -r 48000 -d 20 -c 1 -D plughw:1 pi3.wav")

def pi4():
    os.system("ssh pi4@192.168.137.103 sudo nice -n -20 arecord  -f S16_LE -r 48000 -d 20 -c 1 -D plughw:1 pi4.wav")
    
    
if __name__ == '__main__':
    argumentList = sys.argv[1:]                       #Read in File name
    file_name = argumentList[0] 
    os.system("mkdir "+str(file_name))                #Create directory
    
    #Define parallel processes and target the ssh command
    process1 = multiprocessing.Process(target=pi1)
    process2 = multiprocessing.Process(target=pi2)
    process3 = multiprocessing.Process(target=pi3)
    process4 = multiprocessing.Process(target=pi4)

    #Execute parallel processes to begin recording
    process1.start()
    process2.start()
    process3.start()
    process4.start()
    
    sleep(5)                                           #Wait 5 seconds
    
    #Read in calibration signal
    filename = 'chirp_0_15000_5s.wav' 

    data, fs = sf.read(filename, dtype='float32')      # Extract data and sampling rate from file
    sd.play(data, fs)                                  #Play calibration signal
    status = sd.wait()                                 # Wait until file is done playing
    
    sleep(15)                                          #Wait until recordings fininsh
   

    #Collect recordings from Raspberry Pis
    os.system("scp pi1@192.168.137.240:pi1.wav "+file)
    os.system("scp pi2@192.168.137.71:pi2.wav "+file)
    os.system("scp pi3@192.168.137.198:pi3.wav "+file)
    os.system("scp pi4@192.168.137.103:pi4.wav "+file)
    
    exit
exit
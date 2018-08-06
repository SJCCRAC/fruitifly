'''
Program Description:
1. *Opens dataCollect1.csv
2. *Saves amino-acid-sequences to a dictionary
    2a. Dict format:
        {'1':{'>D1 207.175833 1.000000':['AMINO-ACID-SEQUENCE', 'startBelow', 'endBelow', 'seqBelow', startAbove]}}
3. *Iterates over amino-acid-sequences in dictionary
    3a. Runs sequence through iupred database
    3b. Stores iupred data output as json file to ordered dictionary
    3c. Identifies regions in each sequence above and below 0.5
    3d. Outputs ordered dict with subsets of original amino-acid-sequences
4. Runs disordered regions through Blast
5. Collects subject lines from the first occurance of each fly species (could be other organism type)
6.

Need to add:
1. Function to open dataCollect1.csv > save amino-acid-sequences to dictionary
2. For loop to iterate over dictionary containing amino-acid-sequences
'''


import warnings  # allows program to be able to ignore benign warnings
#####
# IGNORE WARNINGS
#####
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
# chrome headless testing
import os
import time
import requests
import ast  # supposed to help jason -> dict
import collections  # allows orderedDict
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd  # allows for csv r/w
import json
import ijson
from pandas.io.json import json_normalize
import multiprocessing as mp

'''
USEFUL SELENIUM:
-fullscreen_window(self)
    Invokes the window manager-specific 'full screen' operation
-close(self)
    Closes the current window.
    :Usage:
        driver.close()
-back(self)
     |      Goes one step backward in the browser history.
     |
     |      :Usage:
     |          driver.back()

get(self, url)
     |      Loads a web page in the current browser session.

get_window_position(self, windowHandle='current')
     |      Gets the x,y position of the current window.
     |
     |      :Usage:
     |          driver.get_window_position()






'''


class headless:

    # class variables - open to all instances not just individual instance
    def __init__(self, seq):
        self.seq = seq

    def i_scrape(self, q='q'):  # output added for multiprocessing testing
        ######
        # SET OPTIONS
        ######
        chrome_options = Options()
        chrome_options.add_argument("--headless")
        chrome_options.add_argument("--window-size=1920x1080")
        # chrome_options.binary_location = "C:/ProgramFiles/Google/Chrome/chrome.exe" #part of first tutorial not sure what this is for
        chrome_driver = os.getcwd() + "\\chromedriver.exe"

        ######
        # CREATE BROWSER OBJECT > OPEN URL
        ######
        # driver = webdriver.Chrome(executable_path=os.path.abspath("chromedriver"), chrome_options=chrome_options)
        # driver = webdriver.Chrome("C:/Users/SJCCRAC/Downloads/chromedriver_win32/chromedriver.exe")
        driver = webdriver.Chrome(chrome_options=chrome_options, executable_path="C:/Users/SJCCRAC/Downloads/chromedriver_win32/chromedriver.exe")
        driver.implicitly_wait(10)
        driver.get("https://iupred2a.elte.hu/")  # opens iupred2 website

        ######
        # NEED FUNCTION TO GET DATA SET TO SEARCH
        ######

        ######
        # FIND SEARCH FIELD BY ID > CLEAR ANY TEXT > PASTE IN DATASET
        ######
        search_field = driver.find_element_by_id("inp_seq")
        search_field.clear()
        search_field.send_keys(self.seq)

        ######
        # FIND PREDICTION TYPE BUTTON BY XPATH > CLICK SUBMIT
        ######
        predBtn = driver.find_element_by_xpath("/html/body/table/tbody/tr[4]/td[3]/form/table/tbody/tr[@id='tr1'][2]/td[@id='tr1']/fieldset/label[2]/input[@id='a_type']")
        predBtn.click()
        print('predBtn clicked')

        ######
        # FIND SUBJECT BUTTON BY XPATH > CLICK SUBMIT
        ######
        subjectBtn = driver.find_element_by_xpath("/html/body/table/tbody/tr[4]/td[3]/form/table/tbody/tr[@id='tr2']/td[@id='td2']/input[1]")
        subjectBtn.click()
        print('Subject btn clicked')

        ######
        # SAVES NAME/HANDLE OF CURRENT WINDOW (SO WE CAN SWITCH TO NEEDED BROWSER WINDOW)
        ######
        window_before = driver.window_handles[0]

        ######
        # FIND 'DOWNLOAD JSON' BUTTON BY XPATH > CLICK BTN
        ######
        jsonBtn = driver.find_element_by_xpath("/html/body/table/tbody/tr[3]/td[3]/table/tbody/tr/td/form[2]/input[2]")
        jsonBtn.click()

        ######
        # WAIT 2SEC > SAVES NAME OF NEW BROWSER WINDOW HANDLE > SWITCHES FOCUS TO NEW WINDOW
        ######
        time.sleep(2)
        window_after = driver.window_handles[1]
        driver.switch_to_window(window_after)

        ######
        # FIND TEXT BODY BY XPATH > SAVE TO DICT 'IUPRED' > PRINT TO CONSOLE
        ######
        json_data = driver.find_element_by_xpath("/html/body")
        iupred = ast.literal_eval(json_data.text)  # uses ast library to read raw json.. it then sees/saves it as a python dict
        # print(iupred['iupred2']) ##use to console testing

        ######
        # CLOSES ALL WINDOWS PROGRAM OPENED
        ######
        for i in xrange(0, len(driver.window_handles)):
            handle = driver.window_handles[0]
            driver.switch_to_window(handle)
            driver.close()
            print("closed: " + handle)

        # print(iupred)
        print("returned iupred")
        if(q != 'q'):
            q.put(iupred)
            print("q passed")
        return iupred

    def blast_scrape(self, q='q'):  # output added for multiprocessing testing
        ######
        # SET OPTIONS
        ######
        chrome_options = Options()
        # chrome_options.add_argument("--headless")
        chrome_options.add_argument("--window-size=1920x1080")
        # chrome_options.binary_location = "C:/ProgramFiles/Google/Chrome/chrome.exe" #part of first tutorial not sure what this is for
        chrome_driver = os.getcwd() + "\\chromedriver.exe"

        ######
        # CREATE BROWSER OBJECT > OPEN URL
        ######
        # driver = webdriver.Chrome(executable_path=os.path.abspath("chromedriver"), chrome_options=chrome_options)
        # driver = webdriver.Chrome("C:/Users/SJCCRAC/Downloads/chromedriver_win32/chromedriver.exe")
        driver = webdriver.Chrome(chrome_options=chrome_options, executable_path="C:/Users/SJCCRAC/Downloads/chromedriver_win32/chromedriver.exe")
        driver.implicitly_wait(280)
        driver.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome")  # opens iupred2 website

        ######
        # FIND SEARCH FIELD BY ID > CLEAR ANY TEXT > PASTE IN DATASET
        ######
        search_field = driver.find_element_by_id("seq")
        search_field.clear()
        search_field.send_keys(self.seq)

        ######
        # FIND SEARCH FIELD BY ID > CLEAR ANY TEXT > PASTE IN DATASET
        ######
        blastBtn = driver.find_element_by_id("b1")
        blastBtn.click()
        print('blastBtn clicked')

        ######
        # FIND DOWNLOAD BUTTON TO EXPAND BY XPATH
        ######
        dwnLoadBtn1 = driver.find_element_by_xpath("/html/body[@id='type-a']/div[@id='wrap']/div[@id='content-wrap']/div[@id='content']/a[@id='showDownload']/span[@class='ui-ncbitoggler-master-text']")
        dwnLoadBtn1.click()
        print('dwnLoadBtn1 clicked')

        time.sleep(5)

        ######
        # FIND DOWNLOAD BUTTON BY XPATH > CLICK SUBMIT
        ######
        dwnLoadBtn = driver.find_element_by_link_text('Single-file JSON')
        dwnLoadBtn.click()
        print('dwnLoadBtn clicked')

        ######
        # SLEEP FOR TESTING
        ######
        time.sleep(1)

        ######
        # FIND SUBJECT BUTTON BY XPATH > CLICK SUBMIT
        ######
        '''
        subjectBtn = driver.find_element_by_xpath("/html/body/table/tbody/tr[4]/td[3]/form/table/tbody/tr[@id='tr2']/td[@id='td2']/input[1]")
        subjectBtn.click()
        print('Subject btn clicked')
        '''

        ######
        # SAVES NAME/HANDLE OF CURRENT WINDOW (SO WE CAN SWITCH TO NEEDED BROWSER WINDOW)
        ######
        window_before = driver.window_handles[0]

        ######
        # FIND 'DOWNLOAD JSON' BUTTON BY XPATH > CLICK BTN
        ######
        '''
        jsonBtn = driver.find_element_by_xpath("/html/body/table/tbody/tr[3]/td[3]/table/tbody/tr/td/form[2]/input[2]")
        jsonBtn.click()
        '''

        ######
        # WAIT 2SEC > SAVES NAME OF NEW BROWSER WINDOW HANDLE > SWITCHES FOCUS TO NEW WINDOW
        ######
        '''
        time.sleep(2)
        window_after = driver.window_handles[1]
        driver.switch_to_window(window_after)
        '''

        ######
        # FIND TEXT BODY BY XPATH > SAVE TO DICT 'IUPRED' > PRINT TO CONSOLE
        ######
        '''
        json_data = driver.find_element_by_xpath("/html/body")
        iupred = ast.literal_eval(json_data.text)  # uses ast library to read raw json.. it then sees/saves it as a python dict
        # print(iupred['iupred2']) ##use to console testing
        '''

        ######
        # CLOSES ALL WINDOWS PROGRAM OPENED
        ######
        for i in xrange(0, len(driver.window_handles)):
            handle = driver.window_handles[0]
            driver.switch_to_window(handle)
            driver.close()
            print("closed: " + handle)

        iupred = 0

        # print(iupred)
        print("returned iupred")
        if(q != 'q'):
            q.put(iupred)
            print("q passed")
        return iupred

    def testIupredJson(self):
        # function to open iupred.json sample file and reference it in iupredParse()
        data = open('iupredTestFile.json', 'r').read()
        iupred = ast.literal_eval(data)
        return iupred

    def iupredParse(self, userQuery, repThresh=7, mode='iupred'):
        programState = 0  # if zero, program won't run all the way through
        # online or offline test, query

        if(userQuery == 1):
            print("\n---ONLINE MODE ACTIVATED---\n")
            programState += 1
            if(mode == 'iupred'):
                iupred = self.i_scrape()  # calls i_scrape() and sets 'iupred' equal to dict that it returns
            elif(mode == 'blast'):
                iupred = self.blast_scrape()
        elif(userQuery == 2):
            print("\n---OFFLINE MODE ACTIVATED---\n")
            programState += 1
            iupred = self.testIupredJson()  # calls testIupredJson() and sets 'iupred' equal to dict that it returns

        if(programState != 0):
            indexCounter = 0  # counts index position through dictionary
            iCounter = 0  # counts number of x < 0.5 in a row
            jCounter = 0  # counts number of x > 0.5 in a row
            startCounter = 0  # saves beginning of pattern below 0.5
            endCounter = 0  # saves position of end of pattern below 0.5

            startCounterA = 0  # saves beginning of pattern above 0.5
            endCounterA = 0  # saves beginning of pattern above 0.5

            # dicSet = collections.OrderedDict()  # creates empy orderedDict # holds begin/end index in a dictionary
            dicSet = {}
            dicAbc = 97  # chr(dicAbc) will convert to char aka character aka a letter (abc), starting as # will allow to increment
            print('BEGIN TESTING ')
            for i in (iupred['iupred2']):
                if(i < 0.5):

                    # checks to see if this is the first occurance below 0.5
                    if(iCounter == 0):
                        #print("     ---I-COUNTER STARTED")
                        startCounter = indexCounter
                    iCounter += 1

                    # Resets jCounter to zero, we want 5 in a row!
                    if(jCounter >= 1):
                        jCounter = 0
                        #print("     ---J-COUNTER RESET")
                    #print("TESTING:" + str(i) + " iCounter:" + str(iCounter) + " jCounter:" + str(jCounter) + " indexCounter:" + str(indexCounter) + " BELOW!!!!")
                elif(i >= 0.5 and iCounter >= 1):

                    # checks to see if this is the first occurance above 0.5
                    if(jCounter == 0):
                        #print("     ---J-COUNTER STARTED")
                        pass
                    jCounter += 1
                    #print("TESTING:" + str(i) + " iCounter:" + str(iCounter) + " jCounter:" + str(jCounter) + " indexCounter:" + str(indexCounter) + " ABOVE!!!!")

                    # repThresh = 7  # sets number of required repetitions for algorythm

                    if(jCounter >= repThresh):
                        '''
                        if more than 5 in a row are above 0.5 in a set
                        stop counting up iCounter
                        set endCounter = position before start of 5 in a row above 0.5
                        '''

                        if(iCounter >= repThresh):
                            endCounter = indexCounter - jCounter

                            if(len(dicSet) == 0):
                                startCounterA = 0
                                endCounterA = startCounter - 1
                            else:
                                startCounterA = endHolder + 1
                                endCounterA = startCounter - 1

                            endHolder = endCounter

                            ##print(iupred['iupred2'][startCounter:(endCounter + 1)])
                            dicSet.update({chr(dicAbc): [startCounter, endCounter, iupred['seqence'][startCounter:(endCounter + 1)], iupred['iupred2'][startCounter:(endCounter + 1)], startCounterA, endCounterA, iupred['seqence'][startCounterA:(endCounterA + 1)], iupred['iupred2'][startCounterA:(endCounterA + 1)]]})  # updates dictionary with new begin/end vlas
                            '''
                            Each dict entry format:
                            {'a' : [startBelow,endBelow,'amino-acid-sequence[below 0.5]',startAbove,endAbove,'amino-acid-sequence[above 0.5]']}
                            {'b' : [startBelow,endBelow,'amino-acid-sequence[below 0.5]',startAbove,endAbove,'amino-acid-sequence[above 0.5]']}
                            '''

                            dicAbc += 1
                            #print('\n\nTHIS IS THE DICT:\n')
                            #print(str(dicSet) + ": INCLUSIVE")
                            #print('StartNum: ' + str(iupred['iupred2'][startCounter]) + ' EndNum: ' + str(iupred['iupred2'][endCounter]))

                        iCounter = 0
                        jCounter = 0
                        #print("     ---I-COUNTER RESET")
                        #print("     ---J-COUNTER RESET")
                indexCounter += 1
            # In case tail end is all below 0.5, this checks for that and sets end val
            if(iCounter >= repThresh):
                endCounter = (indexCounter - 1)
                ##print("start: (" + str(startCounter) + "): end: " + str(endCounter))
                #print(iupred['iupred2'][startCounter:endCounter + 1])
                dicSet.update({chr(dicAbc): [startCounter, endCounter, str(iupred['seqence'][startCounter:(endCounter + 1)]), str(iupred['iupred2'][startCounter:(endCounter + 1)]), startCounterA, endCounterA, str(iupred['seqence'][startCounterA:(endCounterA + 1)]), str(iupred['iupred2'][startCounterA:(endCounterA + 1)])]})  # updates dictionary with new begin/end vlas
                '''
            Each dict entry format:
            {'a' : [start,end,'amino-acid-sequence[below 0.5]','amino-acid-sequence[above 0.5]']}
            {'b' : [start,end,'amino-acid-sequence[below 0.5]','amino-acid-sequence[above 0.5]']}
            '''
                dicAbc += 1
            print('\n\nTHIS IS THE DICT:\n')
            # print(dicSet)
            if(len(dicSet) <= 0):
                print("\n----Original Sequence Primarily Disordered - None Found")
                return dicSet.update({chr(dicAbc): [1, 1, 1, 1, 1, 1, 1, 1]})
            else:
                return dicSet

    def dicSet2(self):
        # FUNTION TO CREATE AMINO-ACID-SEQUENCE FROM dicSet, ROUGHLY EQUIV TO GETTING INVERSE OF dicSet
        for i in xrange(10):
            i = 0

    # Used in original file before it was a class to run main program
    # Now need to call
    # iupredParse()

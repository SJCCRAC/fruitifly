'''
This program parses a txt file containing proteins to analyse with IUPRED/BLAST/JALVIEW
'''
import warnings  # allows program to be able to ignore benign warnings
#####
# IGNORE WARNINGS
#####
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import requests
import csv  # allows for csv r/w
import pandas as pd  # allows for csv r/w
import json
import mechanize
import webbrowser
import time
import collections  # allows orderedDict
from selenium import webdriver  # for web browser interation/headless browser
from bs4 import BeautifulSoup
import glob
import os.path
import datetime
import urllib2

# Also Using PhantomJS installed via npm (included with NodeJS)


#########################################
############WORKING FUNCTIONS############
#########################################

def parseDataSet(fileName='FruitiData.txt'):
    '''
    Parses original dataSet containing the amino acid sequences of the maternal transcription factors we're interested in.

    Takes in file name as string
    Outputs:
        1. orderdDict
        2. list of dict keys
        3. list of dict vals
    Can call function and set global variable equal to one or all of
    the dataTypes/Sets that this outputs.

    Example:
        variable = parseDataSet()[1]

    This would result in variable being equal to list of all keys in dict
    created.
    '''

    # open dataset text file > create var == to each line as list
    fList = open(fileName).readlines()

    # convert list to dictionary
    fDict = collections.OrderedDict()  # creates empty orderedDict
    ##fDict = {}
    dictVal = ''  # empty string to hold dictVals
    dictKey = ''  # empty string to hold dictKeys
    length = len(fList)

    for line in xrange(0, length):
        #print('inside for')
        #print('line: ' + str(line))
        if(line % 2 == 0):  # if zero or even > use as key
            #print('inside if1')
            dictKey = str(fList[line]).replace('\n', '')
        if(line % 2 != 0):  # if odd > use as value
            #print('inside if2')
            dictVal = str(fList[line]).replace('\n', '')
        if(dictKey != '' and dictVal != ''):
            #print('inside if3')
            fDict.update({dictKey: dictVal})
            dictKey = dictVal = ''

    listFDictKeys = fDict.keys()  # saves dict keys as list
    listFDictVals = fDict.values()  # saves dict vals as list

    # testing prints
    # print(fDict)
    # print(listFDictVals)

    return fDict, listFDictKeys, listFDictVals


# creates timestamp
def timeStamp():
    '''
    returns list = ['mmddyy','hh:mm:ss','Weekday']
    '''
    # ts = time.gmtime()
    ts = time.localtime()
    ts2 = time.strftime('%m%d%y,%H:%M:%S,%d%m%y-%H%M%S,%A', ts)
    ts2 = ts2.split(',')
    return ts2


###############################################
############TESTING BELOW THIS LINE############
###############################################


# creates a csv to write to, add headers row
def csvCreate(listX, listY, csvName='preIupred.csv'):
    '''
    Takes in listFDictKeys, listFDictVals
    '''

    f = csv.writer(open(csvName, 'w'), delimiter=',', lineterminator='\n')
    # f.writerow(['iupred2', 'meta', 'seqence', 'anchor2'])
    f.writerow(['mmddyy', 'hh:mm:ss', 'Key', 'Value', 'example1', 'example2', 'example3'])
    for i in xrange(len(listX)):
        f.writerow((timeStamp()[0], timeStamp()[1], listX[i], listY[i]))


# using 'csv' library open csv > updates specific cell
def csvUpdate():
    '''
    1. Opens preIupred.csv (r)
    2. Opens preIupred.csv (w)
    3. Writes over header names
    4.
    '''

    # read csv file into 'fooReader'
    fooReader = csv.reader(open('preIupred.csv', 'rb'), delimiter=',', lineterminator='\n')
    f = csv.writer(open('preIupred.csv', 'w'), delimiter=',', lineterminator='\n')
    f.writerow(['mmddyy', 'hh:mm:ss', 'Key', 'Value', 'example1', 'example2', 'example3'])
    input = '>Mnt 64.001883 0.822785'

    # read each row in 'fooReader'
    for row in fooReader:
        # define first row column as 'value' for testing
        key = row[2]
        # test if value (1st column) is the same as input (user input)
        if key == input:
            #... if it is then print the 5th column in a certain way
            f.writerow(('FUCKOFF-ITWORKED', '', '', '', '', '', 'hello'))
            #print('this is where the beat drops!')

    '''
    # f.writerow(['iupred2', 'meta', 'seqence', 'anchor2']) #OLD HEADER NAMES, MIGHT USE THEM AGAIN, JUST HERE TO SAVE EM
    # f.writerow(['mmddyy', 'hh:mm:ss', 'Key', 'Value', 'example1', 'example2', 'example3'])
    for i in xrange(5):
        f.writerow(('FUCKOFF-ITWORKED', '', '', '', '', '', 'hello'))
    '''


# using pandas - update csv file at cell level
def csvUpdate2():
    '''
    Pandas Cheatsheet:

    import pandas as pd

    #Open csv and set to var:
    df = pd.read_csv('preIupred.csv')

    #Select single cell by row/column:
    df.iloc([0], [0])
        OR
    df.iat([0], [0])

    #Select single cell by row and column label
    df.loc([0], ['COLUMN-HEADER-NAME'])
        OR
    df.at([0], ['COLUMN-HEADER-NAME'])

    #Select single cell by row and column label
    df.ix[0, 'COLUMN-HEADER-NAME']


    '''
    pd.options.display.max_colwidth = 1000  # sets max string length to display
    df = pd.read_csv('preIupred.csv')  # load csv to var 'df'
    df['example1']  # focuses on column with header 'example1'
    match = df['example1'].str.contains('>Mnt 64.001883 0.822785')
    #print('match: ' + str(match))
    shell = df['Value'][match]
    # print(df)
    # print(df['Key'][match].value_counts())

    # df.set_value(5, 'example1', 'USEFUL-DATA') #updates value of cell at row 5 + header 'Value' to 'CHANGED'
    #df.to_csv('preIupred.csv', index=False)


# creates list holding URLs to visit
def urlCreate():
    pages = []  # empty list to hold URLs to visit

    # create list of urls to visit
    for i in xrange(1, 2):
        url = 'https://iupred2a.elte.hu/'
        # is missing other types of scenarios
        pages.append(url)

    '''
    # opens each URL > sets var to html > sets var to cleaned up html
    for item in pages:
        page = requests.get(item)
        soup = BeautifulSoup(page.text, 'html.parser')
        # print(soup)
    '''

# Demo function


def demo(txtName='FruitiData.txt', csvName='preIupred.csv', dateApndOpt=1):
    if(csvName[-4:] == '.csv'):
        if(dateApndOpt == 1):
            csvNameTime = csvName[:-4] + '_' + timeStamp()[2] + '.csv'
        else:
            csvNameTime = csvName[:-4] + '.csv'
    else:
        if(dateApndOpt == 1):
            csvNameTime = csvName + '_' + timeStamp()[2] + '.csv'
        else:
            csvNameTime = csvName + '.csv'

    listD, listX, listY = parseDataSet(txtName)  # this parses data from file txtName, can insert different file name within same directory
    '''
    1. Calls function to parse data set from FruitiData.txt then saves/outputs as ordered dict
    2. Calls function that takes parsed data from step one and then saves it to a csv 'collectData1.csv'
    '''
    csvCreate(listX, listY, csvNameTime)  # this takes in vars from 'parseDataSet()' > creates/writes to csv

# csvUpdate()

# csvUpdate2()


# csvUpdate()  # uncomment to continue testing this


# csvUpdate2()  # updates csv at cell level using pandas (seems best method)

# demo() # uncomment to run main program


def blastParse(fileName='PFK3E0EY015-Alignment.json', jalName='jalViewFile.fa'):
    with open(fileName) as json_file:
        data = json.load(json_file)
        # print(type(data))
        # print(json.dumps(data, indent=2)) #pretty printed
        # for i in xrange(10):
        #     print(data['BlastOutput2'][0]['report']['results']['search']['hits'][2]['hsps'][i])
        #     print('')
        #     print('')

        dictHolder = {}
        iterMain = data['BlastOutput2'][0]['report']['results']['search']['hits']
        f = open(jalName, 'w')
        f.write('')
        fl = open(jalName, 'a')
        for i in xrange(4):
            print '#########################'
        for item in xrange(len(iterMain)):
            subject = data['BlastOutput2'][0]['report']['results']['search']['hits'][item]['hsps']
            title = data['BlastOutput2'][0]['report']['results']['search']['hits'][item]['description'][0]['title']
            sciName = str(data['BlastOutput2'][0]['report']['results']['search']['hits'][item]['description'][0]['sciname'])
            dictHolder[sciName] = dictHolder.get(sciName, 0) + 1
            if(dictHolder[sciName] == 1):
                fl.write('\n' + '> ' + sciName)
                print("title: " + str(title))
                print("sciname: " + str(sciName))
                subHolder = ''
                for i in xrange(len(subject)):
                    subHolder += str(subject[i]['hseq'])
                    print("index: " + str(i) + " subject: " + str(subject[i]['hseq']))
                print("subjectFull: " + str(subHolder))
                fl.write('\n' + str(subHolder))
                print('\n\n')
        print(dictHolder)
        fl.close()

        # print data['BlastOutput2'][0]['report']['results']['search']['hits'][0]['description'][0]['title']

    # fList = open(fileName).readlines()
    # print fList

    '''
      # open dataset text file > create var == to each line as list
    fList = open(fileName).readlines()

    # convert list to dictionary
    fDict = collections.OrderedDict()  # creates empty orderedDict
    ##fDict = {}
    dictVal = ''  # empty string to hold dictVals
    dictKey = ''  # empty string to hold dictKeys
    length = len(fList)

    for line in xrange(0, length):
        #print('inside for')
        #print('line: ' + str(line))
        if(line % 2 == 0):  # if zero or even > use as key
            #print('inside if1')
            dictKey = str(fList[line]).replace('\n', '')
        if(line % 2 != 0):  # if odd > use as value
            #print('inside if2')
            dictVal = str(fList[line]).replace('\n', '')
        if(dictKey != '' and dictVal != ''):
            #print('inside if3')
            fDict.update({dictKey: dictVal})
            dictKey = dictVal = ''

    listFDictKeys = fDict.keys()  # saves dict keys as list
    listFDictVals = fDict.values()  # saves dict vals as list

    # testing prints
    # print(fDict)
    # print(listFDictVals)

    return fDict, listFDictKeys, listFDictVals
    '''


def openDownloads():
    list_of_files = glob.glob("C:/Users/SJCCRAC/Documents/Python Code")  # * means all if need specific format then *.csv
    latest_file = max(list_of_files, key=os.path.getctime)
    print list_of_files
    print latest_file


# blastParse() #runs blastParse function


def downloadUrl():
    print('Beginning file download with urllib2...')
    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID=P09YHPX0014&FORMAT_TYPE=JSON2_S&FORMAT_OBJECT=Alignment&CMD=Get'

    filedata = urllib2.urlopen(url)
    datatowrite = filedata.read()

    with open('/Users/SJCCRAC/Documents/Python Code/testDownload.json', 'wb') as f:
        f.write(datatowrite)
    print(datatowrite)

# openDownloads() # tests openDownloads() functions


# downloadUrl()


demo('7_proteins.txt', 'preIupred.csv', 1)  # (txtName='FruitiData.txt', csvName='preIupred.csv', apndDate[1=yes, 0=no])
'''
Parses original formatted amino acid sequence data
Outputs is to csv file that you specify, default = 'preIupred.csv'
'''

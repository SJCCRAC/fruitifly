import headlessClass
import FruitiFly as ff
import collections
import csv
import time
import pandas as pd


# emptyDict = {}

# creates timestamp


def timeStamp():
    '''
    returns list = ['mmddyy','hh:mm:ss','Weekday']
    '''
    # ts = time.gmtime()
    ts = time.localtime()
    ts2 = time.strftime('%m%d%y,%H:%M:%S,%A', ts)
    ts2 = ts2.split(',')
    return ts2


def iupredParseToCsv(txtName='7_proteins.txt', repThresh=7):
    print("repThresh: " + str(repThresh))
    userQuery = 1  # int(raw_input('Choose Mode: 1-Online or 2-Offline: '))

    fDict, dictKeys, dictVals = ff.parseDataSet(txtName)

    emptyDict = collections.OrderedDict()

    for i in xrange(len(dictKeys)):
        # print("\n" + str(i + 1) + " " + str(dictVals[i]))
        maxSeq = ''
        maxSeqNum = 0

        emptyDict.update({dictKeys[i]: headlessClass.headless(dictVals[i]).iupredParse(userQuery, repThresh)})
        filename = str(i) + '_' + str(dictKeys[i][1:]) + '_P1.csv'
        print("filename: " + filename)
        f = csv.writer(open(filename, 'w'), delimiter=',', lineterminator='\n')
        f.writerow(['mmddyy', 'hh:mm:ss', 'Seq', 'IUPRED2'])
        dictKey = dictKeys[i]

        print("this is the dict")
        time.sleep(2)
        print emptyDict
        time.sleep(10)

        if(emptyDict[dictKey] == None):
            pass
            print("passed")
        else:
            print(emptyDict[dictKey]['a'])
            indexCounter = 0
            dicAbc = 97  # chr(dicAbc) will convert to char
            for Abc in emptyDict[dictKey]:
                subSet = chr(dicAbc)
                seqLen = len(emptyDict[dictKey][subSet][6])
                if(seqLen > maxSeqNum):
                    maxSeq = emptyDict[dictKey][subSet][6]  # the disordered sequence

                beginA = emptyDict[dictKey][subSet][4]
                endA = emptyDict[dictKey][subSet][5]

                beginEndA = str(str(beginA) + ':' + str(endA))

                f.writerow(('TITLE: ', dictKey, subSet, seqLen))
                f.writerow(('Begin:End', beginA, endA, ''))
                # for group in emptyDict[dictKey][Abc]:

                for item in emptyDict[dictKey][subSet][6][0:]:
                    dictABCVal = emptyDict[dictKey][subSet][6][indexCounter]
                    dictIupred = emptyDict[dictKey][subSet][7][indexCounter]
                    print("Index: " + str(indexCounter) + ' - ' + str(dictABCVal) + ' - ' + str(dictIupred))

                    f.writerow((timeStamp()[0], timeStamp()[1], dictABCVal, dictIupred))
                    indexCounter += 1
                    # print item
                for i in xrange(3):
                    f.writerow(('', '', '', ''))
                indexCounter = 0
                dicAbc += 1
            # for i in emptyDict[dictKey]:
            #     f.writerow((timeStamp()[0], timeStamp()[1], SEQUENCE, IUPRED2))
        f.writerow(('longest seq: ', maxSeq))

# for i in xrange(len(emptyDict[1]['a'][6])):
#     print(str(emptyDict[1]['a'][6][i]) + ' ' + str(emptyDict[1]['a'][7][i]))


def blastSubsetToCsv():
    #subSet = ['SSESEQKDEKEELLCPKPQIDCTNTDLEQSTAIETDTEQVEEKRSNRRKSRRIRNEKFKTETDTLSDHLDAKKAENASLEISMRPKCTLETQQSDPVTAKNKRNSGRLSRKEKSVINAAKSEKDKSPSAISQSTERKQLLNENPSKKDKKTEQSGNKKEAVVGPLDKTETSSSTNIIDKKSNESFDSAMQPSDRLNQKESAFTKLSSISSPKKIMKDQDKDLDALSKGGDSNPTIRDTGEDSRQTDKKHQENDTKHEEEDSSKLKANIDETKSSSEKDAEPISKDSSQDSAKPRLSKPKSRNKRKKNEKKPNDSIAESDIEGGFQVNTETVQATCSTPSESNKKDMVKSDETNEEPNLSETEIGRIRKRGQAFHIENPKDDLHITPQNENQ', 'ARLSDSGTSASGSSSSSSSSSDSAMGGEVVPMPGPGETLQLPGVPAAITTVMRVQPTQSQKAPPSNSVTLTPILPLPTSPKQRQKKPRKKKAITSAPILDSSDDDEPPPKHPGLDHTAVSVQTQPATDTVKKGRGRPRKQQQSGGSGNLSSASAGSSSQTKGPTLTAAKKPLAKTPLAMSRARKREHSSQSSSNGNTPTKKVATPQLVAAPLKPTSNTAGSSSSDEDSSSSAESSSKSSSSSSSSDDTETQNTNCRIVKLNKTGAVQKKALLGSGSSSPSSSGSEAEDQTTRSQVGSGQALAQQLPPYKQLPISQHS', 'ARLSDSGTSASGSSSSSSSSSDSAMGGEVVPMPGPGETLQLPGVPAAITTVMRVQPTQSQKAPPSNSVTLTPILPLPTSPKQRQKKPRKKKAITSAPILDSSDDDEPPPKHPGLDHTAVSVQTQPATDTVKKGRGRPRKQQQSGGSGNLSSASAGSSSQTKGPTLTAAKKPLAKTPLAMSRARKREHSSQSSSNGNTPTKKVATPQLVAAPLKPTSNTAGSSSSDEDSSSSAESSSKSSSSSSSSDDTETQNTNCRIVKLNKTGAVQKKALLGSGSSSPSSSGSEAEDQTTRSQVGSGQALAQQLPPYKQLPISQHS', 'ARLSDSGTSASGSSSSSSSSSDSAMGGEVVPMPGPGETLQLPGVPAAITTVMRVQPTQSQKAPPSNSVTLTPILPLPTSPKQRQKKPRKKKAITSAPILDSSDDDEPPPKHPGLDHTAVSVQTQPATDTVKKGRGRPRKQQQSGGSGNLSSASAGSSSQTKGPTLTAAKKPLAKTPLAMSRARKREHSSQSSSNGNTPTKKVATPQLVAAPLKPTSNTAGSSSSDEDSSSSAESSSKSSSSSSSSDDTETQNTNCRIVKLNKTGAVQKKALLGSGSSSPSSSGSEAEDQTTRSQVGSGQALAQQLPPYKQLPISQHS', 'MALPSNYKQIAVGGQGSATPLQGGGGGSGGGSRSRSSGGGGGGDRNKDQTPIFTHSNYGNPAFTPQKVTKSSSSKNQNESRLPKPPKPPEKPI', 'NSGPLSGHSIDQQQQQVHQADGLGMGGGGGGGVGADGMHCPVTTGLPPISSFRPTSGGIGGPGAGQQAPVNVNVNPPAVFNSPQAHNHNHTVQAQHSALSTAGPLGHHSLNHTPHAHSHTLPLPHALPHGHTLPHPHHSQQNSPAVQSSDA']
    subSet = ['RQRLDEPNANANSAANPAAAAATAAAAPVTRSEEAVKPPTQPGQPAATPAGQEPASAVPAPAAPPKETPPAVKPATLNPTPSSTPTPAPAVHVHETASKTDPEPMDIEPPPKPSVPPPPIKPEKLEMAAALPPQSTLVEPPKTEPAKVVAQPGKVPTPVPTPTPPPEVAPAAGAATAATTTTTTTPTPPPVPQQVPLPQQQGGPAPPPGMPQMHPHPGHPPPGHSLMPPHMGPHQPPPGMPGLPPPPPHTGYANYGGPPHGPPPGPPGGPARPYYQPQYGGHPTPQPYYAPFSPYQQSYGPPPGSHYMSPRPPPPQHNGNPGHPYAPEHGSNPPPPQQQQQQQPPPGHLHEPSGGGPGAPGGGAGA']
    for i in xrange(len(subSet)):
        headlessClass.headless(subSet[i]).blast_scrape()


###
# Main Program Functions Can Be Called By Uncommenting Below:
###

# # Calls main function in this program
# iupredParseToCsv()
# '''
# First step in process. Save proteins that you'd like run through IUPRED to this file: 7_proteins.txt
# then run iupredParseToCsv(). It'll save a csv named protein name + timestamp. Find Longest protein in this csv and copy string.
# '''


# # Processes a string through blast > downloads results
# blastSubsetToCsv()
# '''
# Second step in process. Add proteins, that you want processed through BLAST, as strings to
# subSet variable in blastSubsetToCsv(). Then run blast SubsetToCsv. It'll click download on the BLAST results
# page to download the results as a JSON file.
# '''


# # # Opens sample Blast json output file > parses is for desired data
# # # outputs it to jalViewFile.fa be default
# saveAsList = ['CG1620 84.276617 0.846827.fa']
# fileNamesList = ['CG1620 84.276617 0.846827.json']


# for i in xrange(len(saveAsList)):
#     # ff.blastParse(fileNamesList[i], saveAsList[i])
#     ff.blastParse("PFV3CYBG01R-Alignment.json", "PFV3CYBG01R-Alignment.fa")
'''
LAST STEP IN PROCESS.
'''

ff.blastParse("PFY6TT52015-Alignment.json", "Bap111 71.669550 0.857143.fa")

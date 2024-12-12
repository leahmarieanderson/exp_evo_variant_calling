import sys
import os

def getAverageCoverage(summary):
	## Interpert the summary file and get the total read coverage
	summaryFile = open(summary, 'r')
	for line in summaryFile:
		currentLine = line.strip().split()
		if (currentLine[0] == "Total"):
			averageCoverage = float(currentLine[2])
			# print("Average Coverage =", averageCoverage)
	summaryFile.close()
	
	return averageCoverage	

def makeAllDict(allWig, chrDict, chrDict2):
	## Return a dictionary of the allWig read number with the key of 'Chr_Pos'
	allWigFile = open(allWig, 'r')
	wigFileLines = allWigFile.readlines()
	allWigFile.close()
	
	i = 0
	chr = 1

	EdgeWindow = 0
	
	AllFiltDict = {}
	
	for i in range(len(wigFileLines)):
		if wigFileLines[i].split()[0] == 'track': #Identify the first line and skips it
			i += 2
		elif (wigFileLines[i].split()[0] == 'variableStep') & ((i + 2) < len(wigFileLines)):
			if (wigFileLines[i+2].split()[0] != 'variableStep'):
				chr = chrDict2[wigFileLines[i].split()[1][9:]][0]
		elif wigFileLines[i].split()[0] == 'variableStep':
			i = len(wigFileLines) + 1
		elif chr == 'M':
			continue
		else:
			splitLine = wigFileLines[i].strip().split()
						
			## Skips the ends of chromosomes
			if (int(splitLine[0]) < EdgeWindow) or (int(splitLine[0]) > ((chrDict[str(int(chr) + 1)] - (chrDict[str(int(chr))]) - EdgeWindow))):
				#print "Skipped"
				continue
			AllFiltDict[str(chr) + "_" + splitLine[0].strip()] = splitLine[1]
			#outFile.write(str(chr) + '\t' + str(chrDict[str(chr)] + int(splitLine[0].strip())) + '\t' + str(splitLine[1]) + '\t' +  str(averageCoverage) + '\n')
	return(AllFiltDict)
	
def normalizeReadDepth(averageCoverage, allDict, filtWig, ploidy, out):
## Take in the two wig files and a summary file and output a tab delimited file that has:
## Chrom \t Position \t NormalizedCoverage \t ProportionPassingFilter
	
	filtWigFile = open(filtWig,'r')
	outFile = open(out,'w')
	
	## What to normalize by and filter the ends of the chromosomes by 
	BasePloidy = float(ploidy)
	EdgeWindow = 5000
	CutOff = 0.90
	
	wigFileLines = filtWigFile.readlines()

	i = 0
	chr = 1
	
	for i in range(len(wigFileLines)):
		if wigFileLines[i].split()[0] == 'track': #Identify the first line and skips it
			i += 2
		elif (wigFileLines[i].split()[0] == 'variableStep') & ((i + 2) < len(wigFileLines)):
			if (wigFileLines[i+2].split()[0] != 'variableStep'):
				chr = chrDict2[wigFileLines[i].split()[1][9:]][0]
		elif wigFileLines[i].split()[0] == 'variableStep':
			i = len(wigFileLines) + 1
		elif chr == 'M':
			continue
		else:
			splitLine = wigFileLines[i].strip().split()
			
			## Skips the ends of chromosomes
			if (int(splitLine[0]) < EdgeWindow) or (int(splitLine[0]) > ((chrDict[str(int(chr) + 1)] - (chrDict[str(int(chr))]) - EdgeWindow))):
				#print "Skipped"
				continue
			allDictNumber = allDict[str(chr) + "_" + splitLine[0].strip()]

			#print str(float(splitLine[1]) / float(allDictNumber))
			
			if((float(splitLine[1]) / float(allDictNumber)) > CutOff):
				outFile.write(str(chr) + '\t' + str(chrDict[str(chr)] + int(splitLine[0].strip())) + '\t' + str((float(splitLine[1]) / averageCoverage) * BasePloidy) + '\n')
	filtWigFile.close()
	outFile.close()

### MAIN ###
#summaryFile = open(sys.argv[1], 'r')
#filtWigFile = open(sys.argv[2],'r')
#allWigFile = open (sys.argv[3], 'r')
ploidy = sys.argv[4]
#outFile = open(sys.argv[5],'w')

## Dictionary with global chromosome lengths
chrDict = {'1':0, 
    '2':230218,
    '3':1043402,
    '4':1360022,
    '5':2891955,
    '6':3468829,
    '7':3738990,
    '8':4829930,
    '9':5392573,
    '10':5832461,
    '11':6578212,
    '12':7245028,
    '13':8323205,
    '14':9247636,
    '15':10031969,
    '16':11123260,
    '17':12157105,
    'M':''}
        
## Dictionary with global chromosome start positions       
chrDict2 = {'I':['1',0],
	'II':['2',230218],
	'III':['3',1043402],
	'IV':['4',1360022],
	'V':['5',2891955],
	'VI':['6',3468829],
	'VII':['7',3738990],
	'VIII':['8',4829930],
	'IX':['9',5392573],
	'X':['10',5832461],
	'XI':['11',6578212],
	'XII':['12',7245028],
	'XIII':['13',8323205],
	'XIV':['14',9247636],
	'XV':['15',10031969],
	'XVI':['16',11123260],
	'M':['M',0]}

makeAllDict(sys.argv[3], chrDict, chrDict2)

normalizeReadDepth(getAverageCoverage(sys.argv[1]), makeAllDict(sys.argv[3], chrDict, chrDict2), sys.argv[2], ploidy, sys.argv[5])

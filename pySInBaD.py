#!/usr/bin/env python
"""
SInBaD Galaxy wrapper

CLI Usage:
python pySInBaD.py testdata/TestFile-hs-hg19.csv csv -1 > out_file 

Requirements: 
"""

import os
import sys
import operator
import scipy as sp

FDRdefDIC =    {"c":0.5,"i":0.5,"p":0.5}
FDR1DIC   =    {"c":0.85693359375,"i":0.8359375,"p":0.8790311865}
FDR5DIC   =    {"c":0.732421875,"i":0.458984375,"p":0.784355163574}
FDRallDIC =    {"c":0.0,"i":0.0,"p":0.0}
dir_database = "testdata/TestDB-SInBaD.dat"

def search(FH, data):
    """
    Core search module
    """
    l=0
    h=os.fstat(FH.fileno())[6]
    while l<h:
        mid=(l+int(h))/2
        FH.seek(mid)
        dumpfile=FH.readline()
        midvalue=FH.readline()
        midvalue=midvalue.split("\t")
        try:
	 position=int(midvalue[0])
         if position<int(data[1]):
            l=mid+1
         elif position>int(data[1]):
            h=mid-1
         else:	
            return midvalue
	except ValueError:
	 data=[0,0]
	 return data
    if mid==0:
	FH.seek(0)
        midvalue=FH.readline()
	midvalue=midvalue.split("\t")
	if int(data[1])==int(midvalue[0]):
	    return midvalue
    return data

def sort_data(data,col=0):
    return sorted(data,key=operator.itemgetter(col))

def runwithsearch(searchitems,fn_database,fdrthres):
    """
    The result part
    """
    searchitems_sorted=sort_data(searchitems,1)
    posin={'A':0,'C':1,'G':2,'T':3}   
    FH_database=open(fn_database,'r')
    for item in searchitems_sorted:
        record=search(FH_database,item)
	if len(record)>5:
	    if len(item)>2:
	        refal=item[2].strip(' ')#.strip(' \r')
		altal=item[3].strip(' ')#.strip(' \r')
		if float(record[1+4*posin[refal]+posin[altal]])>=fdrthres[record[len(record)-1].strip("\n")]:
		    print item[0]+","+record[0]+","+refal+","+altal+","+record[1+4*posin[refal]+posin[altal]]+","+record[len(record)-1],
	    else:
	        for refal in posin.keys():
		    for altal in posin.keys():
			if float(record[1+4*posin[refal]+posin[altal]])>=fdrthres[record[len(record)-1].strip("\n")]:
		            print item[0]+","+record[0]+","+refal+","+altal+","+record[1+4*posin[refal]+posin[altal]]+","+record[len(record)-1],
    FH_database.close()

def processFormdata(fn_in):
    """
    Processing CSV format data  
    """

    FH = open(fn_in, 'r')
    datasplit = FH.readlines()
    FH.close()
    datasplit = [x.strip('\n') for x in datasplit]
    databychr={}
    for i in xrange(0,len(datasplit)):
        if datasplit[i][0:1]!="#":
            datasplit[i]=datasplit[i].split(",")
            if len(datasplit[i])>1:
                datasplit[i][1]=int(datasplit[i][1])
		try:
                    databychr[datasplit[i][0]].append(datasplit[i])
    		except KeyError:
		    databychr[datasplit[i][0]]=[datasplit[i]]
    return databychr    

def processFormdatavcf(fn_in):
    """
    Processing VCF format 
    """

    databychr = {}
    FH        = open(filename,'r')
    formatcol = -1
    genotypes = []
    for line in FH:
	line=line.strip("\n")
        if line[0:2]=="##":
	    line=line.strip("#")
	    linesplit=line.split("=")
	    if linesplit[0]=="fileformat":
		version=linesplit[1].strip("VCFv")
		if not float(version)>=4.0:
		    print "Wrong version"
		    return data
	elif line[0:1]=="#":
	    line=line.strip("#")
	    linesplit=line.split("=")
	    for i in xrange(0,len(linesplit)):
		if formatcol!=-1:
		    genotypes.append(linesplit[i])
		if linesplit[i]=="FORMAT":
		    formatcol=i
	else:
	    linesplit=line.split("\t")
	    if linesplit[0][0:3]!="chr":
		linesplit[0]="chr"+linesplit[0]
	    temp=[int(linesplit[1]),linesplit[3],linesplit[4]]
	    try:
		data[linesplit[0]].append(temp)
	    except KeyError:
		data[linesplit[0]]=temp
    FH.close()
    return databychr    

def main():    
    """
    main function 
    """

    try:
        fn_in  = sys.argv[1]
        FILE_TYPE = sys.argv[2]
        FDR = int(sys.argv[3])
    except:
        print __doc__
        sys.exit(-1)

    ## FDR setting 
    if FDR == -1:
        fdrthres = FDRallDIC
    elif FDR == 1:
        fdrthres = FDR1DIC
    elif FDR == 5:
        fdrthres = FDR5DIC
    elif FDR == 0:
        fdrthres = FDRdefDIC
    else:
        print 'FDR not in the default range'
        print 'Cannot continue, exiting...'
        sys.exit(-1)
    
    ## Input file formats 
    if FILE_TYPE == "vcf":
        filedata = processFormdatavcf(fn_in)
    elif FILE_TYPE ==  "csv":
        filedata = processFormdata(fn_in)
    else:
        print 'More Input file formats in TODO list!'
        print 'Currently supports CSV, VCF'
        sys.exit(-1)

    for chr in filedata.keys():    
        name = chr.strip("chr")
        runwithsearch(filedata[chr],dir_database,fdrthres)

if __name__=="__main__":
    main()

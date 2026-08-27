#!/usr/bin/env python

import sys
from collections import defaultdict
#out=open(sys.argv[1].rsplit('.', 1)[0]+'.bed','w')

def parseAtt(attributes):
    #print attributes
    attrib=dict(att.split() for att in attributes.rstrip(';').split(';'))
    #target=attrib['gene_id'].strip('\"').split('=')[0]+'_'+attrib['transcript_id'].strip('\"').split('=')[0] #
    #target=attrib['gene_id'].strip('\"')+'|'+attrib['transcript_id'].strip('\"')
    target=attrib['transcript_id'].strip('\"')
    gene=attrib['gene_id'].strip('\"')#.strip('\"').split('=')[0]+'_'+attrib['transcript_id'].strip('\"').split('=')[0] #
    #name=attrib['gene_name'].strip('\"')#.strip('\"').split('=')[0]+'_'+attrib['transcript_id'].strip('\"').split('=')[0] #
    #btype=attrib['gene_type']
    #if '=' in attrib['gene_id']:
    #    target+='_'+attrib['gene_id'].strip('\"').split('=')[1].replace('Hypothetical_Protein','HyPr')

    #gene=elts[0][1].strip('\"')
    return target,gene#,btype

#GeneNames=dict(line.rstrip().split("\t") for line in open(sys.argv[2]))


exonNb=0
exonStarts=[]
exonEnds=[]
cdsStarts=[]
cdsEnds=[]
item=[]
genb=0
#prev_target=''
targetSet=set()
trCount=defaultdict(int)
out=open(sys.argv[1].rsplit('.',1)[0]+'.bed','w')
for line in open(sys.argv[1]):
    if line.startswith('##') or line=='\n': continue
    (scaf,source,feat,start,end,score,strand,phase,attributes)=line.rstrip().split('\t')
    #print attributes
    if feat=='gene': continue
    target,gene=parseAtt(attributes)
    #if btype=='pseudogene': continue
    #target=parseAtt(attributes)
    #print target
    if not target in targetSet:
        trCount[gene]+=1
        #print target,'-',prev_target,exonNb
        targetSet.add(target)
        if exonNb>0:
            exonSizes=[]
            exonCoords=[]
            item[1]=min(exonStarts)-1
            item[2]=max(exonEnds)
            for eStart,eEnd in zip(exonStarts,exonEnds):
                exonSizes.append(eEnd-eStart+1)
                exonCoords.append(eStart-int(item[1])-1)
            #print exonSizes[0:5]
            #print exonCoords[0:5]
            if len(exonCoords)>1 and exonCoords[0]>exonCoords[1]:
                exonSizes=reversed(exonSizes)
                exonCoords=reversed(exonCoords)
            if not len(cdsStarts+cdsEnds)==0:
                item[6]=min(cdsStarts+cdsEnds)-1
                item[7]=max(cdsStarts+cdsEnds)
                item[8]=0
                item[9]=exonNb
            else:
                item[6]=min(exonStarts+exonEnds)-1
                item[7]=max(exonStarts+exonEnds)
                item[8]=0
                item[9]=exonNb
            #item.extend([min(exonStarts+exonEnds),max(exonStarts+exonEnds),0,exonNb])
            item[10]=','.join(map(str,exonSizes))+','
            item[11]=','.join(map(str,exonCoords))+','
            exonNb=0
        if not cdsStarts==[]: 
             out.write('\t'.join(map(str,item))+'\n')
#Re-initialise
        exonStarts,exonEnds=[],[]
        cdsStarts,cdsEnds=[],[]
        #gid='{0}|{1}'.format(target,name)
        #gid='{0}.{1}'.format(name,trCount[gene])
        item=[scaf,start,end,target,1000,strand,'','','','','','',gene]
        #item=[scaf,start,end,gid,1000,strand,'','','','','','',gene]

    if feat=='exon':
        #print 'exon',start,end
        if not int(start) in exonStarts and not int(end) in exonEnds:
            if not start==end:
                exonStarts.append(int(start))
                exonEnds.append(int(end))
                exonNb+=1

    if feat=='CDS':
        #print 'exon',start,end
        if not start==end:
            cdsStarts.append(int(start))
            cdsEnds.append(int(end))
            #exonNb+=1

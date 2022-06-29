#!/usr/bin/env python
import gzip
import os
import sys

#Input
InputCNVvcf=#Output.vcf.gz from CNV_FACETS
LOHRegion=#Output of LOH region
InputVCF=#Filtered VCF from Small_Nucleotide_Variant_Calling_Filtering.sh
TumorBam=#Tumor Bam File


#Output
FilteredVCFbed=#bed format output of filtered VCF 
LOHCandidateVariant=#Intersected SNV with LOH region
ASLOHOutput=#



#Used Tool
Bedtools=#Path of Bedtools
IGVtools=#Path of igvtools



class SNV:

	def __init__(self):
		self.sChromosome=""
		self.nChromosome=""
		self.nPosition=""
		self.sReferenceAllele=""
		self.sAlternativeAllele=""
		self.sGenotype=""
		
	
	def ParseSNV(self,sLine):
		
		sLine=sLine.strip()
		t=sLine.split("\t")
		
		
		(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
		sSample=t[-1]
		sGenotype=sSample.split(":")[0]
		self.sChromosome=sChr
		
		self.nPosition=int(nPosition)
		self.sReferenceAllele=sRef
		self.sAlternativeAllele=sAlt
		self.sGenotype=sGenotype
		
		
		if sChr=="X":
			self.nChromosome=23
		else:
			self.nChromosome=int(sChr)
	
		





def ExtractingLOHRegion(InputCNVvcf):
	
	fp=gzip.open(InputCNVvcf,"rt")
	fout=open(LOHRegion,"w")
	for sLine in fp.readlines():
		
		if sLine[0]=="#":
			pass
		else:
			sLine=sLine.strip()
			t=sLine.split("\t")
			sFilter=t[6]
			if sFilter=="PASS":
				
				sINFO=t[-1]
				lINFO=sINFO.split(";")
				(sCNVTYPE,nTCM_EM,nEndPosition)=(lCNVINFO[0],lCNVINFO[-3],lCNVINFO[2])
				(sCNVTYPE,nTCM_EM,nEndPosition)=(sSVTYPE.split("=")[1],nTCM_EM.split("=")[1],nEndPosition.split("=")[1])
				if ((sCNVTYPE=="DUP-LOH") or (sCNVTYPE=="HEMIZYG") or (sCNVTYPE=="LOH")):
					fout.write("{0}\n".format("\t".join([t[0].replace("chr",""),t[1],nEndPosition,sCNVTYPE,nTCM_EM])))
				
				
				
			
			
			
	fout.close()
		
	

def LowQuality_Variant_Filtering(FilteredVCF):
	
	fp=open(InputVCF)
	fout=open(FilteredVCFbed,"w")
	lSNV=[]
	for sLine in fp.readlines():
		if sLine[0]=="#":
			pass
		else:
			t=sLine.split("\t")
			if t[6]=="PASS":
				sSample=t[-1]
				sGenotype=sSample.split(":")[0]
				if sGenotype=="0/1":
				
					(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
					if (sChr=="chrM"):
						pass
					else:
						sVariant=SNV()
						sVariant.ParseSNV(sLine)
						lSNV.append(sVariant)
						
	lSNV.sort(key = lambda x: x.nPosition)
	lSNV.sort(key = lambda x: x.nChromosome)
	for sVariant in lSNV:
	
		fout.write("{0}\n".format("\t".join([sVariant.sChromosome,\
		sVariant.nPosition,sVariant.nPosition,sVariant.sReferenceAllele,sVariant.sAlternativeAllele,sVariant.sGenotype])))
	
	fout.close()
	
	
	

def Intersection_for_LOH_Candidate_Variant(FilteredVCFbed,LOHRegion):
	
	
	
	os.system(Bedtools+" intersect -a "+LOHRegion+" -b "+FilteredVCFbed+" -wo > "+LOHCandidateVariant)
	
	
	




def Igvtools(TumorBam,LOHCandidateVariant):
	
	
	sID=sFile.split(".")[0]
	try:
		os.mkdir("IGVtools/")
	except:
		pass
	
	fp=open("Intersection/Intersected_"+sFile)
	
	sBamFile="/mnt/towel/TCGA/Finalpostrecalbam/"+sID+"-T.bam"
	
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition,sRef,sAlt,sSS)=(t[5],t[6],t[8],t[9],t[10])
		sKey="_".join([sChr, nPosition,sRef,sAlt,sSS])
		os.system(IGVtools+" count -w 1 --bases --query chr"+sChr+":"+nPosition+"-"+nPosition+" "+sBamFile+ \
		#print("igvtools count -w 1 --bases --query chr"+sChr+":"+nPosition+"-"+nPosition+" "+sBamFile+ \
	" IGVtools/"+ sKey + ".bam.wig hg19")
		#break





def ParseIGVtools(sFile,dIGVdict):
	
	fp=open(sFile)
	fp.readline()
	fp.readline()
	fp.readline()
	sLine=fp.readline()
	sLine=sLine.strip()
	#17_72436932_G_C_Somatic_het
	sFile=sFile.split("/")[-1]
	sKey=sFile.split(".")[0]
	(sChr,nPosition,sRef,sAlt,sSS,sGenotype)=sKey.split("_")
	try:
		(nPos, nA, nC, nG, nT, nN, nDel, nIns)=sLine.split("\t")
		dIGVdict[sKey]={"A":nA,"G":nG,"C":nC,"T":nT}
	except ValueError:
		dIGVdict[sKey]="Deletion"
	
	
	
	
	
	return dIGVdict


def Allele_specific_LOH_Variant_Calling(TumorBam,LOHCandidateVariant):
	
	Igvtools(TumorBam,LOHCandidateVariant)
	
	dIGVdict=dict()
	#IGVtools/"+sID+"/"+ sKey + ".bam.wig
	lIGVToolsFiles=glob.glob("IGVtools/"+"*.bam.wig")
	
	for sIgvcall in lIGVToolsFiles:
		ParseIGVtools(sIgvcall,dIGVdict)
	#print(dIGVdict["1_17265560_C_T_Germline_het"])
	fp=open(LOHCandidateVariant)
	
	fout=open(ASLOHOutput,"w")
	fout.write("Chr\tPosition\tRef\tAlt\tSS\tRefreadinTumor\tAltreadinTumor\tAARF\tCNV\n")
	#for sLin
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition,sRef,sAlt,sSS)=(t[5],t[6],t[8],t[9],t[10])
		sKey="_".join([sChr, nPosition,sRef,sAlt,sSS])
		
		
		try:
			
			if not dIGVdict[sKey]=="Deletion":
				
				
				(sRefinTumor,sAltinTumor)=(dIGVdict[sKey][sRef],dIGVdict[sKey][sAlt])
				try:
					#fout.write("{0}\t".format("\t".join([sChr, nPosition,sRef,sAlt,sSS])))
					if float(sRefinTumor)==0:
						sAllreads=dIGVdict[sKey].values()
						sTotalRead=sum(map(float,sAllreads))
						nAARF=float(sAltinTumor)/float(sTotalRead)
					else:
						nAARF=float(sAltinTumor)/(float(sRefinTumor)+float(sAltinTumor))
					if "0/1" in sSS:
					
						if nAARF<0.2:
							sSSwithCNV="RLOH"
						elif nAARF>0.8:
							sSSwithCNV="LOH"
						else:
							sSSwithCNV=sSS
					else:
						pass
						#sys.exit()
					fout.write("{0}\t".format("\t".join([sChr, nPosition,sRef,sAlt,sSS])))
					fout.write("{0}\n".format("\t".join([sRefinTumor,sAltinTumor,str(nAARF),sSSwithCNV])))
				except ZeroDivisionError:
					pass	
			else:
				pass
				#fout.write("{0}\n".format("\t".join(["0","0","0",sSS])))
		except KeyError:
			pass
	
	
	






if __name__=="__main__":
	
	Extracting_LOH_Region(CNVOutput)
	LowQuality_Variant_Filtering(FilteredVCF)
	Intersection_for_LOH_Candidate_Variant(FilteredVCFbed,LOHRegion)
	Allele_specific_LOH_Variant_Calling(TumorBam,LOHCandidateVariant)
	
	
	
	
	





#!/usr/bin/env python
import gzip
import os
import sys
import glob
#Input
InputCNVvcf="../../Data/Exampledata_for_ASLOH/Example_CNV.vcf.gz" #Output.vcf.gz from CNV_FACETS
InputVCF="../../Data/Exampledata_for_ASLOH/Example_SNV.vcf.gz" #Filtered VCF from Small_Nucleotide_Variant_Calling_Filtering.sh
TumorBam="../../Data/Exampledata_for_ASLOH/Example_TumorBam.bam" #Tumor Bam File


#Output
sOutputPrefix="Example"
LOHRegion="LOHRegion_"+sOutputPrefix #Output of LOH region
FilteredVCFbed="FilteredBed_"+sOutputPrefix #bed format output of filtered VCF 
LOHCandidateVariant="Candidate_"+sOutputPrefix #Intersected SNV with LOH region
ASLOHOutput="ASLOH_"+sOutputPrefix #



#Used Tool
Bedtools="bedtools" #Path of Bedtools
IGVtools="igvtools" #Path of igvtools



class SNV:

	def __init__(self):
		self.sChromosome=""
		self.nChromosome=""
		self.nPosition=""
		self.sReferenceAllele=""
		self.sAlternativeAllele=""
		self.sGenotype=""
		self.sVariantType=""
	
	def ParseSNV(self,sLine):
		
		sLine=sLine.strip()
		t=sLine.split("\t")
		
		
		(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
		sSample=t[-1]
		sGenotype=sSample.split(":")[0]
		
		self.sChromosome=sChr
		
		if len(sRef)==len(sAlt):
			self.sVariantType="SNP"
			self.sReferenceAllele=sRef
			self.sAlternativeAllele=sAlt
			self.sGenotype=sGenotype
			self.nPosition=int(nPosition)
		
		
		
		elif len(sRef)>len(sAlt):
			self.sVariantType="Del"
			sNewRef=sRef[len(sAlt):]
			sNewAlt="-"
			self.sReferenceAllele=sNewRef
			self.sAlternativeAllele=sNewAlt
			self.sGenotype=sGenotype
			self.nPosition=int(nPosition)+len(sAlt)
		elif len(sRef)<len(sAlt):
			self.sVariantType="Ins"
			sNewRef="-"
			sNewAlt=sAlt[len(sRef):]
			self.sReferenceAllele=sNewRef
			self.sAlternativeAllele=sNewAlt
			self.sGenotype=sGenotype
			self.nPosition=int(nPosition)
		
		
		
		
		
		
		
		
		
		if sChr=="X":
			self.nChromosome=23
		else:
			sChr=sChr.replace("chr","")
			self.nChromosome=int(sChr)
	
		





def Extracting_LOH_Region(InputCNVvcf):
	
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
				lCNVINFO=sINFO.split(";")
				(sCNVTYPE,nTCM_EM,nEndPosition)=(lCNVINFO[0],lCNVINFO[-3],lCNVINFO[2])
				(sCNVTYPE,nTCM_EM,nEndPosition)=(sCNVTYPE.split("=")[1],nTCM_EM.split("=")[1],nEndPosition.split("=")[1])
				if ((sCNVTYPE=="DUP-LOH") or (sCNVTYPE=="HEMIZYG") or (sCNVTYPE=="LOH")):
					fout.write("{0}\n".format("\t".join([t[0].replace("chr",""),t[1],nEndPosition,sCNVTYPE,nTCM_EM])))
				
				
				
			
			
			
	fout.close()
		
	

def LowQuality_Variant_Filtering(InputVCF):
	
	fp=gzip.open(InputVCF,"rt")
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
	
		fout.write("{0}\n".format("\t".join(map(str,[sVariant.sChromosome.replace("chr",""),\
		sVariant.nPosition,sVariant.nPosition,sVariant.sReferenceAllele,sVariant.sAlternativeAllele,sVariant.sGenotype]))))
	
	fout.close()
	
	
	

def Intersection_for_LOH_Candidate_Variant(FilteredVCFbed,LOHRegion):
	
	
	
	os.system(Bedtools+" intersect -a "+LOHRegion+" -b "+FilteredVCFbed+" -wo > "+LOHCandidateVariant)
	
	
	




def Igvtools(TumorBam,LOHCandidateVariant):
	
	
	
	try:
		os.mkdir("IGVtools/")
	except:
		pass
	
	fp=open(LOHCandidateVariant)
	
	sBamFile=TumorBam
	
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition,sRef,sAlt,sSS)=(t[5],t[6],t[8],t[9],t[10])
		if sSS=="0/1":
		
			sKey="_".join([sChr, nPosition,sRef,sAlt])
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
	(sChr,nPosition,sRef,sAlt)=sKey.split("_")
	try:
		(nPos, nA, nC, nG, nT, nN, nDel, nIns)=sLine.split("\t")
		dIGVdict[sKey]={"A":nA,"G":nG,"C":nC,"T":nT,"Del":nDel,"Ins":nIns}
	except ValueError:
		dIGVdict[sKey]="noRead"
	
	
	
	
	
	
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
	fout.write("Chr\tPosition\tRef\tAlt\tSS\tRefreadinTumor\tAltreadinTumor\tAARF\tAllele-specific_LOH\n")
	#for sLin
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition,sRef,sAlt,sSS)=(t[5],t[6],t[8],t[9],t[10])
		sKey="_".join([sChr, nPosition,sRef,sAlt])
		
		if sSS=="0/1":
			try:
				
				if not dIGVdict[sKey]=="noRead":
					
					if sAlt=="-":
						try:
							sAltinTumor=dIGVdict[sKey]["Del"]
							sAllreads=dIGVdict[sKey].values()
							sTotalRead=sum(map(float,sAllreads))
							nAARF=float(sAltinTumor)/float(sTotalRead)
							sRefinTumor=str(float(sTotalRead)-float(sAltinTumor))
							if nAARF<0.2:
								sSSwithCNV="RLOH"
							elif nAARF>0.8:
								sSSwithCNV="LOH"
							else:
								sSSwithCNV=sSS
							fout.write("{0}\t".format("\t".join([sChr, nPosition,sRef,sAlt,sSS])))
							fout.write("{0}\n".format("\t".join([sRefinTumor,sAltinTumor,str(nAARF),sSSwithCNV])))
						except ZeroDivisionError:
							pass	
					
						
					elif sRef=="-":
						try:
							sAltinTumor=dIGVdict[sKey]["Ins"]
							sAllreads=dIGVdict[sKey].values()
							sTotalRead=sum(map(float,sAllreads))
							nAARF=float(sAltinTumor)/float(sTotalRead)
							sRefinTumor=str(float(sTotalRead)-float(sAltinTumor))
							if nAARF<0.2:
								sSSwithCNV="RLOH"
							elif nAARF>0.8:
								sSSwithCNV="LOH"
							else:
								sSSwithCNV=sSS
							fout.write("{0}\t".format("\t".join([sChr, nPosition,sRef,sAlt,sSS])))
							fout.write("{0}\n".format("\t".join([sRefinTumor,sAltinTumor,str(nAARF),sSSwithCNV])))
						except ZeroDivisionError:
							pass	
							
							
							
					else:
						(sRefinTumor,sAltinTumor)=(dIGVdict[sKey][sRef],dIGVdict[sKey][sAlt])
						try:
							if float(sRefinTumor)==0:
								sAllreads=dIGVdict[sKey].values()
								sTotalRead=sum(map(float,sAllreads))
								nAARF=float(sAltinTumor)/float(sTotalRead)
							else:
								nAARF=float(sAltinTumor)/(float(sRefinTumor)+float(sAltinTumor))
							
							
							if nAARF<0.2:
								sSSwithCNV="RLOH"
							elif nAARF>0.8:
								sSSwithCNV="LOH"
							else:
								sSSwithCNV=sSS
							fout.write("{0}\t".format("\t".join([sChr, nPosition,sRef,sAlt,sSS])))
							fout.write("{0}\n".format("\t".join([sRefinTumor,sAltinTumor,str(nAARF),sSSwithCNV])))
						except ZeroDivisionError:
							pass	
					
					
					
				else:
					pass
			except KeyError:
				pass
		
		
		






if __name__=="__main__":
	
	Extracting_LOH_Region(InputCNVvcf)
	LowQuality_Variant_Filtering(InputVCF)
	Intersection_for_LOH_Candidate_Variant(FilteredVCFbed,LOHRegion)
	Allele_specific_LOH_Variant_Calling(TumorBam,LOHCandidateVariant)
	
	
	
	
	





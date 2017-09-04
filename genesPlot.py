import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pysam
import sys

from bioformat import GTF_Line

from genemodel import *
import bamio

#plt.rcParams['axes.labelpad']=10
plt.rcParams['axes.ymargin'] =0.2 
if len(sys.argv) == 4:
    depthplotflag = 1
    sys.stdout.write("Assume you supply a bam file to plot sequence coverage.\n")
else:
    depthplotflag = 0

f=file(sys.argv[1])  ## only support standard gtf format for now
geneOI = sys.argv[2].strip()  ## geneid or genename of interst

for line in f:
    gtfline = GTF_Line(line)
    if gtfline.gene_name() == geneOI or gtfline.gene_id() == geneOI: 
        pass
    else:
        continue
    chromosome = gtfline.seqname()
    start = gtfline.start()
    end = gtfline.end()
    strand = gtfline.strand()
    gene_id = gtfline.gene_id()
    gene_name = gtfline.gene_name()
    transcript_id = gtfline.transcript_id()
    transcript_name = gtfline.transcript_name()
    feature = gtfline.feature()
    if feature == 'gene':
        geneid = Gene(chromosome, start, end, strand, gene_id, name=gene_name,attr={'biotype':gtfline.gene_type()})
    elif feature == "transcript":
        transcriptobj = transcript(chromosome, start, end, strand,feature,transcript_id,{'biotype':gtfline.gene_type()})
        geneid.children[transcript_id]=transcriptobj
    elif feature == "CDS":
        cdsobj = CDS(chromosome, start, end, strand,transcript_id,{})
        geneid.children[transcript_id].addCDS(cdsobj)
    elif feature == "exon":
        exonobj = Exon(chromosome, start, end, strand,transcript_id,{})
        geneid.children[transcript_id].addexon(exonobj)
    elif feature == "UTR":
        utrobj = UTR(chromosome, start, end, strand,transcript_id,{})
        geneid.children[transcript_id].addUTR(utrobj)
    elif feature == "start_codon":
        cobj = StartCodon(chromosome, start, end, strand,transcript_id,{})
        geneid.children[transcript_id].addCodon(cobj)
    elif feature == "stop_codon":
        cobj = StopCodon(chromosome, start, end, strand,transcript_id,{})
        geneid.children[transcript_id].addCodon(cobj)
    else:
        sys.stderr.write("This is an offending line !!!! \n%s"%line)
f.close()

def axaddexons(ax,exons,ycoord,strand="+",color='orange'):
    ax=plt.gca()
    for s,e in exons:
        if strand == "+":
            ax.arrow(s,ycoord,e-s,0,fc=color,ec='gold', lw=0.05,width=0.5,head_width=0.5,head_length=0.3*(e-s),shape='full',length_includes_head=True)
        else:
            ax.arrow(e,ycoord,s-e,0,fc=color,ec='gold', lw=0.05,width=0.5,head_width=0.5,head_length=0.3*(e-s),shape='full',length_includes_head=True)
    return True

def formatxlabels(xcoord):
    return "%.2fkb"%(xcoord/float(1000))




plt.rcParams[u'ytick.major.pad']=10
#plt.tick_params('y',pad=20)   # don't work
fig=plt.figure(figsize=(20,16))
if depthplotflag:
    axdepth = fig.add_axes([0,0.05,1,0.125])
    ax = fig.add_axes([0.1,0.15,0.8,0.85],sharex=axdepth)
else:
    ax=fig.add_subplot(111)

txID = []
for i,transcript in enumerate(geneid.children.keys()):
    txID.append(transcript)
    geneid.children[transcript].inferUTRort()
    geneid.children[transcript].inferUTR()
    geneid.children[transcript].inferCodons()
    geneid.addisoform(geneid.children[transcript])
    axaddexons(ax,geneid.children[transcript].exons,i,strand=geneid.children[transcript].strand)
    axaddexons(ax,geneid.children[transcript].fp_utr,i,strand=geneid.children[transcript],color='darkviolet')
    axaddexons(ax,geneid.children[transcript].tp_utr,i,strand=geneid.children[transcript],color='red')
    ax.hlines(i ,geneid.start-0.1*geneid.length,geneid.end+0.1*geneid.length, linestyle=u'dotted', color='black', linewidth=1)


targetRegion = "-".join(map(str,[geneid.chromosome,geneid.start,geneid.end]))
if depthplotflag:
    samfile = sys.argv[3]
    bamfile = pysam.AlignmentFile(samfile,'rb')
    basecount,Chr_Region,depth_arr = bamio.singlebaseCov(bamfile,targetRegion)
    bamfile.close()    
    axdepth.fill_between(Chr_Region,depth_arr,alpha=0.4)
    axdepth.set_xlim(geneid.start-0.1*geneid.length,geneid.end+0.1*geneid.length)
    axdepth.set_ylim(bottom=0)
    axdepth.spines['top'].set_visible(False)
    axdepth.spines['right'].set_visible(False)
    axdepth.set_yscale('log')
    axdepth.set_ylabel('depth(log10)')
#ax.yaxis.labelpad =1000
ax.set_yticks(range(len(txID)))
ax.set_yticklabels(txID)

xticks = np.linspace(geneid.start-0.1*geneid.length,geneid.end+0.1*geneid.length,6)
xlabels = map(formatxlabels,xticks)
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_facecolor("snow")
#ax.set_axis_bgcolor("snow")
ax.set_xlim(geneid.start-0.1*geneid.length,geneid.end+0.1*geneid.length)
ax.set_ylim(bottom=-1,top=i+1)

plt.title(geneOI)
plt.tight_layout()

plt.savefig("%s.png"%geneOI,format='png',dpi=500)
plt.savefig('%s.svg'%geneOI,format='svg',dpi=500)
plt.clf()
plt.close()

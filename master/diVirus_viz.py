
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import scipy as sp
import matplotlib.lines as mlines

 ### #%% to create cells
############################################################################################################################################
#DEFINE FUNCTIONS and UNIVERSALLY-USED VARIABLES
############################################################################################################################################
seg_dict = {'NC_036476.1':'L1',
            'NC_036468.1':'L2',
            'NC_036477.1':'L3',
            'NC_036469.1':'M1',
            'NC_036470.1':'M2',
            'NC_036471.1':'M3',
            'NC_036472.1':'S1',
            'NC_036473.1':'S2',
            'NC_036474.1':'S3',
            'NC_036475.1':'S4'}

#for bubble plot legend
minColor = 0.03#np.min(shandiv_data)
maxColor = 0.9#np.max(shandiv_data)
sizeScaling = 500
legendBubbles = [0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8]
#[0.2,0.4,0.6,0.8,1.0]

def renameSegments(text):
    """
    Function: replace string with corresponding string in dict
    Variables: text in dict
    Return: corresponding text
    """
    return seg_dict[text]

def getBubbleLegend(legendBubbles,scatterInfo,sizeScaling,minColor,maxColor):
    """
    Function: create dummy plot to retrieve legend from and use on final plot
    Variables: legendBubbles=values in legend
               scatterInfo=information from official/ final scatterplot to be used in replica
               sizeScaling=integer by which bubbles will be scaled
               minColor=retrieve colour at min value
               maxColor=retrieve colour at max value
    Return: legend
    """
    plt.figure(num='dummy')
    for value in legendBubbles:
        thisColor = (value-minColor)/(maxColor-minColor)
        plt.plot([0],[0],'o',color=scatterInfo.cmap(thisColor),markersize=np.sqrt(value*sizeScaling),label='%.2f' %(value),alpha=0.5)
    plt.legend()
    ax_tmp = plt.gca()
    handles, labels = ax_tmp.get_legend_handles_labels()
    plt.close('dummy')
    return handles,labels

def piSites(df):
    """
    Function: interleaves same number of sites between existing positions in file with pi information
    for each segment.
    Variable: df is a data frame constructed from the nucleotide diversity files
    Return: sites left joined with respective pi files and NaN values replaced with 0
    """
    sites_df=pd.DataFrame(np.arange(1, longest_seg+1))
    sites_df.columns=['site']
    sites_df=sites_df.site.astype(object)
    tempidf=pd.merge(sites_df, df, how='left', on=['site']).fillna(0)
    for elem in seg_names:
        if elem in tempidf['product'].tolist():
            tempidf['product']=elem
    return tempidf

def afSites(df):
    """
    Function: sums frequencies at same position then interleaves same number of
    sites between existing positions in file with alt allele freq information
    for each segment.
    Variable: df is a data frame constructed from the alt allele freq files
    Return: sites (longest segment x number of segments in viral genome left
    joined with respective freq files and NaN values replaced with 0
    """
    # for i in range(0,len(df)):
    #     if (df.iloc[i].site==df.iloc[i-1].site):
    #         df.loc[i-1,'freq'] = (df.iloc[i].freq) + (df.iloc[i-1].freq)
    #         df.drop([i], inplace=True)
    #         df.reset_index(drop=True, inplace=True)
    #     elif (df.iloc[i].site==df.iloc[i-1].site):
    #         df.loc[i-1,'freq'] = (df.iloc[i-1].freq) + (df.iloc[i].freq)
    #         df.drop([i], inplace=True)
    #         df.reset_index(drop=True, inplace=True)
    df=df.groupby(['product', 'site'])['freq'].agg('sum')
    df=df.reset_index()
    df['site']=df['site'].astype('str').astype('int') #if Pos column not converted to int, treated as str and sorted lexicographically (i.e. 1,10,100,100,2,20,200,etc..)
    df['product']=df['product'].apply(renameSegments)
    df=df.sort_values(['product', 'site'], ascending=[True, True])
    siteRows=list(range(1,longest_seg+1))*seg_num
    for i in seg_dict.values():
        segRows.extend([i]*longest_seg)
    df2merge=pd.DataFrame(list(zip(segRows, siteRows)),
               columns=['product', 'site'])
    temp_freq=pd.merge(df2merge, df, how='left', on=['product','site'], left_index=False, right_index=False).fillna(0)
    return temp_freq


###########################################################################################################################################
#COVERAGE PLOTS FOR ORF ONLY/ cov files made using bed with specific orf regions only#
###########################################################################################################################################

#load in coverage data for each file - just ORFs

cov_data=[]
row_pos=[0]*11

#exclude the Gblock/ VB seqs. they will inhibit plotting the rest
for file in glob.glob('*_.sort.cov'):
        df=pd.read_csv(file, delim_whitespace=True, header=None)
        df=df.sort_values([0, 1], ascending=[True, True])
        df.reset_index(drop=True, inplace=True)
        df_list = file, df
        cov_data.append(df_list)

orf_len=df.groupby([0]).count()
orf_len=orf_len[1].values.tolist()
orf_len=list(map(int, orf_len))
#seg_len=seg_len.astype('str', inplace=True)

for i in range(1, 11):
	row_pos[i]=orf_len[i-1]+row_pos[i-1]


#segments L1-L3,M1-M3,S1-S4
#figure number, size was 100, 40
fig = plt.figure(110, figsize=(100, 40))
n = len(cov_data)
#c_map = plt.cm.get_cmap('viridis',n)
#cnorm = PowerNorm(0.4, vmin=0, vmax=50)
viridis = plt.cm.get_cmap('viridis', n)
colors=viridis(np.linspace(1,0,n))
seg_names=('L1','L2','L3','M1','M2','M3','S1','S2','S3','S4')
order=(8,0,9,1,2,3,4,5,6,7)
#fig.text(0.5, 0.05, 'Genomic Position (nt)', ha='center', fontsize='80')
#fig.text(0.01, 0.5, 'Sequencing Depth (reads)', va='center', rotation='vertical', fontsize='80')
# for j in range(0, len(cov_data)):
#         for i in range(1,len(order)+1):
#                 plt.subplot(2, 5, i).set_title('Segment '+seg_names[i-1], fontdict={'size':65})
#                 plt.subplots_adjust(left=0.05, wspace=0.15, hspace=0.2)
#                 plt.tick_params('both', labelsize=45)
#                 xi=row_pos[order[i-1]]
#                 xf=row_pos[order[i-1]]+orf_len[order[i-1]]-1
#                 plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color=colors[j], linewidth=4, label=cov_data[j][0])
#                 plt.axhline(y=100, color='r', linestyle='--',  dashes=(3, 6), linewidth=5)
#                 #plt.ylim(1, 8050)
#                 plt.ylim(1, 10000)
#                 plt.xlim(0)
#                 plt.yscale('log')
#         leg=plt.legend(bbox_to_anchor=(1.04,2), loc="upper left", fontsize=40, frameon=False)
#         #fig.tight_layout()
#         for line in leg.get_lines():
#                 line.set_linewidth(6.0)
#                 plt.savefig('110_plot.png') #may wish to add bbox_inches='tight'

# gold100=mlines.Line2D([], [], color='gold', linewidth=0.9,
#                           label='samples with regions of < 100 reads')
# teal500=mlines.Line2D([], [], color='teal', linewidth=0.9,
#                           label='samples with regions of < 500 reads')
# indigoAbove=mlines.Line2D([], [], color='indigo', linewidth=0.9,
#                           label='samples with regions of > 500 reads')
fig = plt.figure(12)
fig.text(0.44, 0.02, 'ORF Position (nt)', fontsize=14)
fig.text(0.005, 0.3, 'Coverage log(reads)', rotation=90, fontsize=14)
for j in range(0, len(cov_data)):
        for i in range(1,len(order)+1):
            plt.subplot(2, 5, i).set_title(seg_names[i-1], fontdict={'size':12})
            plt.subplots_adjust(left=0.05, wspace=0.15, hspace=0.3)
            plt.tick_params('both', labelsize=10)
            xi=row_pos[order[i-1]]
            xf=row_pos[order[i-1]]+orf_len[order[i-1]]-1
            if np.any(np.asarray(cov_data[j][1][2][list(range(xi+5, xf-5))]) < 100):
                plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color='gold', linewidth=0.7, alpha=0.5)
            elif np.any(np.asarray(cov_data[j][1][2][list(range(xi+5, xf-5))]) < 500):
                plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color='teal', linewidth=0.7, alpha=0.5)
            else:
                plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color='indigo', alpha=0.6, linewidth=0.7)
            plt.axhline(y=100, color='gold', linestyle='--', dashes=(3, 6), linewidth=0.7)
            plt.axhline(y=500, color='teal', linestyle='--', dashes=(3, 6), linewidth=0.7)
            #plt.ylim(1, 8050)
            plt.ylim(1, 1000000)
            plt.xlim(0)
            plt.yscale('log')
            ax=plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
       # leg2=plt.legend(handles=[indigoAbove, teal500, gold100], title='Read Depth Threshold', frameon=False, bbox_to_anchor=(1.13,0.45), loc='center right', labelspacing=1.12)
       # ax = plt.gca().add_artist(leg2)




Reds=plt.cm.get_cmap('Reds', n)
reds=Reds(np.linspace(0,1,n))
Magma=plt.cm.get_cmap('magma', n)
magma=Magma(np.linspace(0,1,n))
#individual segment plots with samples that have sites <n reads highlighted in red
for i in range(1, len(order)+1):
    plt.figure()
    for j in range(0, len(cov_data)):
        plt.title('Segment '+seg_names[i-1]+ 'ORF start+100, ORF stop-100', fontsize=16)
        plt.xlabel('Genomic Position (nt)', fontsize=12)
        plt.ylabel('Coverage log(reads)', fontsize=12)
        plt.tick_params('both', labelsize=10)
        xi=row_pos[order[i-1]]
        xf=row_pos[order[i-1]]+orf_len[order[i-1]]-1
        if np.any(np.asarray(cov_data[j][1][2][list(range(xi+500, xf-500))]) < 100):
            plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color='green', linewidth=0.9, label=cov_data[j][0])
        elif np.any(np.asarray(cov_data[j][1][2][list(range(xi+500, xf-500))]) < 500):
            plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color=magma[j], linewidth=0.9, label=cov_data[j][0])
        else:
            plt.plot(cov_data[j][1][1][list(range(xi, xf))], cov_data[j][1][2][list(range(xi, xf))], color='navy', linewidth=0.9)
        plt.axhline(y=100, color='r', linestyle='--', dashes=(3, 6), linewidth=0.8)
        plt.axhline(y=500, color='orange', linestyle='--', dashes=(3, 6), linewidth=0.8)
        #dashes=(3, 6)
        #plt.ylim(1, 8050)
        plt.ylim(1, 1000000)
        plt.xlim(0)
        plt.yscale('log')
        leg=plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=8, frameon=False)
        plt.tight_layout()
        for line in leg.get_lines():
            line.set_linewidth(0.9)
    ax=plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

#between runs or it will append to previous figure (even with diff name..) and be a hot mess
plt.clf()
plt.close()
fig.clf()
###########################################################################################################################################
#INTERHOST SNPS - ORF ONLY#
###########################################################################################################################################
#load breakdown of interhost snps by category
interNon=pd.read_csv('/Users/Jessica/Documents/MASTERS/PRV_Quasi/THESIS_RESULTS/interhost_aa_snp_sites_breakdown.csv')

#fill in sites without variation
aaSites=pd.DataFrame(np.arange(1,7487+124))
aaSites.columns=['aa_pos']
aaSites=aaSites.aa_pos.astype(object)
interNonsyn=pd.merge(aaSites, interNon, how='left', on=['aa_pos']).fillna(0)
reorder=['aa_pos','snp_count','farm_chin','w_15_chin','coho','vent_atl','other_wild_chin','other_atl']
interNonsyn=interNonsyn[reorder]
colors=['dimgrey','deeppink','orange','gold','mediumblue','tomato','lightskyblue']
categoryTotals=[0.01,102,15,15,5,31,19,17]
interNonsynPerc=interNonsyn/categoryTotals*100

#regular counts
fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7,1, sharex='all')

for i,c in zip(range(1,7+1), colors):
    plt.subplot(7,1,i)
    plt.plot(interNonsyn.aa_pos, interNonsyn.iloc[:,i], linestyle='-', color=c)
    ax=plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yscale('log')
    #ax.set_ylim(0, 40)

ax1=plt.subplot(7,1,1)
#ax1.axhline(40, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax2=plt.subplot(7,1,2)
ax2.axhline(15, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax3=plt.subplot(7,1,3)
ax3.axhline(15, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax4=plt.subplot(7,1,4)
ax4.axhline(5, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax4.set_ylabel('log(Number of Nonsynonymous SNPs)', fontsize=12)
ax5=plt.subplot(7,1,5)
ax5.axhline(31, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax6=plt.subplot(7,1,6)
ax6.axhline(19, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax7=plt.subplot(7,1,7)
ax7.axhline(17, lw=0.5, linestyle='--', color='black', alpha=0.4)
ax7.set_xlabel('Genomic Position (aa)', fontsize=12)

newax=ax7.twiny()
newax.set_frame_on(True)
newax.patch.set_visible(False)
newax.xaxis.set_ticks_position('bottom')
newax.xaxis.set_label_position('bottom')
newax.spines['bottom'].set_position(('outward', 40))
newax.xaxis.set_ticks([-379.45, 0, 1283, 2574, 3861, 4622, 5310, 6063, 6394,
                       6518,6939,7294, 7610, 7990.45])
newax.set_xticklabels('')


#proportion of each category with SNP
# fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7,1, sharex='all')

# for i,c in zip(range(1,7+1), colors):
#     plt.subplot(7,1,i)
#     plt.plot(interNonsyn.aa_pos, interNonsynPerc.iloc[:,i], linestyle='-', color=c)
#     ax=plt.gca()
#     ax.legend([''], frameon=False, loc='upper right')
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)

plt.figure()
plt.ylabel('Number of Nonsynonymous SNPs')
plt.xlabel('Genomic Position (aa)')

newax=ax7.twiny()
newax.set_frame_on(True)
newax.patch.set_visible(False)
newax.xaxis.set_ticks_position('bottom')
newax.xaxis.set_label_position('bottom')
newax.spines['bottom'].set_position(('outward', 40))
newax.xaxis.set_ticks([-370, 0,1283, 2574, 3861, 4622, 5310, 6063, 6394,
                       6815, 7170, 7486, 7600])
newax.set_xticklabels('')




###########################################################################################################################################
#NUCLEOTIDE DIVERSITY BUBBLE PLOTS - ORF ONLY#
###########################################################################################################################################


longest_seg=3934
seg_num= 11
pi_filenames=[]
pidf=[]
tempidf=[]
final=[]
final_df=[]
pi_prelim=[]
pi_data=[]
samples=[]
col_names=[]
temp_pidf=[]
pi_bubble=[]
sample_names=[]
seg_names=('L1','L2','L3','M1','M2','M3','S1','S1b','S2','S3','S4')
start=[19, 19, 8, 22, 27, 84, 29, 108, 22, 29, 39]
stop=[3867, 3891,  3868, 2304, 2090, 2342, 1021, 482, 1284, 1093, 986]
scatterInfo=[]


#read in pi data files
for pi_file in glob.glob('*_ntdiv.txt'):
    pi_filenames.append(pi_file)
    pi_filenames.sort()

for pi_file in pi_filenames:
    pidf=pd.read_csv(pi_file, sep='\t', header=0)
    pidf.drop(pidf[pidf['product'] == 'noncoding'].index, inplace=True)
    pidf['class_vs_ref'].replace(['NONREF_NONPOLY'], 0, regex=True, inplace=True)
    tempidf=piSites(pidf)
    final.append(tempidf)
    final_df=pd.concat(final, axis=0)

#create df containing pi values (col=samples, rows=sites)..*note: every segment is filled to the longest segment with 0's
pi_prelim=pd.DataFrame(np.asarray(final_df).reshape(longest_seg*seg_num, -1, order='F'))
pi_prelim.drop(pi_prelim.iloc[:,1:int(len(pi_prelim.columns)/2)-1], axis=1, inplace=True)

#create column names for subsequent df's
samples=pi_filenames[0::seg_num]

for name in samples:
    new_file=name.split("-")[0]
    col_names.append(new_file)


col_names.insert(0,'product')
col_names.insert(0,'site')
col_names[52]='NGS268'

sample_names=col_names[2:len(col_names)]


#df containing info on whether SNV causes syn or nonsyn mutation
pi_aa=pi_prelim.drop(pi_prelim.iloc[:, 2:int(len(pi_prelim.columns)/2+1)], axis=1, inplace=False)
#pi_aa=pi_aa.rename(columns={pi_aa.columns[0]: "site", pi_aa.columns[1]: "product"})
pi_aa.columns=col_names

#eliminate the aa change info from the last columns now that they have been transferred to pi_aa df
pi_data=pi_prelim.drop(pi_prelim.iloc[:, int(len(pi_prelim.columns)/2+1):int(len(pi_prelim.columns))], axis=1, inplace=False)
#pi_data=pi_data.rename(columns={pi_data.columns[0]: "site", pi_data.columns[1]: "product"})
pi_data.columns=col_names

seg_ihDi=np.zeros(int(len(pi_filenames)/seg_num))
synCount=np.zeros(int(len(pi_filenames)/seg_num))
nonsynCount=np.zeros(int(len(pi_filenames)/seg_num))
ambigCount=np.zeros(int(len(pi_filenames)/seg_num))
genPi=np.zeros(int(len(pi_filenames)/seg_num))
piSegSnvCounts=np.zeros(int(len(pi_filenames)/seg_num))
GenomeSyn=np.zeros(int(len(pi_filenames)/seg_num))
GenomeNonsyn=np.zeros(int(len(pi_filenames)/seg_num))
GenomeAmbiguous=np.zeros(int(len(pi_filenames)/seg_num))
piGenome=np.zeros(int(len(pi_filenames)/seg_num))
#allMuts=np.zeros(int(len(pi_filenames)/seg_num))
total_nt=0
#colNums=np.asarray(range(204, 306))

phylOrder=['NGS140','NGS108', 'NGS195','NGS265','NGS127','NGS138',
           'NGS124','NGS133','NGS260','NGS84','NGS261','NGS264','NGS131',
           'NGS263','NGS129','NGS162','NGS309','NGS145','NGS146','NGS190',
           'NGS171','NGS181','NGS167','NGS166','NGS104','NGS164','NGS186',
           'NGS134','NGS103','NGS99','NGS174','NGS194','NGS144','NGS107',
           'NGS329','NGS327','S1003','NGS328','S1053','NGS330','H1056',
           'NGS76','NGS326','NGS312','NGS311','NGS314','NGS313','NGS317',
           'NGS251','NGS207','NGS211','NGS213','NGS200','NGS212', 'NGS249',
           'NGS176','NGS325','NGS318','NGS319', 'NGS321','NGS280','NGS302',
           'NGS282','NGS315','NGS266','NGS281','NGS301','NGS73','NGS324',
           'NGS322','NGS150','NGS243','NGS71','NGS74','NGS323','NGS235',
           'NGS231','NGS75','NGS292','NGS271','NGS272','NGS273','NGS274',
           'NGS275','NGS276','NGS277','NGS278','NGS279','NGS286','NGS269',
           'NGS291','NGS268','NGS227','NGS285','NGS234','NGS226','NGS225',
           'NGS294','NGS270','NGS70','NGS283','NGS267']

#legend for marker shape

triangle = mlines.Line2D([], [], color='grey', marker='^', linestyle='None',
                          markersize=10, label='Nonsynonymous')
circle = mlines.Line2D([], [], color='grey', marker='o', linestyle='None',
                          markersize=10, label='Synonymous')
square = mlines.Line2D([], [], color='grey', marker='s', linestyle='None',
                          markersize=10, label='Ambiguous')
#plot
for xi,xf,seg in zip(start,stop,seg_names):
    temp_pidf=pi_data.loc[(pi_data['product'] == seg) & (pi_data['site'] > xi-1) & (pi_data['site'] < xf+1)]
    temp_pidf.drop(columns=['product','site'], axis=1, inplace=True)
    #reorder the way you want
    temp_pidf=temp_pidf[phylOrder]
    temp_pidf.columns=np.arange(len(temp_pidf.columns))
    #geometric mean
    pi_array=temp_pidf.to_numpy(dtype='float64')
    mask_pi_array=np.ma.masked_where(pi_array==0, pi_array)
    pi_geom_mean=sp.stats.mstats.gmean(mask_pi_array, axis=0)
    genPi=np.vstack((genPi, pi_array))
    #put ORFs with individual site pi values into df
    piGenome=np.vstack((piGenome, pi_array))
    #number sites with snv in orf
    snv_count=np.ma.MaskedArray.count(mask_pi_array, axis=0)
    #each row = seg, col=sample
    piSegSnvCounts=np.vstack((piSegSnvCounts, snv_count.data))
    #total number nt in seg
    total_seg_nt=len(mask_pi_array)
    #total nt in genome ORFs
    total_nt=total_nt+total_seg_nt
    #orf ihDi
    ihDi=(snv_count/total_seg_nt)*pi_geom_mean
    #array containing values (col=sample, rows=seg)
    seg_ihDi=np.vstack((seg_ihDi, ihDi.data))
    pi_bubble = temp_pidf.unstack().reset_index()
    pi_bubble.columns = list("YXZ")
    pi_bubble.Z=pi_bubble.Z.astype('float')
    temp_aadf=pi_aa.loc[(pi_aa['product'] == seg) & (pi_aa['site'] > xi-1) & (pi_aa['site'] < xf+1)]
    temp_aadf.drop(columns=['product','site'], axis=1, inplace=True)
    #reorder the way you want
    temp_aadf=temp_aadf[phylOrder]
    #assign int values as column names in place of string; will allow plotting
    #in consistent, predefined order
    temp_aadf.columns=np.arange(len(temp_aadf.columns))
    #convert to array for counting # syn vs. non
    aa_array=temp_aadf.to_numpy(dtype='str')
    #row=orf, column=sample, value= count
    segSyn= np.count_nonzero(aa_array=='Synonymous',axis=0)
    GenomeSyn=np.vstack((GenomeSyn, segSyn ))
    normsegSyn=segSyn/total_seg_nt
    #normsegSyn=segSyn/1000
    segNonsyn= np.count_nonzero(aa_array=='Nonsynonymous',axis=0)
    GenomeNonsyn=np.vstack((GenomeNonsyn, segNonsyn ))
    normsegNonsyn=segNonsyn/total_seg_nt
    segAmbiguous=np.count_nonzero(aa_array=='Ambiguous',axis=0)
    GenomeAmbiguous=np.vstack((GenomeAmbiguous, segAmbiguous ))
    normsegAmbiguous=segAmbiguous/total_seg_nt
    #normsegNonsyn=segNonsyn/1000
    #all orfs stacked in array
    synCount=np.vstack((synCount, normsegSyn))
    nonsynCount=np.vstack((nonsynCount,normsegNonsyn))
    ambigCount=np.vstack((ambigCount,normsegAmbiguous))
    aa_bubble = temp_aadf.unstack().reset_index()
    aa_bubble.columns = list("yxz")
    #pi_bubble['samp_ord']=pd.Categorical(pi_bubble['Y'], sample_order)
    #pi_bubble.sort_values(by='samp_ord')
    #aa_bubble['samp_ord']=pd.Categorical(aa_bubble['y'], sample_order)
    #aa_bubble.sort_values(by='samp_ord')
    plt.figure(figsize=(100,150))
    for mutation, mark in zip(['Nonsynonymous','Synonymous','Ambiguous'], ['^','o','s']):
        markevery=aa_bubble.index[aa_bubble.z==mutation].tolist()
        scatterInfo=plt.scatter(x=pi_bubble['X'].iloc[markevery], y=pi_bubble['Y'].iloc[markevery],
                                marker=mark, s=pi_bubble['Z'].iloc[markevery]*sizeScaling,
                                c=pi_bubble['Z'].iloc[markevery], vmin=minColor,vmax=maxColor,
                                cmap="viridis", alpha=0.5)
    leg2=plt.legend(handles=[triangle, circle, square], title='SNV Type', frameon=False,
                    bbox_to_anchor=(1.13,0.45), loc='center right', labelspacing=1.12)
    ax = plt.gca().add_artist(leg2)
    handles,labels = getBubbleLegend(legendBubbles,scatterInfo,sizeScaling,minColor,maxColor)
    plt.legend(handles,labels,title='Nucleotide (\u03C0) Diversity', frameon=False,
               bbox_to_anchor=(1.13,0.8), loc='center right', labelspacing=1.12)
    plt.xticks(np.linspace(np.min(pi_bubble['X']),np.max(pi_bubble['X']),5),np.linspace(xi,xf+1,5).astype(int))
    plt.xlim([np.min(pi_bubble['X'])-10,np.max(pi_bubble['X'])+10])
    allZeros=np.asarray(temp_pidf.loc[:, (temp_pidf == 0).all()].columns)
    #unnecessary: nonZeros=np.setxor1d(colNums, allZeros)
    #names
    indices1 = np.where(np.in1d(np.arange(len(temp_pidf.columns)), allZeros))[0]
    #indices1 = np.where(np.in1d(colNums, allZeros))[0]
    #accessed_mapping = map(sample_names.__getitem__, indices1)
    #allZeros_names=list(accessed_mapping)
    plt.yticks(ticks=np.arange(len(temp_pidf.columns)), labels=phylOrder, color='black', fontsize=5)
    #plt.yticks(ticks=colNums, labels=sample_names, color='black', fontsize=5)
    #plt.axhline(y=allZeros, color='red', linestyle='--', alpha=0.3, linewidth=1)
    plt.title('Nucleotide (\u03C0) Diversity per Position of Segment %s ORF for PRV Isolated from Salmon' %seg, fontsize=16)
    plt.ylabel('Salmon Host ID', fontsize=14)
    plt.xlabel('Genomic Position (nt)', fontsize=14, labelpad=10.0)
    ax1 = plt.gca()
    for y in indices1:
        ax1.get_yticklabels()[y].set_color("red")
    #for tick in ax1.yaxis.get_major_ticks()[1::2]:
        #tick.set_pad(15)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    #plt.box(False)


#delete 0th row in syn/ nonsyn counts
synCount=np.delete(synCount, 0, 0)
nonsynCount=np.delete(nonsynCount, 0, 0)
ambigCount=np.delete(ambigCount, 0, 0)
genPi=np.delete(genPi, 0, 0)
GenomeNonsyn=np.delete(GenomeNonsyn, 0, 0)
GenomeSyn=np.delete(GenomeSyn, 0, 0)
GenomeAmbiguous=np.delete(GenomeAmbiguous, 0, 0)
piGenome=np.delete(piGenome, 0, 0)

#how many variable sites
sitesPerGen=pi_data.astype(bool).sum(axis=0)
sitesPerGen=sitesPerGen[2:len(sitesPerGen)]
avgsitesPerGen=sitesPerGen/total_nt
percVariableGen=avgsitesPerGen*100


#sum of total non vs. syn
np.sum(GenomeNonsyn)
np.sum(GenomeSyn)
np.sum(GenomeAmbiguous)

#Geometric mean and Intrahost diversity index for entire genome (ORFs) calc from pi values
total_pi_snv=np.sum(piSegSnvCounts, axis=0)
normPiSnv=total_pi_snv/total_nt
maskGenPi=np.ma.masked_where(genPi==0, genPi)
full_piGeomean=sp.stats.mstats.gmean(maskGenPi, axis=0)
piGenIhdi=normPiSnv*full_piGeomean

#total syn nonsyn counts in genome (normalized by genome concatenated orf length)
normGenomeSyn=np.sum(GenomeSyn, axis=0)/total_nt
normGenomeNonsyn=np.sum(GenomeNonsyn, axis=0)/total_nt
normGenomeAmbig=np.sum(GenomeAmbiguous, axis=0)/total_nt

#is data normally distributed?

#histogram to view intrahost distribution amongst fish per segment
#rows=segments, col=salmon samples
#skip row 0 because this was empty array used to initialize for subsequent vstack
plt.figure()
for seg, color in zip(range(1,len(seg_ihDi)), ['red', 'blue','green','purple',
                                              'magenta','grey','orange','yellow',
                                              'pink','aqua','black']):
    #plt.figure()
    plt.hist(seg_ihDi[seg], bins=9, align='left', orientation='vertical',
              alpha=0.65, edgecolor=color, color=color, linewidth=1.2)
    plt.xlabel('Intrahost Diversity Index')
    plt.ylabel('Number of Salmon')
    plt.legend()


#Goodness of fit test:
alpha=0.05
for i, seg in zip(range(len(synCount)), seg_names):
    synStat, synP= sp.stats.shapiro(synCount[i])
    print('(', synStat,',', synP,'):')
    if synP > alpha:
        print('Synonymous samples for segment %s look Gaussian (fail to reject H0). PARAMETRIC DATA.\n' %seg)
    else:
        print('Synonymous samples for segment %s don\'t look Gaussian (reject H0). NONPARAMETRIC DATA.\n' %seg)
    nonsynStat, nonsynP= sp.stats.shapiro(nonsynCount[i])
    print('(', nonsynStat,',', nonsynP,'):')
    if nonsynP > alpha:
        print('Nonsynonymous samples for segment %s look Gaussian (fail to reject H0). PARAMETRIC DATA.\n' %seg)
    else:
        print('Nonsynonymous samples for segment %s don\'t look Gaussian (reject H0). NONPARAMETRIC DATA.\n' %seg)
#if alpha is 0.05...
#H0: the data and Gaussian ditr are the same - they're both normal.
#p <= alpha: reject H0, NOT NORMAL. NONPARAMETRIC.
#p > alpha: fail to reject H0, NORMAL. PARAMETRIC.

#syn/ nonsyn bar plot for PARAMETRIC distribution
# plt.figure()
# width = 0.35
# x=np.arange(len(seg_names))
# #error bars
# synErr95=np.std(synCount, axis=1)/np.sqrt(int(len(synCount[0])))*1.96
# nonsynErr95=np.std(nonsynCount, axis=1)/np.sqrt(int(len(nonsynCount[0])))*1.96
# #remove 1.96 to get standard error of mean (SE) values intead of 95% confidence interval (CI)
# plt.bar(x - width/2, np.mean(synCount, axis=1).T, width, label='Synonymous', color='lightgray', yerr=synErr95, capsize=2)
# plt.bar(x + width/2, np.mean(nonsynCount, axis=1).T, width, label='Nonsynonymous', color='dimgray', yerr=nonsynErr95, capsize=2)
# plt.xticks(range(0, seg_num+1), seg_names,fontsize=12.5)
# plt.xlabel('ORF', fontsize=14)
# plt.ylabel('Number SNPs/ Kb', fontsize=14)
# plt.title('Normalized Number of Nonsynonymous and Synonymous SNPs per Open Reading Frame in PRV', fontsize=16)
# plt.legend(title='Mutation Type', frameon=False, bbox_to_anchor=(1.13,0.8), loc='center right', labelspacing=1.12)
# ax=plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)



#for NONPARAMETRIC data
syn_df=[]
non_df=[]
#ambig_df=[]
syn_non_df=[]

#transform skewed data log(x+0.001)
tran_synCount=np.log(synCount+0.001)
tran_nonsynCount=np.log(nonsynCount+0.001)


#df for swarmplot and pointplot
syn_df=pd.DataFrame(tran_synCount.T)
syn_df.columns=list(seg_names)
syn_df['mut_type']='Synonymous'
nonsyn_df=pd.DataFrame(tran_nonsynCount.T)
nonsyn_df.columns=list(seg_names)
nonsyn_df['mut_type']='Nonsynonymous'
#ambig_df=pd.DataFrame(ambigCount.T)
#ambig_df.columns=list(seg_names)
#ambig_df['mut_type']='Ambiguous'
#non_syn_df=pd.concat([syn_df, nonsyn_df, ambig_df])
non_syn_df=pd.concat([syn_df, nonsyn_df])
long_non_syn_df=pd.melt(non_syn_df, id_vars='mut_type')

#syn/ nonsyn plot for NONPARAMETRIC DATA:
# plt.figure()
# sns.swarmplot(x=long_non_syn_df['variable'], y=long_non_syn_df['value'], hue=long_non_syn_df['mut_type'], palette={'Synonymous': 'g', 'Nonsynonymous':'blue'}, marker='o', dodge=0.45, alpha=0.25)
# sns.pointplot(x=long_non_syn_df['variable'], y=long_non_syn_df['value'], hue=long_non_syn_df['mut_type'], palette={'Synonymous': 'g', 'Nonsynonymous':'blue'}, size=4, ci=95, n_boot=1000, estimator=np.median, join=True, linestyles=['--', ':'], scale=0.5, dodge=0.4, errwidth=0.8, capsize=0.035)
# plt.xticks(range(0, seg_num+1),['Core Shell','Capping Enzyme', 'RdRp','Core NTPase',
#                                 'Outer Shell', 'NS-Factory','Outer Clamp','Cytotoxic Protein',
#                                 'Core Clamp', 'NS-RNA Protein', 'Attachment Fiber'], rotation=55, ha='right', fontsize=12)
# plt.xlabel('ORF', fontsize=14)
# plt.ylabel('Number SNPs/ Kb', fontsize=14)
# plt.title('Normalized Number of Nonsynonymous and Synonymous SNPs per Open Reading Frame in PRV', fontsize=16)
# plt.legend(title='Mutation Type', frameon=False, bbox_to_anchor=(1,0.8), loc='center right', labelspacing=1.12)
# ax=plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

#box and stripplot
# plt.figure()
# sns.stripplot(x='variable',y='value',data=long_non_syn_df, hue='mut_type', palette={'Synonymous': 'darkturquoise', 'Nonsynonymous':'tomato'}, alpha=0.25, dodge=True, jitter=0.3)
# ax1=sns.boxplot(x='variable',y='value',data=long_non_syn_df,hue='mut_type', color='white', showfliers = False)
# plt.xticks(range(0, seg_num+1),['Core Shell','Capping Enzyme', 'RdRp','Core NTPase',
#                                 'Outer Shell', 'NS-Factory','Outer Clamp','Cytotoxic Protein',
#                                 'Core Clamp', 'NS-RNA Protein', 'Attachment Fiber'], rotation=55, ha='right', fontsize=12)
# for i,box in enumerate(ax1.artists):
#     box.set_edgecolor('black')
#     box.set_facecolor('white')
#     for j in range(i*5,i*5+5):
#         line = ax1.lines[j]
#         line.set_color('black')
#         line.set_mfc('black')
#         line.set_mec('black')
# plt.xlabel('ORF', fontsize=14)
# plt.ylabel('Mean SNP Number', fontsize=14)
# plt.title('Mean Nonsynonymous and Synonymous SNPs per Open Reading Frame in PRV', fontsize=16)
# ax=plt.gca()
# handles, labels = ax.get_legend_handles_labels()
# l = plt.legend(handles[2:4], labels[2:4], bbox_to_anchor=(1,0.8), loc='center right', title='Mutation Type', frameon=False)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

#USE THIS syn/ nonsyn plot for NONPARAMETRIC DATA:
plt.figure()
ax=sns.boxenplot(x='variable',y='value',data=long_non_syn_df,hue='mut_type', palette={'Synonymous': 'darkturquoise', 'Nonsynonymous':'tomato'})
sns.stripplot(x='variable',y='value',data=long_non_syn_df, hue='mut_type', palette={'Synonymous': 'darkturquoise', 'Nonsynonymous':'tomato'}, size=3, linewidth=0.7, edgecolor='black', alpha=0.2, dodge=True, jitter=0.3)
for l in ax.lines:
    l.set_linestyle('-')
    l.set_color('black')
    l.set_alpha(1)
plt.xticks(range(0, seg_num+1),['L1','L2', 'L3','M1',
                                'M2', 'M3','S1','S1b',
                                'S2', 'S3', 'S4'], rotation=0, fontsize=12)
#['L1:Core Shell','L2:Capping Enzyme', 'L3:RdRp','M1:Core NTPase',
#                                'M2:Outer Shell', 'M3:NS-Factory','S1:Outer Clamp','S1b:Cytotoxic Protein',
#                                'S2:Core Clamp', 'S3:NS-RNA Protein', 'S4:Attachment Fiber'], rotation=55, ha='right'
plt.xlabel('ORF', fontsize=14)
plt.ylabel('ln(Mean SNP Number + 0.001)', fontsize=14)
plt.title('Mean Nonsynonymous and Synonymous SNPs per Open Reading Frame in PRV', fontsize=16)
ax1=plt.gca()
handles, labels = ax1.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], bbox_to_anchor=(1.1,0.8), loc='center right', title='Mutation Type', frameon=False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

#nonparametric significance tests
#Friedman test (non-parametric ANOVA) ranking segments to see if there exists
#significant difference in normalized syn counts or norm non counts
#between any, but not which is most/ least, etc...for this you
#require post hoc analysis
sp.stats.kruskal(synCount[0], synCount[1],  synCount[2],  synCount[3],
                            synCount[4],  synCount[5], synCount[6], synCount[7],
                            synCount[8], synCount[9], synCount[10] )
##KruskalResult(statistic=22.54081066658539, pvalue=0.012574488129942078)
sp.stats.kruskal(nonsynCount[0], nonsynCount[1],  nonsynCount[2],  nonsynCount[3],
                            nonsynCount[4],  nonsynCount[5], nonsynCount[6], nonsynCount[7],
                            nonsynCount[8], nonsynCount[9], nonsynCount[10] )
##KruskalResult(statistic=20.027918943943487, pvalue=0.028989724934722644)
#sp.stats.friedmanchisquare(synCount[0], synCount[1],  synCount[2],  synCount[3],
#                            synCount[4],  synCount[5], synCount[6], synCount[7],
#                            synCount[8], synCount[9], synCount[10] )
#sp.stats.friedmanchisquare(nonsynCount[0], nonsynCount[1],  nonsynCount[2],  nonsynCount[3],
#                            nonsynCount[4],  nonsynCount[5], nonsynCount[6], nonsynCount[7],
#                            nonsynCount[8], nonsynCount[9], nonsynCount[10] )

#post-hoc analysis (Nemenyi) ONLY if Friedman is significant (<0.05 alpha):
import scikit_posthocs as skph
segSynSig=skph.posthoc_dunn(synCount)
segNonsynSig=skph.posthoc_dunn(nonsynCount)
#segSynSig=skph.posthoc_nemenyi_friedman(synCount.T)
#segNonsynSig=skph.posthoc_nemenyi_friedman(nonsynCount.T)


from matplotlib import colors

Nonsyn_heat_vals=np.asarray(segNonsynSig)
Syn_heat_vals=np.asarray(segSynSig)

fig, (ax, ax2) = plt.subplots(1,2)
cmap = colors.ListedColormap(["white", "gold", "teal", "purple"])
bounds = [-1, 0, 0.049,0.05, 1]
norm = colors.BoundaryNorm(bounds, cmap.N)
mask1 = np.zeros_like(Syn_heat_vals, dtype=np.bool)
mask1[np.triu_indices_from(mask1)] = True
mask1_text=np.ma.masked_array(Syn_heat_vals, mask1)
im = ax.imshow(np.ma.masked_array(Syn_heat_vals, mask1), cmap=cmap, norm=norm)
mask2 = np.zeros_like(Nonsyn_heat_vals, dtype=np.bool)
mask2[np.triu_indices_from(mask2)] = True
im2 = ax2.imshow(np.ma.masked_array(Nonsyn_heat_vals, mask2), cmap=cmap,
                 norm=norm)
for i in range(len(seg_names)):
    for j in range(len(seg_names)):
        if (Syn_heat_vals[i, j] < 0.049) & (Syn_heat_vals[i, j] >0):
            ax.text(j, i, '%.3f' % Syn_heat_vals[i, j], ha="center",
                va="center", color="black", fontsize=7)
        else:
            ax.text(j, i, '%.3f' % Syn_heat_vals[i, j], ha="center",
                va="center", color="w", fontsize=7)
        if (Nonsyn_heat_vals[i, j] < 0.049) & (Nonsyn_heat_vals[i, j] >0):
            ax2.text(j, i, '%.3f' % Nonsyn_heat_vals[i, j], ha="center",
                va="center", color="black", fontsize=7)
        else:
            ax2.text(j, i, '%.3f' % Nonsyn_heat_vals[i, j], ha="center",
                va="center", color="w", fontsize=7)
ax.title.set_text('Synonymous SNPs')
ax.set_xticks(np.arange(len(seg_names)))
ax.set_yticks(np.arange(len(seg_names)))
ax.set_xticklabels(seg_names)
ax.set_yticklabels(seg_names)
ax.set_ylim(len(seg_names)-0.5, -0.5)
ax.xaxis.tick_top()
ax.tick_params(length=0)
ax2.title.set_text('Nonsynonymous SNPs')
ax2.set_xticks(np.arange(len(seg_names)))
ax2.set_yticks(np.arange(len(seg_names)))
ax2.set_xticklabels(seg_names)
ax2.set_yticklabels(seg_names)
ax2.set_ylim(len(seg_names)-0.5, -0.5)
ax2.xaxis.tick_top()
ax2.tick_params(length=0)
fig.colorbar(im2, ticks=[-1, 0, 0.049, 0.05, 1])
#plt.tight_layout()


#Wilcoxon test (nonparametric equivalent of paired t test
# for sig within segments between synonymous and nonsynonymous
for i in range(0, len(synCount)):
    print('seg', i, ':', sp.stats.wilcoxon(synCount[i], nonsynCount[i],
                                           alternative='less',
                                           zero_method='pratt'))

#nonsyn/syn ratio plot of individual points
plt.figure()
plt.scatter(np.log(normGenomeSyn +0.001), np.log(normGenomeNonsyn+0.001),
            alpha=0.3)
x=np.linspace(-7,-2,100)
#x=np.linspace(-3,0,102)
#line of best-fit
np.polyfit(np.log(normGenomeSyn +0.001), np.log(normGenomeNonsyn+0.001), deg=1)
#following depends on above output
y=x
#y=1.1535344*x+1.20277292
plt.plot(x,y,color='red', linestyle='--', label='y=x')
#label='y=1.1535344*x+1.20277292'
plt.ylabel('ln(Mean Genomic Nonsynonymous SNPs + 0.001)')
plt.xlabel('ln(Mean Genomic Synonymous SNPs + 0.001)')
plt.legend(loc='upper left', frameon=False)

ratio=np.log((normGenomeNonsyn+0.001)/np.log(normGenomeSyn +0.001))
#which samples are what
morenonsynSamps=[]
moresynSamps=[]
samesynnonSamps=[]
#more nonsyn than syn samples:
for i in np.where(ratio<1)[0]:
    morenonsynSamps.append(sample_names[i])

#more syn than nonsyn samples:
for i in np.where(ratio>1)[0]:
    moresynSamps.append(sample_names[i])

#more nonsyn than syn samples:
for i in np.where(ratio==1)[0]:
    samesynnonSamps.append(sample_names[i])




#syn/ nonsyn bar plot for nonparametric distribution..bootstrap calculation of
#confidence intervals..flexible - can use mean, median, geometric mean..etc.

#1. Empirical data
#2. Data median (or other measure)
#3. Median of Bootstrap data. Repeat this step n times. The larger the better.
# #4. Calculate differences between bootstrap medians and data median
# synBoots=[]
# nonBoots=[]
# synLower=[]
# synUpper=[]
# nonLower=[]
# nonUpper=[]
# numIterate=10000
# sampSize=int(len(synCount[0]))
#array of medians (cols=orf)
#estimator in pointplot does this step..
# synMed=np.median(synCount, axis=1)
# nonsynMed=np.median(nonsynCount, axis=1)

# #bootstrapping
# for i in range(0, len(seg_names)):
#     synMed[i]
#     nonsynMed[i]
#     for j in range(0, numIterate):
#         synBoots.append(np.random.choice(synCount[i],size=sampSize))
#         nonBoots.append(np.random.choice(nonsynCount[i],size=sampSize))
#     synBootsMed=np.median(synBoots, axis=1)
#     nonBootsMed=np.median(nonBoots, axis=1)
#     synMedDiffs=synBootsMed-synMed[i]
#     nonsynMedDiffs=nonBootsMed-nonsynMed[i]
#     synMedDiffs.sort()
#     nonsynMedDiffs.sort()
#     synLower.append(synMedDiffs[240])
#     synUpper.append(synMedDiffs[9740])
#     nonLower.append(nonsynMedDiffs[240])
#     nonUpper.append(nonsynMedDiffs[9740])

#asymmetric error bars
# asymSynErr=[np.absolute(synLower), np.absolute(synUpper)]
# asymNonErr=[np.absolute(nonLower), np.absolute(nonUpper)]


# plt.figure()
# plt.errorbar(x - 0.25, np.median(synCount, axis=1).T, yerr=asymSynErr, fmt='d', label='Synonymous', color='lightgray', capsize=2, zorder=0)
# plt.errorbar(x + 0.22, np.median(nonsynCount, axis=1).T, yerr=asymNonErr, fmt='d', label='Nonsynonymous', color='dimgray', capsize=2, zorder=1)
# sns.swarmplot(x=long_non_syn_df['variable'], y=long_non_syn_df['value'], hue=long_non_syn_df['mut_type'], dodge=True, alpha=0.6, zorder=50)
# plt.xticks(range(0, seg_num+1), seg_names,fontsize=12.5)
# plt.xlabel('ORF', fontsize=14)
# plt.ylabel('Number SNPs/ Kb', fontsize=14)
# plt.title('Normalized Number of Nonsynonymous and Synonymous SNPs per Open Reading Frame in PRV', fontsize=16)
# plt.legend(title='Mutation Type', frameon=False, bbox_to_anchor=(1.13,0.8), loc='center right', labelspacing=1.12)
# ax=plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

#pi_aa=pd.DataFrame(0, index=range(longest_seg*seg_num),columns=range(len(col_names)))

###########################################################################################################################################
#SHANNON DIVERSITY BUBBLE PLOTS - ORF ONLY#
###########################################################################################################################################

seg_names=('L1','L2','L3','M1','M2','M3','S1','S1','S2','S3','S4')
start=[19, 19, 8, 22, 27, 84, 29, 108, 22, 29, 39]
stop=[3867, 3891,  3868, 2304, 2090, 2342, 1021, 482, 1284,  1093, 986]
sdf_data=[]
filenames=[]
col_names=[]
temp_sdf=[]
shandiv_data=[]
totalGenNt=0

for file in glob.glob('*_shandiv500.csv'):
    filenames.append(file)
    sdf=pd.read_csv(file, sep=',', header=0)
    sdf['Pos']=sdf['Pos'].astype('str').astype('int') #if Pos column not converted to int, treated as str and sorted lexicographically (i.e. 1,10,100,100,2,20,200,etc..)
    sdf['RefGenome']=sdf['RefGenome'].apply(renameSegments)
    sdf=sdf.sort_values(['RefGenome', 'Pos'], ascending=[True, True])
    sdf.reset_index(drop=True, inplace=True)
    sdf_data.append(sdf)
    shandiv_data=pd.concat(sdf_data, axis=1, ignore_index=True)

#seg_names=shandiv_data[0].unique().tolist()

#initialize arrays based
shanSegihDi=[0]*int(len(filenames))
totalShanDiv=[0]*int(len(filenames))
totalShanSnvCount=[0]*int(len(filenames))
shan=[0]*int(len(filenames))

for file in filenames:
    new_file=file.split("-")[0] #change depending on names of files..or comment out if full file names desired
    col_names.append(new_file)

#plot
for xi,xf,seg in zip(start,stop,seg_names):
    temp_sdf=shandiv_data.loc[(shandiv_data[0] == seg) & (shandiv_data[1] > xi-1) & (shandiv_data[1] < xf+1), 2::3]
    temp_sdf.columns=col_names
    temp_sdf=temp_sdf[phylOrder]
    #filter out 'noise'
    temp_sdf[temp_sdf<0.04]=0
    #shan diversities only
    shan=np.vstack((shan,temp_sdf))
    #geometric mean
    shan_array=temp_sdf.to_numpy(dtype='float64')
    mask_shan_array=np.ma.masked_where(shan_array==0, shan_array)
    segShanGeomean=sp.stats.mstats.gmean(mask_shan_array, axis=0)
    #shan div values for all ORFs
    totalShanDiv=np.vstack((totalShanDiv, shan_array))
    #number sites with snv in orf
    shanOrfSnvCount=np.ma.MaskedArray.count(mask_shan_array, axis=0)
    #total number variant sites
    totalShanSnvCount=np.vstack((totalShanSnvCount, shanOrfSnvCount))
    #total number nt in seg
    totalSegNt=len(mask_shan_array)
    #total nt in genome ORFs
    totalGenNt=totalGenNt+totalSegNt
    #ihDi for ORF
    seg_ihDi=(shanOrfSnvCount/totalSegNt)*segShanGeomean
    #array containing values (col=sample, rows=seg)
    shanSegihDi=np.vstack((shanSegihDi, seg_ihDi.data))
    bdf = temp_sdf.unstack().reset_index()
    bdf.columns = list("YXZ")
    #samples with no snvs
    zero_shan_samps=np.asarray(temp_sdf.loc[:, (temp_sdf == 0).all()].columns)
    #names
    indices5 = np.where(np.in1d(phylOrder, zero_shan_samps))[0]
    plt.figure(figsize=(100,250))
    scatterInfo = plt.scatter(x="X", y="Y", data = bdf, s=bdf.Z*sizeScaling, c="Z", vmin=minColor,vmax=maxColor,cmap="inferno", alpha=0.5)
    handles,labels = getBubbleLegend(legendBubbles,scatterInfo,sizeScaling,minColor,maxColor)
    plt.legend(handles,labels,title='Shannon Diversity', frameon=False, bbox_to_anchor=(1.13,0.8), loc='center right', labelspacing=1.12)
    plt.xticks(np.linspace(np.min(bdf['X']),np.max(bdf['X']),5),np.linspace(xi,xf+1,5).astype(int))
    plt.xlim([np.min(bdf['X'])-10,np.max(bdf['X'])+10])
    plt.title('Shannon Diversity per Position of Segment %s ORF for PRV Isolated from Salmon' %seg, fontsize=16)
    plt.ylabel('Salmon Host ID', fontsize=14)
    plt.xlabel('Genomic Position (nt)', fontsize=14, labelpad=10.0)
    #plt.box(False)
    ax=plt.gca()
    ax.set_yticklabels(phylOrder, fontsize=5.5)
    for y in indices5:
        ax.get_yticklabels()[y].set_color("red")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

#drop row 0 that was used for initialization
shanSegihDi=np.delete(shanSegihDi, 0, 0)
totalShanDiv=np.delete(totalShanDiv, 0, 0)
shan=np.delete(shan, 0, 0)

#shan Div DF
shanOnly=pd.DataFrame(shan)

#Geometric mean and Intrahost diversity index for entire genome (ORFs) calc from shann div values
totalShanSnvCount=np.sum(totalShanSnvCount, axis=0)
normShanSnv=totalShanSnvCount/totalGenNt
maskGenShan=np.ma.masked_where(totalShanDiv==0, totalShanDiv)
full_shanGeomean=sp.stats.mstats.gmean(maskGenShan, axis=0)
shanGenIhdi=normShanSnv*full_shanGeomean

#bbox_to_anchor is bounding box whose bottom left corner begins at position x, y. The bounds of a plot are 1, 1.
shanBreakdown=[]
allShan=pd.melt(shanOnly).value
maskShanData=np.ma.masked_where(allShan==0, allShan)
for i in np.linspace(0.0000178, 0.85, 51):
    result=np.count_nonzero((maskShanData > i) & (maskShanData<=(i+0.01)))
    shanBreakdown=np.append(shanBreakdown, result)


#zero counts
2328966-np.count_nonzero(allShan)
#percent of non zero data (sum)
np.count_nonzero(allShan)/2328966*100
#array of percent data breakdown for each snp freq
percShanBreakdown=shanBreakdown/2328966*100

#plot shan distribution
fig, ax1 = plt.subplots()
ax1.hist(maskShanData, bins=50, color= 'purple', alpha=0.8, edgecolor='black')
ax1.set_xlabel('Shannon Diversity')
ax1.set_ylabel('Number of Sites', color='purple')
ax2 = ax1.twinx()
ax2.plot(np.linspace(0.0000178, 0.85, 51), percShanBreakdown, color='orange', linestyle='--')
ax2.set_ylim(0, 0.83)
ax2.set_ylabel('Percent of Sites (%)', color='orange', rotation=-90, labelpad=14)





###########################################################################################################################################
#ALTERNATE ALLELE BUBBLE PLOTS - ORF ONLY#
###########################################################################################################################################

siteRows=[]
segRows=[]
seg_num=10
longest_seg=3934
freq_df=[]
freq_filenames=[]
final=[]
altFreq=[]
colNames=[]
sampName=[]
segFreq=[]
freqBubble=[]
orf_names=('L1','L2','L3','M1','M2','M3','S1','S1','S2','S3','S4')
start=[19, 19, 8, 22, 27, 84, 29, 108, 22, 29, 39]
stop=[3867, 3891,  3868, 2304, 2090, 2342, 1021, 482, 1284,  1093, 986]

#read in %non-ref bases/site data files
freq_filenames = glob.glob('*_lofreq_af')
freq_filenames.sort()
for f in freq_filenames:
    freq_df = pd.read_csv(f, header=None, delim_whitespace=True,
                          names=['product','site','freq'])
    #add line if only sites with >1% variance is desired: freq_df.drop(freq_df[freq_df['freq'] <= 0.01].index, inplace=True)
    temp_freq = afSites(freq_df)
    final.append(temp_freq)
    altFreq=pd.concat(final, axis=1)
    sampName=f.split("-")[0] #change depending on names of files..or comment out if full file names desired
    colNames.append(sampName)

#initialize array with length equal to number of samples
segIhdi=[0]*int(len(freq_filenames))
segSnv=[0]*int(len(freq_filenames))
genFreqs=[0]*int(len(freq_filenames))
pcaFreq=[0]*int(len(freq_filenames))
totalNt=0

#data manipulation
altFreq=altFreq.iloc[:,np.r_[0, 1,  2:int(len(altFreq.columns)):3]]
colNames.insert(0,'site')
colNames.insert(0,'product')
colNames[52]='NGS268'
altFreq.columns=colNames



#"sample_names" variable defined in previous pi diversity section. Sort columns
#predefined order
#segFreq=segFreq[phylOrder]


#plot & analysis
for xi,xf,seg in zip(start,stop,orf_names):
    segFreq=altFreq.loc[(altFreq['product'] == seg) & (altFreq['site'] > xi-1) & (altFreq['site'] < xf+1)]
    segFreq.drop(columns=['product','site'], axis=1, inplace=True)
    segFreq=segFreq[phylOrder]
    #df for PCA
    pcaFreq=np.vstack((pcaFreq,segFreq))
    #geometric mean
    freq_array=segFreq.to_numpy(dtype='float64')
    mask_freq_array=np.ma.masked_where(freq_array==0, freq_array)
    freqGeomean=sp.stats.mstats.gmean(mask_freq_array, axis=0)
    #variant sites (0 sites masked) for each segment appended to genSnv array to be used outside
    #loop to calc geometric mean for full array of ORFs (genome)
    genFreqs=np.vstack((genFreqs, freq_array))
    #number sites with snv in orf
    snvCount=np.ma.MaskedArray.count(mask_freq_array, axis=0)
    #rows=orf, cols=sample name (for summing total genomes later)
    segSnv=np.vstack((segSnv, snvCount.data))
    #total number nt in seg
    segNt=len(mask_freq_array)
    #total number nt in genome
    totalNt=totalNt+segNt
    #ihDi for proportion non-ref bases
    freq_ihDi=(snvCount/segNt)*freqGeomean
    #array containing values (col=sample, rows=seg)
    segIhdi=np.vstack((segIhdi, freq_ihDi.data))
    freqBubble = segFreq.unstack().reset_index()
    freqBubble.columns = list("YXZ")
    freqBubble.Z=freqBubble.Z.astype('float')
    #samples with no snvs
    zeroSamps=np.asarray(segFreq.loc[:, (segFreq == 0).all()].columns)
    #names
    indices3 = np.where(np.in1d(phylOrder, zeroSamps))[0]
    plt.figure(figsize=(100,250))
    scatterInfo=plt.scatter(x=freqBubble['X'], y=freqBubble['Y'],
                                marker='o', s=freqBubble['Z']*sizeScaling,
                                c=freqBubble['Z'], vmin=minColor,vmax=maxColor,
                                cmap="plasma", alpha=0.5)
    #may add to change shape dependeing on mutation type, but need to make sure
    #1. sample order (columns) are same between pi_aa and freq bubble
    #2. snv threshold is set to same level..i.e. if using cutoff of 0.05% in lofreq,
    #change SNPGenie settings to 0.05%, too. then:
    #for mutation, mark in zip(['Nonsynonymous','Synonymous','Ambiguous', 0], ['^','o','x','>']):
    #     markevery=aa_bubble.index[aa_bubble.z==mutation].tolist()
    #     scatterInfo=plt.scatter(x=pi_bubble['X'].iloc[markevery], y=pi_bubble['Y'].iloc[markevery],
    #                             marker=mark, s=pi_bubble['Z'].iloc[markevery]*sizeScaling,
    #                             c=pi_bubble['Z'].iloc[markevery], vmin=minColor,vmax=maxColor,
    #                             cmap="viridis", alpha=0.5)
    # leg2=plt.legend(handles=[triangle, circle], title='SNV Type', frameon=False, bbox_to_anchor=(1.13,0.45), loc='center right', labelspacing=1.12)
    # ax = plt.gca().add_artist(leg2)
    # #
    handles, labels = getBubbleLegend(legendBubbles,scatterInfo,sizeScaling,minColor,maxColor)
    plt.legend(handles, labels, title='SNV Frequency', frameon=False, bbox_to_anchor=(1.13,0.8), loc='center right', labelspacing=1.12)
    plt.xticks(np.linspace(np.min(freqBubble['X']),np.max(freqBubble['X']),5),np.linspace(xi,xf+1,5).astype(int))
    plt.xlim([np.min(freqBubble['X'])-10,np.max(freqBubble['X'])+10])
    plt.title('Proportion of Non-Reference Bases per Variable Position of Segment %s ORF for PRV Isolated from Salmon' %seg, fontsize=16)
    plt.ylabel('Salmon Host ID', fontsize=14)
    plt.xlabel('Genomic Position (nt)', fontsize=14, labelpad=10.0)
    ax=plt.gca()
    ax.set_yticklabels(phylOrder, fontsize=5)
    for y in indices3:
        ax.get_yticklabels()[y].set_color("red")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

#how many genomes with snvs
snvDf=altFreq.iloc[:,2:len(altFreq)]
snvsPerGen=snvDf.astype(bool).sum(axis=0)
avgSnvsPerGen=snvsPerGen/totalNt
percSnvsPerGen=avgSnvsPerGen*100
percSnvsPerGen=percSnvsPerGen[phylOrder]

#snvsPerGen100 was made using above 5 lines of code but with 100rd data

#plot comparison of snvs/genome with 100 rd vs. 500rd
snvsPerGenReordered=snvsPerGen[phylOrder]
snvsPerGen100Reordered=snvsPerGen100[phylOrder]
a, =plt.plot(phylOrder, np.log10(snvsPerGen100Reordered.values), '-o',
                color='gold', alpha=0.4, label='min 100 reads')
b, =plt.plot(phylOrder, np.log10(snvsPerGenReordered.values), '--o',
                color='teal', alpha=0.4, label='min 500 reads')
locs, labels= plt.xticks()
plt.xticks(locs, phylOrder, rotation=90, fontsize=6)
plt.ylabel('log(SNPs per Genome)', fontsize=12)
plt.xlabel('Host ID', fontsize=12, labelpad=10)
plt.legend(title="Read Depth Threshold", handles=[a, b], frameon=False,
           bbox_to_anchor=(1.02, 1.07))
ax=plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#SNP/genome distribution
fig, (ax, ax2) = plt.subplots(2,1)
sns.boxplot(percSnvsPerGen, color='lightblue', ax=ax)
ax.set(xlabel='Percentage Genome Containing SNPs (%)', ylabel='PRV-1 Samples')
ax2.hist(snvsPerGenReordered, bins=50, alpha=0.55, edgecolor='black')
ax2.set(xlabel='SNPs per Genome', ylabel='Number of Samples')


#snv freq only DF
freqOnly=pd.DataFrame(pcaFreq)

genFreqs=np.delete(genFreqs, 0, 0)
segIhdi=np.delete(segIhdi, 0, 0)
pcaFreq=np.delete(pcaFreq, 0, 0)

#vertically sum #snvs per orf in genome then normalize by total length of orfs
#NOTE: if u don't want to include orf p13, delete row=8 (starting at zero row)
#IntraHost Diversity Index of all ORFs..full genome ihDi calculated with alternate allele freqs
genSnv=np.sum(segSnv, axis=0)
normGensnv=genSnv/totalNt
maskGenFreqs=np.ma.masked_where(genFreqs==0, genFreqs)
full_FreqGeomean=sp.stats.mstats.gmean(maskGenFreqs, axis=0)
freqGenIhdi=normGensnv*full_FreqGeomean

#swarmplot of ihDi per sample (y) per segment (x)

viridis=sns.color_palette(palette='viridis', n_colors=seg_num+1, desat=None)

plt.figure()
sns.swarmplot(data=(segIhdi.T)*10000, palette='viridis')
sns.boxplot(data=(segIhdi.T)*10000, whis=np.inf, color='white')
plt.title('Alternate Allele Frequency Intrahost Diversity Index per Segment in PRV' , fontsize=15)
plt.ylabel('Intrahost Diversity Index', fontsize=14)
plt.xlabel('Segment:ORF', fontsize=14, labelpad=10.0)
plt.xticks(range(0, seg_num+1), ['L1:RdRp','L2:capping','L3:core','M1:','M2:','M3:',
                               'S1:sigma4','S1:p13','S2:','S3:','S4:attachment'],
           fontsize=12.5)

#histogram to view intrahost distribution amongst fish per segment
#rows=segments, col=salmon samples
#skip row 0 because this was empty array used to initialize for subsequent vstack
plt.figure()
for seg, color in zip(range(0,len(segIhdi)), ['red', 'blue','green','purple',
                                              'magenta','grey','orange','yellow',
                                              'pink','aqua','black']):
    #plt.figure()
    plt.hist(segIhdi[seg]*10000, bins=4, align='left', orientation='vertical',
              alpha=0.65, edgecolor=color, color=color, linewidth=1.2)
    plt.title('Sample Distribution of Intrahost Diversity for Segment %s' %seg)
    plt.xlabel('Intrahost Diversity Index')
    plt.ylabel('Number of Salmon')

#Friedman test (non-parametric ANOVA) ranking segments to see if there exists
#significant difference between any, but not which is most/ least, etc...for this you
#require post hoc analysis
sp.stats.friedmanchisquare(segIhdi[0], segIhdi[1], segIhdi[2], segIhdi[3],
                           segIhdi[4], segIhdi[5],segIhdi[6],segIhdi[7],
                           segIhdi[8],segIhdi[9],segIhdi[10] )
#Example output for n=5. Note: should be >6
#FriedmanchisquareResult(statistic=14.410909090909115, pvalue=0.10844080924440128)

#1 column freq
freqBreakdown=[]
allFreq=pd.melt(freqOnly).value
maskFreqData=np.ma.masked_where(allFreq==0, allFreq)
for i in np.linspace(0.01, 0.51, 51):
    result=np.count_nonzero((maskFreqData > i) & (maskFreqData<=(i+0.01)))
    freqBreakdown=np.append(freqBreakdown, result)
#zero counts
2328966-np.count_nonzero(allFreq)
#percent of non zero data (sum)
np.count_nonzero(allFreq)/2328966*100
#array of percent data breakdown for each snp freq
percFreqBreakdown=freqBreakdown/2328966*100

#plot snv freq distribution for all 102 samples
fig, ax1 = plt.subplots()
ax1.hist(maskFreqData, bins=50, color= 'mediumblue', alpha=0.8, edgecolor='black')
ax1.set_xlabel('SNP Frequency')
ax1.set_ylabel('Number of Sites', color='mediumblue')
ax2 = ax1.twinx()
ax2.plot(np.linspace(0.01, 0.51, 51), percFreqBreakdown, color='red', linestyle='--')
ax2.set_ylim(0, 0.83)
ax2.set_ylabel('Percent of Sites (%)', color='red', rotation=-90, labelpad=14)
#2.218999999999994% variant sites

###########################################################################################################################################
#PCA (with alt freq)
###########################################################################################################################################
#freqy=altFreq.iloc[:,2:len(altFreq)]
from sklearn import preprocessing
from sklearn.decomposition import PCA

scaledFreqy=preprocessing.scale(pcaFreq)
pca=PCA()
pca.fit(scaledFreqy)
pca_data=pca.transform(scaledFreqy)
per_var=np.round(pca.explained_variance_ratio_*100, decimals=1)
labels=['PC' + str(x) for x in range(1, len(per_var)+1)]
#scree plot
plt.figure()
plt.bar(x=range(1, len(per_var)+1), height=per_var, tick_label=labels)
plt.ylabel('Percentage of Explained Variance (%)')
plt.xlabel('Principle Component')
plt.xticks(range(1, len(per_var)+1), labels, rotation=90, fontsize=6)
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#pca plot
plt.figure()
pca_df=pd.DataFrame(pca_data, index=range(1,22834), columns=labels)
#index=colNames[2::]
plt.scatter(pca_df.PC1,pca_df.PC2)
plt.ylabel('PC2 - {0}%'.format(per_var[1]))
plt.xlabel('PC1 - {0}%'.format(per_var[0]))

#pca annotate
for sample in pca_df.index:
    plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]),
                 fontsize=6)

#top contributors
loading_scores=pd.Series(pca.components_[0], index=range(1,2283))
sorted_loading_scores=loading_scores.abs().sort_values(ascending=False)
top_10_samp_contributors=sorted_loading_scores[0:10].index.values
loading_scores[top_10_samp_contributors]
#very similar means many samples played a role in separation not just one or two


###########################################################################################################################################
#INDIVIDUAL SAMPLE SYN/NON RATIO
###########################################################################################################################################

plt.figure()
#total syn nonsyn counts in genome (normalized by genome concatenated orf length)
plt.scatter(normGenomeNonsyn, normGenomeSyn, alpha=0.4)

GenomeSyn
GenomeNonsyn

#syn/ nonsyn bar plot for nonparametric distribution..bootstrap calculation of
#confidence intervals..flexible - can use mean, median, geometric mean..etc.

#1. Empirical data
#2. Data median (or other measure)
#3. Median of Bootstrap data. Repeat this step n times. The larger the better.
# #4. Calculate differences between bootstrap medians and data median
# synBoots=[]
# nonBoots=[]
# synLower=[]
# synUpper=[]
# nonLower=[]
# nonUpper=[]
# numIterate=10000
# sampSize=int(len(GenomeSyn[0]))
#array of medians (cols=orf)
#estimator in pointplot does this step..
# synMed=np.median(synCount, axis=1)
# nonsynMed=np.median(nonsynCount, axis=1)

# #bootstrapping
# for i in range(0, len(seg_names)):
#     normGenomeSyn[i]
#     normGenomeNonsyn[i]
#     for j in range(0, numIterate):
#         synBoots.append(np.random.choice(GenomeSyn[i],size=sampSize))
#         nonBoots.append(np.random.choice(GenomeNonsyn[i],size=sampSize))
#     synBootsCount=np.sum(synBoots, axis=0)/total_nt
#     nonBootsCount=np.sum(nonBoots, axis=0)/total_nt
#     synCountDiffs=synBootsCount-normGenomeSyn[i]
#     nonsynCountDiffs=nonBootsCount-normGenomeNonsyn[i]
#     synCountDiffs.sort()
#     nonsynCountDiffs.sort()
#     synLower.append(synCountDiffs[240])
#     synUpper.append(synCountDiffs[9740])
#     nonLower.append(nonsynCountDiffs[240])
#     nonUpper.append(nonsynCountDiffs[9740])

#asymmetric error bars
# asymSynErr=[np.absolute(synLower), np.absolute(synUpper)]
# asymNonErr=[np.absolute(nonLower), np.absolute(nonUpper)]


# plt.figure()
# plt.errorbar(x, np.median(synCount, axis=1).T, yerr=asymSynErr, fmt='d', label='Synonymous', color='lightgray', capsize=2)
# plt.errorbar(x, np.median(nonsynCount, axis=1).T, xerr=asymNonErr, fmt='d', label='Nonsynonymous', color='dimgray', capsize=2)


###########################################################################################################################################
#ihDi and nonsyn comparisons by ORF
###########################################################################################################################################

shanSegihDi
seg_ihDi
segIhdi
synCount
nonsynCount

shanIhDi=[]
piIhDi=[]
altIhDi=[]
syn_df=[]
non_df=[]
intraHostIndex_df=[]

#df for swarmplot and pointplot
shanIhDi=pd.DataFrame(shanSegihDi.T)
shanIhDi.columns=list(seg_names)
shanIhDi['index_type']='shanIhDi'
piIhDi=pd.DataFrame(seg_ihDi.T)
piIhDi.columns=list(seg_names)
piIhDi['index_type']='piIhDi'
altIhDi=pd.DataFrame(segIhdi.T)
altIhDi.columns=list(seg_names)
altIhDi['index_type']='altIhDi'
syn_df=pd.DataFrame(synCount.T)
syn_df.columns=list(seg_names)
syn_df['index_type']='Synonymous'
nonsyn_df=pd.DataFrame(nonsynCount.T)
nonsyn_df.columns=list(seg_names)
nonsyn_df['index_type']='Nonsynonymous'
intraHostIndex_df=pd.concat([syn_df, nonsyn_df,shanIhDi, piIhDi, altIhDi])
long_ihDi_df=pd.melt(intraHostIndex_df, id_vars='index_type')

plt.figure()
sns.pointplot(x=long_ihDi_df['variable'], y=long_ihDi_df['value'], hue=long_ihDi_df['index_type'],
              palette={'Synonymous': 'g', 'Nonsynonymous':'blue', 'shanIhDi':'mediumblue', 'piIhDi': 'magenta', 'altIhDi': 'orange'},
              size=4, ci=95, n_boot=1000, estimator=np.median, join=True, markers=['o', '^', 'x', 'D', 's'], scale=0.5, dodge=0.4,
              errwidth=0.8, capsize=0.035)
plt.xticks(range(0, seg_num+1),['Core Shell','Capping Enzyme', 'RdRp','Core NTPase',
                                'Outer Shell', 'NS-Factory','Outer Clamp','Cytotoxic Protein',
                                'Core Clamp', 'NS-RNA Protein', 'Attachment Fiber'], rotation=55, fontsize=12)
plt.xlabel('ORF', fontsize=14)
plt.ylabel('Number SNPs/ Kb', fontsize=14)
plt.title('Intrahost Diversity Indices Compared to Normalized Number of Nonsynonymous and Synonymous Mutations \nper Open Reading Frame in PRV', fontsize=16)
plt.legend(title='Intrahost Diversity Index', frameon=False, bbox_to_anchor=(1,0.8), loc='center right', labelspacing=1.12)
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

###########################################################################################################################################
#interhost snp heatmap PAIRWISE identity heatmap
###########################################################################################################################################

interPairwise=pd.read_csv('Thesis_MUSCLE_aln_matrix.csv')
interPairwise.set_index('Unnamed: 0', inplace=True)

plt.figure()
plt.imshow(interPairwise, cmap='viridis')
plt.tight_layout()
plt.colorbar().set_label('Pairwise Nucleotide Identity (%)', rotation=270,
                         labelpad=10, fontsize=12)
ax=plt.gca()
ax.set_xticks(np.arange(len(interPairwise)))
ax.set_yticks(np.arange(len(interPairwise)))
ax.set_xticklabels(interPairwise.columns, rotation=90, fontsize=5.5)
ax.set_yticklabels(interPairwise.columns, fontsize=5.5)
ax.xaxis.tick_top()
ax.set_xlim(-0.5,len(interPairwise)-0.5)
ax.set_ylim(len(interPairwise)-0.5, -0.5)
ax.set_label('Pairwise Nucleotide Identity (%)')
###########################################################################################################################################
#PCoA (with genome iSNV number to see if cluster by disease state)
#####################################################################################
##MAKES HORSESHOE>>>can't use..use dist matrix in agglom clustering

from sklearn.metrics.pairwise import pairwise_distances
from skbio.stats.ordination import pcoa
genomeSnvs=np.asarray(snvsPerGenReordered).reshape(-1,1)
#or without .reshape(-1,1) can do this: genomeSnvs=np.vstack(genomeSnvs)
pw_bc_snv=pairwise_distances(genomeSnvs, metric='braycurtis')

pcoaSnv=pcoa(pw_bc_snv)
pcoaSNV=pcoaSnv.samples
plt.figure()
plt.scatter(pcoaSNV.PC1, pcoaSNV.PC2, alpha=0.3)

###########################################################################################################################################
#heirarchical clustering (agglomerative)
###########################################################################################################################################

import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering


interDist=np.asarray(interPairwise)
interDist0=np.nan_to_num(interDist)

plt.figure()
dend=sch.dendrogram(sch.linkage(sp.spatial.distance.pdist(interDist0), method='ward',
                                optimal_ordering=True), orientation='left')
#dendro: distance_sort='ascending'
plt.axvline(x=141.36, color='r', linestyle='--')
plt.show()
#ax=plt.gca()
#ax.set_ylim(0,200)

cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
#cluster.fit(sp.spatial.distance.pdist(interDist0).reshape(-1,1))
cluster.fit_predict(sp.spatial.distance.pdist(interDist0).reshape(-1,1))

plt.scatter(interDist0, interDist0)

#retrieve labels
ylab=[]
labelNames=[]

locs,labels=plt.yticks()
for ind in labels:
    ylab.append(int(ind.get_text()))
for label2Ind in ylab:
    labelNames.append(interPairwise.columns[label2Ind])

plt.yticks(locs, labelNames, fontsize=5.5)


######################using iSNPs######################################################
from sklearn import preprocessing
plt.figure()

#scale/ center around 0..
X=preprocessing.scale(genomeSnvs)

dend=sch.dendrogram(sch.linkage(sp.spatial.distance.pdist(X, 'minkowski'), method='complete',
                                optimal_ordering=True), orientation='left')


Z=sch.linkage(sp.spatial.distance.pdist(X, 'minkowski'), method='complete')

#calculate cophenetic coefficient
c, coph_dists = sch.cophenet(Z, sp.spatial.distance.pdist(X, 'minkowski'))
c


#retrieve labels
ylab=[]
labelNames=[]

locs,labels=plt.yticks()
for ind in labels:
    ylab.append(int(ind.get_text()))
for label2Ind in ylab:
    labelNames.append(phylOrder[label2Ind])

plt.yticks(locs, labelNames, fontsize=5.5)

plt.show()



###########################################################################################################################################
#PCoA (with pi sites)
###########################################################################################################################################
from skbio.diversity import beta_diversity

ids=list(range(1, 22834))
bc_dm=sp.spatial.distance.pdist(piGenome, metric='braycurtis')
bc_dm=sp.spatial.distance.squareform(bc_dm)
bc_dm=np.nan_to_num(bc_dm)
pcoa=pcoa(bc_dm)
pcoaSamp=pcoa.samples
plt.scatter(pcoaSamp.PC1, pcoaSamp.PC2, alpha=0.3)

###########################################################################################################################################

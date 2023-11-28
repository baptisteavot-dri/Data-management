import pandas as pd

# import SampleSheet.csv files
keys=list()
for i in range(1,56):
    X=pd.read_csv('keys/key_'+str(i)+'.csv')
    keys.append(X)

# concatenate all files within the array

Z=pd.concat(keys,axis=0)
Z.to_csv('keys.csv')

igfq0067=pd.read_csv('samplesheet_IGFQ001167.csv')
print(igfq0067.columns)
igfq0067['CaseID']=igfq0067['CASE ID']
Z=pd.concat([Z,igfq0067],axis=0)

# remove last part of IGF ID
Z[['Sample_ID','misc','misc_1']]=Z['Sample_ID'].str.split('_',expand=True)
Z.drop(['misc','misc_1'],axis=1,inplace=True)


# remove duplicates
Z.drop_duplicates('Sample_ID',inplace=True)

# merge keys with snRNAseq and bulkRNAseq metadata
def merge(Z,data_type):

    # keeping only relevant columns from key file
    keys=Z[['Sample_ID','Sample_Name']]

    # extracting brain region from sample name
    keys['Brain region'] = keys.Sample_Name.str.extract(pat='(EC|SSC|mTemp|MTEMP|SOM|DPFC|ITG)',expand=False)

    # keeping original sample name
    keys['Sample_Name_original']=keys['Sample_Name']

    # normalizing sample name
    keys['Sample_Name_normalized']=keys['Sample_Name'].str.replace('-','_')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('.','_')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace(' ','_')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('/','_')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_EC','')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_SSC','')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_mTemp','')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_MTEMP','')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_SOM','')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_DPFC','')
    keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_ITG','')

    # merging specific to snRNAseq
    if data_type=='snRNAseq':
        
        # importing snRNAseq metadata
        snRNA=pd.read_csv('snRNAseq_10x.csv')

        # normalizing BrainBankID
        snRNA['Brain Bank ID (normalized)']=snRNA['Brain Bank ID'].str.replace('-','_')
        snRNA['Brain Bank ID (normalized)']=snRNA['Brain Bank ID (normalized)'].str.replace('.','_')
        snRNA['Brain Bank ID (normalized)']=snRNA['Brain Bank ID (normalized)'].str.replace('/','_')

        # setting index for merging
        snRNA.set_index("Brain Bank ID (normalized)",inplace=True)
        keys.set_index('Sample_Name_normalized',inplace=True,drop=False)

        # merging keys with snRNAseq metadata
        snRNA=snRNA.join(keys['Sample_ID'])

        return snRNA
    
    # merging specific to bulkRNAseq
    else:

        # importing AD database
        map_database=pd.read_csv('Map_AD_Database.csv')

        # normalizing caseID that will be used for merging
        map_database['CaseID (normalized)']=map_database['CaseID'].str.replace('-','_')
        map_database['CaseID (normalized)']=map_database['CaseID (normalized)'].str.replace('.','_')
        map_database['CaseID (normalized)']=map_database['CaseID (normalized)'].str.replace('/','_')

        # normalizing BrainBankNetworkID
        map_database['BrainBankNetworkID_normalized']=map_database['BrainBankNetworkID '].str.replace('-','_')
        map_database['BrainBankNetworkID_normalized']=map_database['BrainBankNetworkID_normalized'].str.replace('.','_')
        map_database['BrainBankNetworkID_normalized']=map_database['BrainBankNetworkID_normalized'].str.replace('/','_')

        # setting index for merging
        map_database.set_index('CaseID (normalized)',inplace=True, drop=False)
        keys.set_index('Sample_Name_normalized',inplace=True,drop=False)

        # dropping rows with missing BrainBankID
        map_database.dropna(subset=['BrainBankID'],inplace=True)

        # keeping original BrainBankID
        map_database['BrainBankNetworkID_original']=map_database['BrainBankNetworkID ']
        map_database_1=map_database.set_index('BrainBankNetworkID_normalized',drop=False)

        # merging keys with AD database
        bulkRNA_1=map_database.join(keys[['Sample_ID','Brain region','Sample_Name_normalized','Sample_Name_original']])
        bulkRNA_2=map_database_1.join(keys[['Sample_ID','Brain region','Sample_Name_normalized','Sample_Name_original']])

        # concatenating the two dataframes

        bulkRNA=pd.concat([bulkRNA_1,bulkRNA_2],axis=0)

        return bulkRNA

# apply merge function to snRNAseq and bulkRNAseq
snRNA=merge(Z,'snRNAseq')
bulkRNA=merge(Z,'bulkRNAseq')

# opening text file with all fastq paths
all_fastqs = open("all_fastqs.txt", "r")
all_fastqs_filtered = open("all_fastqs_filtered.csv", "w+")

# filtering fastq paths and creating a csv file
for line in all_fastqs.readlines():

    # filtering fastq paths
    if ('IGF' in line) & (('igf/' in line)|('seqdata/' in line)):
        dics=line.split('/')
        # iterating through all directories in the path
        for i in dics:
            # extracting IGF ID
            if ("IGF" in i) & (len(i)==9):
                line=line[1:]
                # writing in csv format
                all_fastqs_filtered.write(i+ ',' + line)

# closing files
all_fastqs.close()
all_fastqs_filtered.close()

# opening paths csv file and aggregating paths for each IGF ID
paths=pd.read_csv('all_fastqs_filtered.csv',header=None)

# setting base path
paths['base']=' /rds/general/project/ukdrmultiomicsproject/live'
paths[1]=paths['base']+paths[1]

counts=paths[0].value_counts()

# aggregating paths for each IGF ID
paths = paths.groupby(0)[1].agg(';'.join)
paths = paths.to_frame()

# merging paths with snRNAseq metadata
snRNA.set_index('Sample_ID',inplace=True)
snRNA=snRNA.join(paths)

# splitting paths into separate columns
S=snRNA[1].str.split(';',expand=True)
snRNA.drop([1],axis=1,inplace=True)

# joining path dataframe with snRNAseq metadata
snRNA=snRNA.join(S)

# selecting only TREM2 project samples
snRNA=snRNA.loc[snRNA['Study']=='TREM2']
snRNA.dropna(inplace=True,subset=[0])

# dropping unnecessary columns
snRNA.drop(['Sorted?','User','Comment',
                   'FACSSeparation','Kit',
                   'CellRecovery ','cDNAPCRCycles',
                   'cDNABioanalyserDate',
                   'cDNAAveragePeak',
                   'QubitcDNAconcng/ul',
                   'NgtoLibraryPrep',
                   'IndexWell','FinalPCRcycles',
                   'FinalLibBioAnalyserDate',
                   'LibAveragePeak','QbitLibConcNg/Ul',
                   'LibConcnM','Pass/Fail','Comments',
                   'Sequence?','TakenForSeqUl',
                   'SequencingDate','QuoteNumber',
                   'SampleSequencingID','PooledVolume',
                   'Dilutions50/libConcnM',
                   'AllSamplesIds',
                   'StudyID'],axis=1,inplace=True)

# extracting study ID from paths
snRNA['Study_ID']=snRNA[0].str.extract(r'(IGFQ[0-9]{6})')
snRNA.insert(0,'Study_ID',snRNA.pop('Study_ID'))

# keeping original brain region acronyms
bulkRNA['Brain region original']=bulkRNA['Brain region']

# renaming brain regions to standard names
bulkRNA['Brain region']=bulkRNA['Brain region'].str.replace('mTemp','MTG')
bulkRNA['Brain region']=bulkRNA['Brain region'].str.replace('MTEMP','MTG')
bulkRNA['Brain region']=bulkRNA['Brain region'].str.replace('SOM','SSC')

# creating new standardized sample names
bulkRNA['Sample_Name']=bulkRNA['BrainBankNetworkID_normalized']+'_'+bulkRNA['Brain region']+'_S1'

# adding sequencing batch column
bulkRNA['Sequencing_batch']=1

# replacing diagnosis column name
bulkRNA['NeuropathologicalDiagnosis']=bulkRNA['AD/CTRL']

# mapping variant ID to variant names
bulkRNA['TREM2Variant']=bulkRNA['TREM2VariantID'].map({2:'CV',8:'R62H',7:'R47H',1:'N/A',3:'D87N'})
bulkRNA['TREM2Variant']=bulkRNA['TREM2Variant'].fillna('N/A')

# keeping only relevant metadata columns for bulkRNAseq
bulkRNA=bulkRNA[['Sample_ID','Sample_Name','Sample_Name_original',
                 'Sequencing_batch','BrainBankNetworkID_original',
                 'BrainBank','CaseID','Study','Brain region',
                 'Brain region original','Braak','Sex','Age',
                 'NeuropathologicalDiagnosis','TREM2Variant',
                 'APOE','PostMortemDelayHours',
                 'BrainBankNetworkID_normalized']]

# setting the right indeces before merging
bulkRNA.set_index('Sample_ID',inplace=True,drop=False)

# dropping duplicate IGF IDs
bulkRNA.drop_duplicates('Sample_ID',inplace=True)

# joining paths with bulkRNAseq metadata
bulkRNA=bulkRNA.join(paths)

# splitting paths into separate columns
S=bulkRNA[1].str.split(';',expand=True)
bulkRNA.drop([1],axis=1,inplace=True)
bulkRNA=bulkRNA.join(S)

# removing rows that do not contain fastq paths
bulkRNA.dropna(inplace=True,subset=[0])
bulkRNA.drop(['Sample_ID'],axis=1,inplace=True)

# extracting study ID from paths
bulkRNA['Study_ID']=bulkRNA[0].str.extract(r'(IGFQ[0-9]{6})')

# bringing study ID to the first column
bulkRNA.insert(0,'Study_ID',bulkRNA.pop('Study_ID'))

studyid2datatype={'IGFQ000852':'snRNAseq','IGFQ000883':'snRNAseq',
                  'IGFQ000949':'bulkRNAseq','IGFQ000955':'snRNAseq',
                  'IGFQ001110':'snRNAseq','IGFQ001167':'bulkRNAseq',
                  'IGFQ001254':'bulkRNAseq','IGFQ001346':'snRNAseq',
                  'IGFQ001404':'bulkRNAseq','IGFQ001462':'bulkRNAseq',
                  'IGFQ001500':'snRNAseq','IGFQ001509':'bulkRNAseq',
                  'IGFQ001529':'snRNAseq','IGFQ001546':'snRNAseq',
                  'IGFQ001613':'snRNAseq','IGFQ001637':'snRNAseq',
                  'IGFQ001651':'snRNAseq','IGFQ001663':'snRNAseq'}

bulkRNA['data_type']=bulkRNA['Study_ID'].map(studyid2datatype)

pathdatatype={'snRNAseq_TREM2':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/snRNAseq_cortical_tissue/Raw_FASTQ/',
              'bulkRNAseq_TREM2':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/',
              'snRNAseq_MAP':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/MAP/Transcriptomics/snRNAseq_cortical_tissue/Raw_FASTQ/',
              'bulkRNAseq_MAP':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/MAP/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/',
              'bulkRNAseq_Tissue Quality Control':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Associated_tissue_studies/Understanding_post_mortem_effects/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/',
              'snRNAseq_Tissue Quality Control':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Associated_tissue_studies/Understanding_post_mortem_effects/Transcriptomics/snRNAseq_cortical_tissue/Raw_FASTQ/'}


bulkRNA['data_type_project']=bulkRNA['data_type']+'_'+bulkRNA['Study']
bulkRNA['path_datatype']=bulkRNA['data_type_project'].map(pathdatatype)

# creating new names for fastq files
bulkRNA['New_name_R1']=bulkRNA['path_datatype']+ \
    bulkRNA['BrainBankNetworkID_normalized']+ \
    '_'+bulkRNA['Brain region']+ \
    '_S1_L001_R1_001.fastq.gz'
bulkRNA['New_name_R2']=bulkRNA['path_datatype']+ \
    bulkRNA['BrainBankNetworkID_normalized']+ \
    '_'+bulkRNA['Brain region']+ \
    '_S1_L001_R2_001.fastq.gz'

# selecting only one study ID for testing
bulkRNA=bulkRNA.loc[ (bulkRNA['Study_ID']=='IGFQ001509') | (bulkRNA['Study_ID']=='IGFQ001462')]
#bulkRNA=bulkRNA.loc[(bulkRNA['Study_ID']=='IGFQ001167') | (bulkRNA['Study_ID']=='IGFQ000949') | (bulkRNA['Study_ID']=='IGFQ001404')]

# import CD33 genotype data
cd33=pd.read_csv('cd33_genotype.csv')
cd33['individual']=cd33['individual'].str.replace('-','_')
cd33.set_index('individual',inplace=True)

# formatting sample names for merging
bulkRNA['index']=bulkRNA['Sample_Name_original'].str.replace('-','_')
bulkRNA['index']=bulkRNA['index'].str.replace('_EC','')
bulkRNA['index']=bulkRNA['index'].str.replace('_SSC','')
bulkRNA['index']=bulkRNA['index'].str.replace('_mTemp','')

# merging bulkRNAseq metadata with CD33 genotype data
bulkRNA.reset_index(inplace=True,drop=False)
bulkRNA.set_index('index',inplace=True,drop=False)
bulkRNA=bulkRNA.join(cd33)
bulkRNA.set_index('Sample_ID',inplace=True,drop=False)

bulkRNA=bulkRNA.join(counts)

# saving test metadata and paths to csv file
bulkRNA.to_csv('sampleInfo_bulkRNAseq_trem2_internal.csv')
bulkRNA['BrainBankNetworkID']=bulkRNA['BrainBankNetworkID_original']

# saving sample sheet and metadata to csv file
bulkRNA[['Study', 'Study_ID','BrainBankNetworkID',
        'Brain region','Sample_Name']]. \
                to_csv('samplesheet_bulkRNAseq.csv',index=False)

bulkRNA[['BrainBankNetworkID',
                'BrainBank','Braak','Sex','Age',
                'NeuropathologicalDiagnosis',
                'TREM2Variant','APOE','CD33',
                'PostMortemDelayHours','CD33_group']].drop_duplicates('BrainBankNetworkID'). \
                to_csv('metadata_bulkRNAseq.csv', index=False)

bulkRNA[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,'New_name_R1','New_name_R2']]. \
    to_csv('paths_bulkRNAseq.csv',index=False,header=None)

paths=open('paths_bulkRNAseq.csv','r')
R1=open('paths_R1.csv','w+')
R2=open('paths_R2.csv','w+')

with open('paths_bulkRNAseq.csv','r') as paths:
    for line in paths.readlines():
        dic=line.split(',')
        for i in dic:
            if 'R1' in i:
                R1.write(i+',')
            elif 'R2' in i:
                R2.write(i+',')
        R1.write('\n')

R1.close()
R2.close()

bulkRNA['New_name_R1']=bulkRNA['New_name_R1'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')
bulkRNA['New_name_R2']=bulkRNA['New_name_R2'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')

uploaded=pd.read_csv('uploaded.txt',sep='\s+',header=None)
all=set(pd.concat([bulkRNA['New_name_R1'],bulkRNA['New_name_R2']],axis=0).to_list())
uploaded=set(uploaded[1].to_list())
to_upload=all.difference(uploaded)
print(to_upload)
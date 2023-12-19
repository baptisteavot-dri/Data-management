import pandas as pd

# import SampleSheet.csv files
keys=list()
for i in range(1,142):
    X=pd.read_csv('keys/key_'+str(i)+'.csv')
    keys.append(X)

# concatenate all files within the array

Z=pd.concat(keys,axis=0)

igfq0067=pd.read_csv('samplesheet_IGFQ001167.csv')
igfq0067['CaseID']=igfq0067['CASE ID']
Z=pd.concat([Z,igfq0067],axis=0)

# remove last part of IGF ID
Z[['Sample_ID','misc','misc_1']]=Z['Sample_ID'].str.split('_',expand=True)
Z.drop(['misc','misc_1'],axis=1,inplace=True)


# remove duplicates
Z.drop_duplicates('Sample_ID',inplace=True)

# keeping only relevant columns from key file
keys=Z[['Sample_ID','Sample_Name']]
keys['Sample_Name']=keys['Sample_Name'].str.replace('DPFC','FC')
keys['Sample_Name']=keys['Sample_Name'].str.replace('PFC','FC')

# extracting brain region from sample name
keys['Brain region'] = keys.Sample_Name.str.extract(pat='(EC|SSC|mTemp|MTEMP|MTG|SOM|FC|ITG)',expand=False)
keys['Enrichment'] = keys.Sample_Name.str.extract(pat='(UNS|ENR|NEG)',expand=False)
keys['Enrichment'] = keys['Enrichment'].str.replace('NEG','GEN')
keys['Enrichment'] = keys['Enrichment'].str.replace('ENR','GEN')
keys['Brain region'] = keys['Brain region'].str.replace('FC','PFC')

# keeping original sample name
keys['Sample_Name_original']=keys['Sample_Name']

# normalizing sample name

keys['Sample_Name_normalized']=keys['Sample_Name'].str.replace('-','_')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('.','_')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace(' ','_')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('/','_')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_UNS','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_ENR','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_NEG','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_EC','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_SSC','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_mTemp','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_MTEMP','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_MTG','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_SOM','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_DPFC','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_ITG','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('EC','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('SSC','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('mTemp','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('MTEMP','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('SOM','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('DPFC','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('ITG','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('MTG','')
keys['Sample_Name_normalized']=keys['Sample_Name_normalized'].str.replace('_[0-9]$','',regex=True) 

  
# importing AD database
map_database=pd.read_csv('Map_AD_Database.csv')

# normalizing caseID that will be used for merging
map_database['CaseID (normalized)']=map_database['CaseID'].str.replace('-','_')
map_database['CaseID (normalized)']=map_database['CaseID (normalized)'].str.replace('.','_')
map_database['CaseID (normalized)']=map_database['CaseID (normalized)'].str.replace('/','_')
map_database['CaseID (normalized)']=map_database['CaseID (normalized)'].str.replace('*','')

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
map_database_2=map_database.set_index('BrainBankNetworkID_normalized',drop=False)

keys.to_csv('keys.csv')

# merging keys with AD database
keys_merged_1=map_database.join(keys[['Sample_ID','Brain region','Enrichment','Sample_Name_normalized','Sample_Name_original']])
keys_merged_2=map_database_2.join(keys[['Sample_ID','Brain region','Enrichment','Sample_Name_normalized','Sample_Name_original']])

# concatenating the two dataframes

all_data=pd.concat([keys_merged_1,keys_merged_2],axis=0)

# opening text file with all fastq paths
all_fastqs = open("all_fastqs.txt", "r")
all_fastqs_filtered = open("all_fastqs_filtered.csv", "w+")

# filtering fastq paths and creating a csv file
for line in all_fastqs.readlines():

    # filtering fastq paths
    if ('IGF' in line) & (('archived_samples_do_not_use' not in line)):
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
paths['base']=' /'
paths[1]=paths['base']+paths[1]
paths[1]=paths[1].str.replace(' ','')

counts=paths[0].value_counts()

# aggregating paths for each IGF ID
paths = paths.groupby(0)[1].agg(';'.join)
paths = paths.to_frame()

# keeping original brain region acronyms
all_data['Brain region original']=all_data['Brain region']

# renaming brain regions to standard names
all_data['Brain region']=all_data['Brain region'].str.replace('mTemp','MTG')
all_data['Brain region']=all_data['Brain region'].str.replace('MTEMP','MTG')
all_data['Brain region']=all_data['Brain region'].str.replace('SOM','SSC')

# replacing diagnosis column name
all_data['NeuropathologicalDiagnosis']=all_data['AD/CTRL']

# mapping variant ID to variant names
all_data['TREM2Variant']=all_data['TREM2VariantID'].map({2:'CV',8:'R62H',7:'R47H',1:'N/A',3:'D87N'})
all_data['TREM2Variant']=all_data['TREM2Variant'].fillna('N/A')

# keeping only relevant metadata columns for all_dataseq
all_data=all_data[['Sample_ID','Sample_Name_original',
                 'BrainBankNetworkID_original','BrainBank',
                 'Enrichment','CaseID','Study','Brain region',
                 'Brain region original','Braak','Sex','Age',
                 'NeuropathologicalDiagnosis','TREM2Variant',
                 'APOE','PostMortemDelayHours','CaseID (normalized)',
                 'Sample_Name_normalized',
                 'BrainBankNetworkID_normalized']]

# setting the right indeces before merging
all_data.set_index('Sample_ID',inplace=True,drop=False)

# dropping duplicate IGF IDs
all_data.drop_duplicates('Sample_ID',inplace=True)

# joining paths with all_dataseq metadata
all_data=all_data.join(paths)

# splitting paths into separate columns
S=all_data[1].str.split(';',expand=True)
all_data.drop([1],axis=1,inplace=True)
all_data=all_data.join(S)

# removing rows that do not contain fastq paths
all_data.dropna(inplace=True,subset=[0])
all_data.drop(['Sample_ID'],axis=1,inplace=True)

# extracting study ID from paths
all_data[0]=all_data[0].astype(str)
all_data['Study_ID']=all_data[0].str.extract(r'(IGFQ[0-9]{6})')

# bringing study ID to the first column
all_data.insert(0,'Study_ID',all_data.pop('Study_ID'))

# creating dictionary with study ID mapped to data type
studyid2datatype={'IGFQ000852':'snRNAseq','IGFQ000883':'snRNAseq',
                  'IGFQ000949':'bulkRNAseq','IGFQ000955':'snRNAseq',
                  'IGFQ001110':'snRNAseq','IGFQ001167':'bulkRNAseq',
                  'IGFQ001254':'bulkRNAseq','IGFQ001346':'snRNAseq',
                  'IGFQ001404':'bulkRNAseq','IGFQ001462':'bulkRNAseq',
                  'IGFQ001500':'snRNAseq','IGFQ001509':'bulkRNAseq',
                  'IGFQ001529':'snRNAseq','IGFQ001546':'snRNAseq',
                  'IGFQ001613':'snRNAseq','IGFQ001637':'snRNAseq',
                  'IGFQ001651':'snRNAseq','IGFQ001663':'snRNAseq'}

studyid2enrichment={'IGFQ000852':'GEN','IGFQ000883':'UNS',
                    'IGFQ000955':'UNS','IGFQ001110':'UNS',
                    'IGFQ001346':'GEN'}

# adding data type column
all_data['data_type']=all_data['Study_ID'].map(studyid2datatype)

# add enrichment information for specific samples
all_data.loc[ all_data["Enrichment"].isnull(),'Enrichment'] = all_data.loc[all_data["Enrichment"].isnull(),"Study_ID"].map(studyid2enrichment)


# creating dictionary with target directory mapped to data type and project
pathdatatype={'snRNAseq_UNS_TREM2':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/snRNAseq_cortical_tissue/TotalPopulation/Raw_FASTQ/',
              'snRNAseq_GEN_TREM2':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/snRNAseq_cortical_tissue/GliaEnriched/Raw_FASTQ/',
              'bulkRNAseq_TREM2':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/',
              'snRNAseq_MAP':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/MAP/Transcriptomics/snRNAseq_cortical_tissue/Raw_FASTQ/',
              'bulkRNAseq_MAP':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/MAP/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/',
              'bulkRNAseq_Tissue Quality Control':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Associated_tissue_studies/Understanding_post_mortem_effects/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/',
              'snRNAseq_Tissue Quality Control':'/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Associated_tissue_studies/Understanding_post_mortem_effects/Transcriptomics/snRNAseq_cortical_tissue/Raw_FASTQ/'
              }

# creating dictionary with synapse directory ID mapped to data type and project
pathsynapse={'snRNAseq_UNS_TREM2':'syn53165789',
             'snRNAseq_GEN_TREM2':'syn53165787',
              'bulkRNAseq_TREM2':'syn52943643',
              'snRNAseq_MAP':'syn53061538',
              'bulkRNAseq_MAP':'syn53061537',
              'bulkRNAseq_Tissue Quality Control':'syn53061539',
              'snRNAseq_Tissue Quality Control':'syn53061540'
              }

all_data['data_type_project']=all_data['data_type']+'_'+all_data['Study']

# dealing specifically with snRNAseq sample names from the TREM2 study
all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'), 'data_type_project']= \
    all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'data_type']+'_'+ \
    all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'Enrichment']+'_'+ \
    all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'Study']

# mapping target directory with data type and project
all_data['path_datatype']=all_data['data_type_project'].map(pathdatatype)

# mapping synapse directory ID with data type and project
all_data['synapse_dir_id']=all_data['data_type_project'].map(pathsynapse)

# import CD33 genotype data
cd33=pd.read_csv('cd33_genotype.csv')
cd33['individual']=cd33['individual'].str.replace('-','_')
cd33.set_index('individual',inplace=True)

# import RIN data
bulkrna=pd.read_csv('BulkRNA.csv')

# normalizing caseID that will be used for merging
bulkrna['CaseID']=bulkrna['CaseID'].str.replace('-','_')
bulkrna['CaseID']=bulkrna['CaseID'].str.replace('/','_')
bulkrna['CaseID']=bulkrna['CaseID'].str.replace('*','')

# setting index for merging
bulkrna.set_index('CaseID',inplace=True)



# import TissueQC metadata
tissueqc=pd.read_csv('TQC_metadata.csv', encoding='unicode_escape')

"""
tissueqc_key=pd.read_csv('key_tissueqc.csv')
tissueqc_dict=pd.Series(tissueqc_key.Sample_ID.values,index=tissueqc_key.Sample_Name_original).to_dict()
tissueqc.loc[tissueqc['Brain_Bank']=='McGill','Sample_ID']=tissueqc.loc[tissueqc['Brain_Bank']=='McGill','BBN_ID'].map(tissueqc_dict)
print(tissueqc.loc[tissueqc['Brain_Bank']=='McGill','Sample_ID'])
tissueqc.to_csv('TQC_metadata.csv')
"""

# keeping only samples from McGill and Netherlands
tissueqc=tissueqc.loc[(tissueqc['Brain_Bank']=='Netherlands') | (tissueqc['Brain_Bank']=='McGill')]

# normalizing caseID that will be used for merging
tissueqc['Case_ID']=tissueqc['Case_ID'].str.replace('-','_')
tissueqc['Case_ID']=tissueqc['Case_ID'].str.replace('/','_')
tissueqc['Case_ID']=tissueqc['Case_ID'].str.replace('*','')

# setting index for merging
tissueqc.set_index('Sample_ID',inplace=True)

# formatting sample names for merging
all_data['index']=all_data['Sample_Name_original'].str.replace('-','_')
all_data['index']=all_data['index'].str.replace('_EC','')
all_data['index']=all_data['index'].str.replace('_SSC','')
all_data['index']=all_data['index'].str.replace('_mTemp','')

# merging all_dataseq metadata with CD33 genotype data
all_data.reset_index(inplace=True,drop=False)
all_data.set_index('index',inplace=True,drop=False)
all_data=all_data.join(cd33)


# merging all_data with RIN metadata
all_data.set_index('CaseID (normalized)',inplace=True,drop=False)
all_data=all_data.join(bulkrna[['RIN','NanodropRNAConcentration_ng_ul']])

# merging all_data with counts data
all_data.set_index('Sample_ID',inplace=True,drop=False)
all_data=all_data.join(counts)
all_data=all_data.join(tissueqc[['Storage_Temp','Time_Delay']])

# dropping duplicate IGF IDs
all_data.drop_duplicates('Sample_ID',inplace=True)

# creating new standardized sample names
all_data['Sample_Name']= all_data['BrainBankNetworkID_normalized']+'_'+all_data['Brain region']

all_data['Brain region'].fillna('N/A',inplace=True)
all_data['Enrichment'].fillna('N/A',inplace=True)

# dealing specifically with snRNAseq sample names from the TREM2 study
all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'Sample_Name']= \
    all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'BrainBankNetworkID_normalized']+ \
    '_'+all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'Brain region']+ \
    '_'+all_data.loc[(all_data['data_type']=='snRNAseq') & (all_data['Study']=='TREM2'),'Enrichment']

# dealing specifically with samples from IGFQ001254 sequencing project
all_data.loc[all_data['Study_ID']=='IGFQ001254','Sample_Name_original']=all_data.loc[all_data['Study_ID']=='IGFQ001254','Sample_Name_original'].str.replace('-','_')
all_data.loc[all_data['Study_ID']=='IGFQ001254','Sample_Name']=all_data.loc[all_data['Study_ID']=='IGFQ001254','Sample_Name_original']

# adding a column with the number of files for specific sample
all_data['number_of_duplicates']=all_data.groupby(['Sample_Name'])['Sample_Name'].transform('count')

# adding sequencing batch number to sample name
all_data['Sequencing_batch']=all_data.groupby(['Sample_Name']).cumcount()+1
all_data['Sequencing_batch']=all_data['Sequencing_batch'].round(decimals=0).astype(str)
all_data['Sequencing_batch']=all_data['Sequencing_batch'].str.replace('.0','')
all_data['Sequencing_batch']='S'+all_data['Sequencing_batch']
all_data['Sample_Name']=all_data['Sample_Name']+'_'+all_data['Sequencing_batch']

# creating new names for fastq files
all_data['New_name_R1']=all_data['path_datatype']+ \
    all_data['Sample_Name']+'_L001_R1_001.fastq.gz'
all_data['New_name_R2']=all_data['path_datatype']+ \
    all_data['Sample_Name']+'_L001_R2_001.fastq.gz'


# extracting duplicate samples
all_data_duplicates = all_data[all_data.duplicated(['BrainBankNetworkID_original','Brain region','Study_ID'])]
all_data_duplicates=all_data_duplicates.loc[all_data_duplicates['data_type']=='bulkRNAseq']
all_data_duplicates.to_csv('duplicates.csv')

# saving all metadata and paths to csv file
all_data.loc[all_data['Study']=='TREM2'].to_csv('sampleInfo_all_dataseq_trem2_internal.csv', index=False)
all_data.loc[all_data['Study']=='MAP'].to_csv('sampleInfo_all_dataseq_map_internal.csv', index=False)
all_data.loc[all_data['Study']=='Tissue Quality Control'].to_csv('sampleInfo_all_dataseq_tissueQC_internal.csv', index=False)
all_data['BrainBankNetworkID']=all_data['BrainBankNetworkID_original']

# saving sample sheet to csv file
all_data.loc[all_data['Study']=='MAP'][['Study', 'Study_ID','BrainBankNetworkID',
        'Brain region','Sample_Name','data_type']]. \
                to_csv('MAP/samplesheet_MAP.csv',index=False)

# writing MAP metadata to csv file
all_data.loc[all_data['Study']=='MAP'][['BrainBankNetworkID',
                'BrainBank','Braak','Sex','Age',
                'NeuropathologicalDiagnosis',
                'TREM2Variant','APOE',
                'PostMortemDelayHours','RIN']].drop_duplicates('BrainBankNetworkID'). \
                to_csv('MAP/donor_metadata.csv', index=False)

# writing TREM2 sample sheet to csv file
all_data.loc[all_data['Study']=='TREM2'][['Study', 'Study_ID','BrainBankNetworkID',
        'Brain region','Sample_Name','data_type','Enrichment']]. \
                to_csv('TREM2/samplesheet_TREM2.csv',index=False)

# writing TREM2 metadata to csv file
all_data.loc[all_data['Study']=='TREM2'][['BrainBankNetworkID',
                'BrainBank','Braak','Sex','Age',
                'NeuropathologicalDiagnosis',
                'TREM2Variant','APOE','CD33',
                'PostMortemDelayHours','CD33_group','RIN']].drop_duplicates('BrainBankNetworkID'). \
                to_csv('TREM2/donor_metadata.csv', index=False)

# writing TissueQC sample sheet to csv file
all_data.loc[all_data['Study']=='Tissue Quality Control'][['Study', 'Study_ID','BrainBankNetworkID',
        'Brain region','Sample_Name','data_type','Storage_Temp','Time_Delay']]. \
                to_csv('TissueQC/samplesheet_TissueQC.csv',index=False)

# writing TissueQC metadata to csv file
all_data.loc[all_data['Study']=='Tissue Quality Control'][['BrainBankNetworkID',
                        'Brain region','Sample_Name']]. \
                        drop_duplicates('BrainBankNetworkID'). \
                to_csv('TissueQC/donor_metadata.csv', index=False)

all_data.loc[all_data['Study']=='Tissue Quality Control','Sample_Name_original'].to_csv('key_tissueqc.csv')

# select samples from a specific study
all_data=all_data.loc[(all_data['Study_ID']=='IGFQ000883') | (all_data['Study_ID']=='IGFQ000955') | (all_data['Study_ID']=='IGFQ001110') | (all_data['Study_ID']=='IGFQ001346') | (all_data['Study_ID']=='IGFQ001500') | (all_data['Study_ID']=='IGFQ001651') | (all_data['Study_ID']=='IGFQ000852') ]
all_data[['New_name_R1','New_name_R2','synapse_dir_id']].to_csv('to_upload.csv',index=False,header=None)

# function that write the paths of the files to concatenate in a csv file
def write_paths(all_data,data_type):

    # selecting only samples from a specific data type
    all_data=all_data.loc[(all_data['data_type']==data_type)]
                          
    if data_type=='snRNAseq':
        all_data=pd.concat([all_data[range(0,46)],all_data[['New_name_R1','New_name_R2']]],axis=1)
        all_data. \
        to_csv('paths_%s.csv' % (data_type),index=False,header=None)
        
        paths=open('paths_%s.csv' % (data_type),'r')
        R1=open('paths_%s_R1.csv' % (data_type),'w+')
        R2=open('paths_%s_R2.csv' % (data_type),'w+')

        with open('paths_%s.csv' % (data_type),'r') as paths:
            for line in paths.readlines():
                dic=line.split(',')
                for i in dic:
                    if ('R2' in i) & ('synapse_mirror' in i):
                        R2.write(i)
                    elif ('R1' in i) & ('synapse_mirror' in i):
                        R1.write(i)
                    elif 'R2' in i:
                        R2.write(i+',')
                    elif 'R1' in i:
                        R1.write(i+',')
                R1.write('\n')

            R1.close()
            R2.close()
    else:
        all_data=pd.concat([all_data[range(0,25)],all_data[['New_name_R1','New_name_R2']]],axis=1)
        all_data. \
        to_csv('paths_%s.csv' % (data_type),index=False,header=None)

        # opening csv file for both R1 and R2
        paths=open('paths_%s.csv' % (data_type),'r')
        R1=open('paths_%s_R1.csv' % (data_type),'w+')
        R2=open('paths_%s_R2.csv' % (data_type),'w+')

        # writing paths for R1 and R2 in separate files
        with open('paths_%s.csv' % (data_type),'r') as paths:
            for line in paths.readlines():
                dic=line.split(',')
                for i in dic:
                    if ('R2' in i) & ('synapse_mirror' in i):
                        R2.write(i)
                    elif ('R1' in i) & ('synapse_mirror' in i):
                        R1.write(i)
                    elif 'R2' in i:
                        R2.write(i+',')
                    elif 'R1' in i:
                        R1.write(i+',')
                R1.write('\n')

            R1.close()
            R2.close()
    
write_paths(all_data,'bulkRNAseq')
write_paths(all_data,'snRNAseq')

all_data=all_data.loc[(all_data['data_type']=='bulkRNAseq')]
all_data['New_name_R1']=all_data['New_name_R1'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/MAP/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')
all_data['New_name_R2']=all_data['New_name_R2'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/MAP/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')
all_data['New_name_R1']=all_data['New_name_R1'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Associated_tissue_studies/Understanding_post_mortem_effects/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')
all_data['New_name_R2']=all_data['New_name_R2'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Associated_tissue_studies/Understanding_post_mortem_effects/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')
all_data['New_name_R1']=all_data['New_name_R1'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')
all_data['New_name_R2']=all_data['New_name_R2'].str.replace('/rds/general/project/ukdrmultiomicsproject/live/synapse_mirror/Genetically_stratified_cohorts/TREM2/Transcriptomics/BulkRNAseq_cortical_tissue/Raw_FASTQ/','')

uploaded=pd.read_csv('uploaded.txt',sep='\s+')
name2synapseid=pd.Series(uploaded.synapse_file_id.values,index=uploaded.filename).to_dict()
all_data['synapase_file_id_R1']=all_data['New_name_R1'].map(name2synapseid)
all_data['synapase_file_id_R2']=all_data['New_name_R2'].map(name2synapseid)
all_data[['synapase_file_id_R1','synapse_dir_id']].to_csv('synapse_file_id_R1.csv',index=False,header=None)
all_data[['synapase_file_id_R2','synapse_dir_id']].to_csv('synapse_file_id_R2.csv',index=False,header=None)

all=set(pd.concat([all_data['New_name_R1'],all_data['New_name_R2']],axis=0).to_list())
uploaded=set(uploaded.filename.to_list())
missing=all.difference(uploaded)
print(missing)
missing=pd.DataFrame(missing)
missing.to_csv('missing.txt',index=False,header=None,sep='\t')
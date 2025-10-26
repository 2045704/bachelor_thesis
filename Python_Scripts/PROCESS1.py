import R4
from R4 import Fastq
import pandas as pd
import os

# #CREATED THE SALMON INDEX

# wd='/data01/users_space/maria/data/Thesis'
# suffix='.fastq'
# index_file = 'data/genome/Homo_sapiens.GRCh38.cdna.all.fa.gz'

# obj = Fastq(wd,suffix)
# index = obj.Salmon_index(index_file)


srat = pd.read_csv('data/table/GSE171117.SraRunTable.csv')
entries = srat['Run'].tolist()
#entries = ['SRR14102277']

wd='/mnt/alessandro/Volume/Maria/data/Thesis'

suffix='.fastq'
index_file = 'data/genome/Homo_sapiens.GRCh38.cdna.all.fa.gz'

obj = Fastq(wd,suffix)
#index = obj.Salmon_index(index_file)
index = '/mnt/alessandro/Volume/Maria/data/fastq/Salmon_Index'


print(f'THIS IS ALL THE FILES IN THE WD {obj.list_files()}')

#THESE ARE TRUSEQ ADAPTERS
ADAPT1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' 
ADAPT2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' 
cut_opt = "--nextseq-trim=20 --trim-n -a A{100} -A A{100} -a G{100} -A G{100} -a T{100} -A T{100} -a C{100} -A C{100}"



processed =  os.listdir('/mnt/alessandro/Volume/Maria/data/Thesis/Salmon_Out')
unprocessed = [x for x in entries if x not in processed]

print(f'Unprocessed entries: {unprocessed} \n \n \n \n')

for entry in unprocessed:
    obj = Fastq(wd,suffix)
    print(f'Started processing for {entry}')
    try:
        obj.FastqDump(entries=[entry], single=False, paired=True, Skip_Prefetch=True, out_wd=wd)
        
        obj.FastQC()

        obj.Cutadapt(entry, ADAPT1, ADAPT2, cut_opt=cut_opt,zipper=True)
        obj = Fastq(wd,'.fastq.gz')
        obj.tidy_cutadapt(entry)

        obj.FastQC()

        out_dir = os.path.join(wd, "Salmon_Out", entry)
        obj.Run_Salmon(index, libtype='ISR', single=False, paired=True, out_wd=out_dir, options='--quiet')

        for file in obj.list_files():
            os.remove(file) #removing the files to avoid taking up too much space in the server

    except Exception as err:
        print(f"❌ Error processing {entry}: {err}")





# try:
#     obj.MultiQC(folder_title='MQC_OMO')
# except Exception as err:
#     print(f'Something went wrong with {entry}: {err}')







import os
import re
import subprocess
import shutil
import pandas as pd

cellranger_bin='/mnt/alessandro/Volume/Maria/bin/cellranger-9.0.1/cellranger'
sra_bin = '/mnt/alessandro/Volume/Maria/bin/sratoolkit.3.2.1-ubuntu64/bin/'
conda_bin = '/home/maria/miniconda3/bin/'
p_bin='/usr/bin/'

class Fastq():
    def __init__(self,wd,suffix):

        self.wd = wd
        self.suffix=suffix
    
    def list_files(self, name=False, direct=True, quoted=False):
        if name or quoted:
            direct = False
        res = []
        for file in os.listdir(self.wd):
            if file.endswith(self.suffix) and name:
                res.append(file)  # filenames only
            if file.endswith(self.suffix) and direct:
                res.append(self.wd + "/" + file)  # full path
            if file.endswith(self.suffix) and quoted:
                res.append(f'"{self.wd}/{file}"')  # full path with quotes
        return res

    def FastqDump(self,entries=[],Table_wd='',single=False,paired=False,Skip_Prefetch=False,out_wd='',options_fqdump=''):
        
        if single==False and paired==False:
            return print('Please specify single or paired ended files')
            
        if Table_wd!='': #meaning table with accesions is provided 
            SRATable = pd.read_csv(Table_wd)
            entries = SRATable['Run'].tolist()
        
        if entries==[]:
            return print("Please provide a list of entries, what am I supposed to download?")
            
        print(f'Downloading this list of entries: {entries}')

        
        for entry in entries:
            
            if out_wd!='':
                output=' -O '+out_wd
            else:
                output=''
                
            if Skip_Prefetch==False:
                try:
                    prefetch=sra_bin+'prefetch '+entry
                    #prefetch --no-verify <SRA_ID>
    
                    print(prefetch)
                    os.system(prefetch)
                    print(f'\n ...Done prefetching for {entry}\n')
                except Exception as err:
                    print(f"❌ Error prefetching {entry}: {err}")
                    
            if single==True:
                #command='wsl -e '+FatstqDump_bin+' '+entry
                command=sra_bin+'fasterq-dump '+entry+output+options_fqdump
                
            if paired==True:
                
                command=sra_bin+'fasterq-dump '+entry+' --split-files '+output+' '+options_fqdump
            
            print(f'Running the following command: {command}')
            
            os.system(command)
            
            if output!='':
                print(f'Finished downloading for {entry} at this directory: {output}')
            else:
                print(f'Finished downloading for {entry} at this directory: {wd}')
    
    def FastQC(self,FastQC=True, out_wd='', options_fastqc='',options_multiqc=''):
        res=[]
        if out_wd=='':
            out_wd = self.wd + "/FastQC"
        os.makedirs(out_wd, exist_ok=True) #CREATING THE OUT_WD IF IT DOES NOT EXIST

        #fastqc_bin = p_bin+"fastqc"
        fastqc_bin = '/mnt/alessandro/Volume/Maria/bin/FastQC/fastqc'
        
        print(self.wd)
        print(f"Files in directory: {os.listdir(self.wd)}, looking for suffix: {self.suffix}")

        for file in os.listdir(self.wd):
            if file.endswith(self.suffix):
                res.append(f'"{self.wd}/{file}"')  # Full directory with quotes
        in_wd = " ".join(res)#.replace("\\", "\\\\")
        out_wd=f'"{out_wd}"'
        
        print(f'set as input working directory {in_wd}')

        command = fastqc_bin + " --nogroup --quiet -o " + out_wd + options_fastqc + " " + in_wd

        print(f'Running the following command: {command}')
        os.system(command)
    
    def MultiQC(self,folder_title='', in_wd='', options_multiqc=''):
        if in_wd == '':
            in_wd = self.wd  # Fix: use self.wd instead of undefined 'wd'
    
        out_mqc = '/home/maria/MQC_OUTS'
        
        # Run MultiQC
        command = f'multiqc {options_multiqc} -o {out_mqc}/{folder_title} {in_wd}'
        print(command)
        os.system(command)
    


    
    def Cutadapt(self, entry, ADAPT1, ADAPT2='', out_wd='', zipper=True, cut_opt=''):
        if out_wd == '':
            out_wd = self.wd + "/Cutadapt"
        os.makedirs(out_wd, exist_ok=True)
    
        suffix_zip = '.gz' if zipper else ''
    
        res = [x for x in self.list_files() if entry in x]
        names = [x for x in self.list_files(name=True) if entry in x]
        
        # PAIR ENDED MODE
        if ADAPT2 != '':
            print('CUTADAPT --PAIRED ENDED MODE')
            for i in range(0, len(res), 2):
                base_f = names[i].replace(self.suffix, "")
                base_r = names[i+1].replace(self.suffix, "")
    
                out1 = os.path.join(out_wd, base_f + "_forward_trimmed" + self.suffix)
                out2 = os.path.join(out_wd, base_r + "_reverse_trimmed" + self.suffix)
                pair = res[i:i+2]
                in_wd = " ".join(pair)
                
                command = f"cutadapt -a {ADAPT1} -A {ADAPT2} {cut_opt} -o {out1}{suffix_zip} -p {out2}{suffix_zip} {in_wd}"
                print(f'Running the following command: {command}')
                os.system(command)
            print("Done!")
        # SINGLE ENDED MODE
        else:
            print('CUTADAPT --SINGLE ENDED MODE')
            for count, in_wd in enumerate(res):
                base_name = names[count].replace(self.suffix, "")
                out = os.path.join(out_wd, base_name + "_trimmed" + self.suffix)
                command = f"cutadapt -a {ADAPT1} {cut_opt} -o {out}{suffix_zip} {in_wd}"
                print(f'Running the following command: {command}')
                os.system(command)
            print("Done!")


    def tidy_cutadapt(self, entry):
        cutadapt_dir = os.path.join(self.wd, "Cutadapt")
        target_dir = self.wd
        
        # Remove original files that match the entry
        for file in os.listdir(target_dir):
            if file.startswith(entry) and file.endswith(self.suffix):
                file_path = os.path.join(target_dir, file)
                os.remove(file_path)
                print(f"🗑️ Deleted original: {file}")
    
        # Move and rename files back
        for file in os.listdir(cutadapt_dir):
            if file.endswith(self.suffix):
                src_path = os.path.join(cutadapt_dir, file)
                cleaned_name = file.replace('_trimmed', '').replace('_forward_trimmed', '').replace('_reverse_trimmed', '')
                dst_path = os.path.join(target_dir, cleaned_name)
                shutil.move(src_path, dst_path)
                print(f"✅ Moved: {file} → {cleaned_name}")


    
    def bamtofastq(bam_file, out_wd, options=None):
        if options is not None:
            command = cellranger_bin + 'bamtofastq ' + options + ' ' + bam_file + ' ' + out_wd
        else:
            command = cellranger_bin + ' bamtofastq ' + bam_file + ' ' + out_wd
        print(f'Running the following command: {command}')
        os.system(command)


    def cellranger_count(ID, create_bam='false', options=None):
        if options is None:
            command = cellranger_bin + ' count --id ' + ID + ' --create-bam ' + create_bam
        else:
            command = cellranger_bin + ' count --id ' + ID + ' --create-bam ' + create_bam + options
            
        print(f'Running the following command: {command}')
        os.system(command)
      
    
    
    def Salmon_index(self,index_file,index_opt=''):
        if index_opt!='':
            index_opt=' '+index_opt
        index=self.wd+'/Salmon_Index'
        command = p_bin+'salmon index -t '+index_file+' -i '+index+index_opt
        print(f'Running the following command: {command}')
        os.system(command)
        return index
        #example command = ./bin/salmon index -t transcripts.fa -i transcripts_index
    
    def Run_Salmon(self, index, libtype='', single=False, paired=False, out_wd='', options=''):
        if out_wd == '':
            out_wd = self.wd + '/Salmon_Out'
        os.makedirs(out_wd, exist_ok=True)
    
        if options != '':
            options = ' ' + options
    
        if not (single or paired):
            print('Please specify single or paired ended files')
            return
    
        if libtype == '':
            print('Please specify libtype')
            return
    
        res = self.list_files()
        print(f"{res}")
    
        if single:
            res = ' '.join(res)
            command = p_bin + f'salmon quant -i {index} -l {libtype} -r {res} --validateMappings -o {out_wd}{options}'
            print(f'Running: {command}')
            os.system(command)
    
        if paired:
            # improved logic: handle _1 and _2 fastq files
            forward = sorted([i for i in res if '_1' in i and i.endswith(self.suffix)])
            reverse = sorted([i for i in res if '_2' in i and i.endswith(self.suffix)])
    
            if not forward or not reverse:
                print("No forward/reverse files detected.")
                return
    
            forward = ','.join(forward)
            reverse = ','.join(reverse)
            print(f"Forward Reads: {forward}\nReverse Reads: {reverse}")
    
            command = p_bin + f'salmon quant -i {index} -l {libtype} -1 {forward} -2 {reverse} --validateMappings -o {out_wd}{options}'
            print(f'Running: {command}')
            os.system(command)

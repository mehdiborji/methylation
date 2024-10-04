import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
parts_batch = args.parts_batch

######################################################

print("current dir:", os.getcwd())
os.chdir(indir)
print("current dir changed to:", os.getcwd())
# ref_dir = "/n/scratch/users/m/meb521/GRCh38_v44/"

multiplier_per_task = 5

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

for p in range(start, end):
    
    p = str(p).zfill(3)

    bam_suffix = 'clean_trim_bismark_bt2_pe.bam'
    bam = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_{bam_suffix}"
    
    report_suffix = 'clean_trim_bismark_bt2_pe.non-conversion_filtering.txt'
    report = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_{report_suffix}"
    
    command = "/home/meb521/Bismark-0.24.2/filter_non_conversion"
    options = "--paired --threshold 5 --consecutive"
    
    submit = f"{command} {options} {bam}"
    
    if os.path.exists(report):
        file_size = os.path.getsize(report)
    
        if file_size == 0:
            # Run your command if the file size is zero
            print(f"The file '{report}' is empty. Running the command...")
            print(submit)
            subprocess.run(submit, shell=True)
        else:
            print(f"The file '{report}' is not empty. Size: {file_size} bytes")
    else:
        print(f"The file '{report}' does not exist.")
        print(submit)
        subprocess.run(submit, shell=True)
    
    
    
    

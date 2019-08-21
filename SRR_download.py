import sys
import os

#samples.txt output dir
def dump_paired(SRR, args):
	os.system("fastq-dump -A " + SRR + " -O " + args[1] + " --split-files -gzip" )

def dump_single(SRR, args):
	os.system("fastq-dump -A " + SRR + " -O " + args[1] + " -gzip") 

def main():
	args = sys.argv[1:]
	if len(args) < 1:
		print("----------------INSTRUCTIONS-----------------")
		print("Before trying to download data, create a samples.txt file with header:")
		print("run type(paired_or_single) experimental details \n")
		print("1. Save Your samples.txt file in the same directory as SRR_download.py")
		print("2. In Cygwin type: python3 SRR_download.py samples.txt output_directory")
	else:
		des = list(open(args[0]))
		print("Downloading Data...")
		for i in range(1,len(des)):
			info = des[i].split("\t")
			SRR = info[0]
			print(SRR)
			if info[1] == "paired":
				dump_paired(SRR, args)
			else:
				dump_single(SRR, args) 

if __name__ == "__main__":
	main()
import os
import sys
cwd = os.getcwd()

def quant_paired(SRR,samples, species_index, output_dir):
	output_dir = SRR+"/"
	pairs = samples[SRR]
	os.system("kallisto quant -i " +species_index+" -o "+output_dir+" "+pairs[0]+" "+pairs[1])

def quant_single(SRR,samples, species_index, output_dir):
	output_dir = SRR+"/"
	reads = samples[SRR]
	os.system("kallisto quant -i "+species_index+" -o "+output_dir + " --single -l 180 -s 20" +" "+reads)


def main():
	#samples.txt, ExperimentID, index
	args = sys.argv[1:] 
	if len(args) < 1:
		print("----------------INSTRUCTIONS-----------------")
		print("Before trying to map data:")
		print("1. Save Your samples.txt file in the same directory as Fastq_Kallisto.py")
		print("2. Save Your Fastq files and indeces of interest in the same directory as Fastq_Kallisto.py")
		print("3. In Cygwin type: python3 Fastq_Kallisto.py samples.txt ExperimentID Index (e.g. Mus_musculus.idx)")
	else:
		des = list(open(args[0]))
		SRR_to_type = {}
		
		#SRR to single or paired for downstream
		for i in range(1,len(des)):
			SRR = des[i].split("\t")[0]
			p_or_s = des[i].split("\t")[1]
			SRR_to_type[SRR] = p_or_s

		Experiment_ID = args[1]
		species_index = args[2]
		samples = {}
		for file in os.listdir(cwd):
			if "fastq" in file:
				if "_" in file: #then data is paired end
					sample = file.split("_")[0]
					if sample not in samples.keys():
						samples[sample] = [file]
					else:
						samples[sample].append(file) 
				else: #then data is single end
					sample = file.split(".")[0]
					samples[sample] = file    
	

		#Runnining Kallisto quant.
		for SRR in samples.keys():
			if SRR_to_type[SRR].lower() == "paired":
				quant_paired(SRR,samples, species_index, Experiment_ID)
			elif SRR_to_type[SRR].lower() == "single":
				quant_single(SRR,samples, species_index, Experiment_ID)
			os.system("mv " + "testread.sam " + SRR)

		#Organize results for downstream DEseq
		os.system("mkdir "+Experiment_ID+"_results")
		for SRR in samples:
			os.system("mv "+SRR+" "+ Experiment_ID+"_results")


if __name__ == "__main__":
	main()


import statistics
import os
import sys
cwd = os.getcwd()

def main():
	args = sys.argv[1:]
	if len(args) < 1:
		print("Inputs: genome (mm10 or hg38)")
		print("This tool takes ChIP-seq peak.bed files and uses HOMER2 to perform annotation.")
		print("After annotation, the fraction of peaks that overalp with repeats are counted")
		print("and a Z-score transformation is used to rank to repeats from most to least likely bound")
	else:
		genome = args[0]
		files = [f for f in os.listdir(cwd) if "_mESC" in f]
		print(files)
		for f in files:
			outfile = f + "_annotated.bed"
			os.system("annotatePeaks.pl " + f + " " + genome + " > " + outfile)
			peaks = list(open(outfile))
			#count overlaps
			overlaps = {}
			for i in range(1,len(peaks)):
				curr = peaks[i].split("\t")
				anno = curr[8]
				if "|" in anno:
					sub_fam = anno
					if sub_fam not in overlaps.keys():
						overlaps[sub_fam] = 1
					else:
						overlaps[sub_fam] += 1

			#get mean and standard_dev of overlaps
			fracs = []
			for k in overlaps.keys():
				overlaps[k] = overlaps[k] / len(peaks)
				fracs.append(float(overlaps[k]))
			avg = statistics.mean(fracs)
			std_dev = statistics.stdev(fracs)

			# #Transform to Z-score and write to file
			results = open(outfile.split("annotated")[0]+"_ranked_repeats.txt", "w")#edit this line
			results.write("Repeat" + "\t" + "Zscore" + "\n")
			for k in overlaps.keys():
				overlaps[k] = (overlaps[k] - avg) / std_dev
				results.write(k + "\t" + str(overlaps[k]) + "\n")

			print("Complete with no issues!")

if __name__ == "__main__":
	main()

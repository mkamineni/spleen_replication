import pandas as pd
import glob
import os 

# Load GWAS file

path = "/data/cteu/2_users/zyu_lab/mkaminen/data/gwas_results/"
#path = "/data/cteu/2_users/zyu_lab/mkaminen/data/old_gwas_results/"
#path = "/data/cteu/2_users/zyu_lab/mkaminen/data/tr_nolog_gwas_results/"

for gwas_stats in glob.glob(os.path.join(path, "*metal1.txt")):
	filename = os.path.basename(gwas_stats)
	#filename = filename.replace(".metal1.txt", "_tr")
	df = pd.read_csv(gwas_stats, sep="\t")

	# Split MarkerName if formatted like '1:123456:A:G'
	split = df['MarkerName'].str.split(':', expand=True)
	df['chr'] = split[0]
	df['pos'] = split[1].astype(int)

	# Create BED columns (0-based start)
	df['start'] = df['pos'] - 1
	df['end'] = df['pos']

	# Keep relevant columns
	bed = df[['chr', 'start', 'end', 'MarkerName', 'Effect', 'StdErr', 'Direction']]

	# Write BED file
	bed.to_csv(filename+".bed", sep='\t', header=False, index=False)


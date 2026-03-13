import pandas as pd
import requests
import numpy as np

p_val = 0.05

datapath = "../data/"

sigsnps_file = "sigsnps_ukbiobank_spleen_wo_heme_wpval.txt"

def read_bed_hg19(path):
    # BED columns: chr, start, end, SNP_ID
    df = pd.read_csv(path, sep="\t", header=None, names=["chr","start","end","SNP_ID"])
    df["pos"] = df["start"] + 1   # Convert 0-based BED to 1-based position
    return df

def split_marker(df):
    split_cols = df["MarkerName"].str.split(":", expand=True)
    df["ref"] = split_cols[2]
    df["alt"] = split_cols[3]  
    return df

def allele_match(row):
    return (row.ref_ext == row.ref_sig and row.alt_ext == row.alt_sig)

def allele_match_flipped(row):
    return(row.ref_ext == row.alt_sig and row.alt_ext == row.ref_sig)

#with open(datapath+sigsnps, 'r') as file:
	#data_string = file.read()
	#sigsnps = eval(data_string)
	#print(sigsnps)

sigsnps = pd.read_csv(datapath+sigsnps_file)
#sigsnps = sigsnps.groupby('Feature')['MarkerName'].apply(list).to_dict()
sigsnps = (
    sigsnps.groupby('Feature')
           .apply(lambda df: list(zip(df['MarkerName'], df['BETA'])))
           .to_dict()
)

rep_sig_snps_output = pd.DataFrame(columns=["Feature", "DiscoverySNP", "MarkerName", "Pvalue"])

for feature,snps in sigsnps.items():
	print(feature)
	rep_gwas_path = datapath+"gwas_results/"+feature+"_metal1.txt"
	rep_gwas = pd.read_csv(rep_gwas_path, sep = "\t")

	# Read the hg19 BED file from liftover
	hg19_bed_path = feature+"_metal1.txt_hg19.bed"
	bed = pd.read_csv(hg19_bed_path, index_col=False, sep="\t", header=None, names=["chr","pos","end","MarkerName","other1","other2"])
	
	# Keep only necessary columns for merging
	bed = bed[["chr","end","MarkerName"]]
	bed["pos"] = bed["end"].astype(int)

	# Merge BED coordinates with original replication GWAS on MarkerName
	rep_gwas_hg19 = pd.merge(rep_gwas, bed, on="MarkerName")

	rep_gwas_sig = rep_gwas_hg19#[rep_gwas_hg19["P-value"] < p_val]

	rep_gwas_sig_split = split_marker(rep_gwas_sig)
	df_sig = pd.DataFrame([s.split(":") + [s, beta] for s,beta in snps], columns=["chr","pos","ref","alt","DiscoverySNP", "BETA"])
	df_sig["pos"] = df_sig["pos"].astype(int)
	df_sig["Direction_disc"] = df_sig.apply(lambda row: '+' if row['BETA'] >= 0 else '-', axis=1)
		
	# Merge on chr + pos
	merged = pd.merge(rep_gwas_sig_split, df_sig, on=["chr","pos"], suffixes=("_ext","_sig"))

	merged["allele_match"] = merged.apply(allele_match, axis=1)
	merged["allele_match_flipped"] = merged.apply(allele_match_flipped, axis=1)

	replicated = merged[merged["allele_match"] | merged["allele_match_flipped"]]
	replicated["Direction_rep"] = np.where(replicated["Effect"] >= 0, "+", "-")

	print(replicated.Direction_disc.value_counts())
	print(replicated.Direction.value_counts())
	print(replicated.allele_match.value_counts())
	print(replicated.allele_match_flipped.value_counts())
	#replicated_dir_agree = replicated[(replicated["allele_match"] & replicated["Direction_disc"]==replicated["Direction"]) | (replicated["allele_match_flipped"] & replicated["Direction_disc"]!=replicated["Direction"])]
	# Flip replication direction if alleles are flipped
	replicated["Direction_aligned"] = replicated.apply(
    		lambda row: row["Direction_rep"] if row["allele_match"] else ('+' if row["Direction_rep"]=='-' else '-'),
    	axis=1
	)


	# Now compare discovery vs aligned replication
	replicated_dir_agree = replicated[replicated["Direction_disc"] == replicated["Direction_aligned"]]
	
	#replicated_dir_agree = replicated[
   	#	((replicated["allele_match"]) & (replicated["Direction_disc"] == replicated["Direction"])) |
    	#	((replicated["allele_match_flipped"]) & (replicated["Direction_disc"] != replicated["Direction"]))
	#]

	if not replicated.empty:
		rep_sig_snps_output = pd.concat([
			replicated_dir_agree,
			replicated_dir_agree[["DiscoverySNP","MarkerName","P-value"]].assign(Feature=feature)
		], ignore_index=True)

	if len(snps)>0:
		print("Overalll Sig SNPs in discovery "+str(len(snps)))
		print("Overlap between cohorts: "+str(len(replicated)))
		print("SNPs with directional agreement: "+str(len(replicated_dir_agree)))
		print("Fraction of SNPs with directional agreement: "+str(len(replicated_dir_agree)/len(replicated)))

rep_sig_snps_output.to_csv(datapath+"rep_sigsnps_directional_final.csv", index=False)

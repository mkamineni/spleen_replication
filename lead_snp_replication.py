import pandas as pd
import requests
from io import StringIO

p_val = 0.05

datapath = "../data/"

sigsnps_file = "sigsnps_ukbiobank.txt"

LDLINK_TOKEN = "1b005bd89aff"

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
    return (row.ref_ext == row.ref_sig and row.alt_ext == row.alt_sig) or \
           (row.ref_ext == row.alt_sig and row.alt_ext == row.ref_sig)

with open(datapath+sigsnps_file, 'r') as file:
	data_string = file.read()
	sigsnps = eval(data_string)
	#print(sigsnps)

rep_sig_snps_output = pd.DataFrame(columns=["Feature", "DiscoverySNP", "MarkerName", "Pvalue"])

for feature,snps in sigsnps.items():
	if "variance" in feature:
		continue

	print(feature)
	#print(snps)
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

	rep_gwas_sig = rep_gwas_hg19[rep_gwas_hg19["P-value"] < p_val]
	#rep_gwas_sig = rep_gwas_hg19

	rep_gwas_sig_split = split_marker(rep_gwas_sig)
	#print(rep_gwas_sig_split.head(5))
	#df_sig = pd.DataFrame([s.split(":") + [s] for s in snps], columns=["chr","pos","ref","alt","DiscoverySNP"])
	df_sig = pd.read_csv(datapath+"fuma_results/"+feature+"/leadSNPs.txt", sep="\t")
	#print(df_sig.head(5))
	new_columns = df_sig['uniqID'].str.split(':', expand=True)
	new_columns.columns = ["chr", "pos", "ref", "alt"]
	df_sig = pd.concat([df_sig[['uniqID']], new_columns], axis=1)
	df_sig['chr'] = df_sig['chr'].apply(lambda x: 'chr' + str(x))
	df_sig.rename(columns = {'uniqID':'DiscoverySNP'}, inplace=True)

	df_sig["pos"] = df_sig["pos"].astype(int)

	#print(df_sig.head(5))
	# Merge on chr + pos
	merged = pd.merge(rep_gwas_sig_split, df_sig, on=["chr","pos"], suffixes=("_ext","_sig"))
	#print(merged.head(5))

	merged["allele_match"] = merged.apply(allele_match, axis=1)

	replicated = merged[merged["allele_match"]]

	if not replicated.empty:
		rep_sig_snps_output = pd.concat([
			rep_sig_snps_output,
			replicated[["DiscoverySNP","MarkerName","P-value"]].assign(Feature=feature)
		], ignore_index=True)
	
	if len(df_sig)>0:
                print(len(replicated))
                print(len(df_sig))
                print(len(replicated)/len(df_sig))

	'''
	# Find SNPs in df_sig not in merged (i.e., not in replication cohort)
	missing_snps = df_sig.loc[~df_sig['DiscoverySNP'].isin(merged['DiscoverySNP']), 'DiscoverySNP'].tolist()
	missing_snps = [":".join(string.split(":", 2)[:2]) for string in missing_snps]

	
	if missing_snps:
    		# Batch request LDproxy for all missing SNPs at once
		r2_threshold = 0.8
		pop = "EUR"
		batch_results = []

	for snp in missing_snps:
		url = "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var={snp}&pop={pop}&r2_d=r2&token=1b005bd89aff"
		resp = requests.get(url)
		if resp.status_code == 200:
			resp = requests.get(url)
			print(resp.status_code)
			print(resp.text[:500])  # first 500 characters to see what’s returned
			df_ld = pd.read_csv(StringIO(resp.text), sep="\t")
		
		# Filter SNPs with r2 >= threshold
		print(df_ld.columns)  # see what column contains R² values
		r2_col = 'R2' if 'R2' in df_ld.columns else 'r2'
		df_ld = df_ld[df_ld[r2_col] >= r2_threshold]
		#df_ld = df_ld[df_ld['R2'] >= r2_threshold]
		if not df_ld.empty:
			df_ld['DiscoverySNP'] = snp
			batch_results.append(df_ld)

	if batch_results:
		ld_snps_df = pd.concat(batch_results, ignore_index=True)	
		# Merge these LD SNPs with rep_gwas_sig_split to see if any are in replication cohort
		ld_merged = pd.merge(rep_gwas_sig_split, ld_snps_df, left_on='MarkerName', right_on='RS_Number', how='inner')
		ld_merged["allele_match"] = ld_merged.apply(allele_match, axis=1)
		replicated_ld = ld_merged[ld_merged["allele_match"]]

	if not replicated_ld.empty:
		rep_sig_snps_output = pd.concat([rep_sig_snps_output, replicated_ld[["DiscoverySNP", "MarkerName", "P-value"]].assign(Feature=feature)], ignore_index=True)

	if len(df_sig) > 0:
		direct_count = len(replicated) if 'replicated' in locals() else 0
		proxy_count = len(replicated_ld) if ('replicated_ld' in locals() and not replicated_ld.empty) else 0
		total_replicated = direct_count + proxy_count

		print("Replicated directly: "+str(direct_count))        
		print("Replicated via proxy/LD: "+str(proxy_count))
		print("Total replication rate = {}/{} ({:.2f})".format(total_replicated, len(df_sig), total_replicated/len(df_sig)))
	else:
		print("No discovery SNPs found for this feature, skipping replication summary.")
	'''

rep_sig_snps_output.to_csv(datapath+"rep_sigsnps_lead_final.csv", index=False)

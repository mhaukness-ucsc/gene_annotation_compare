# gff3_compare.py
# compare two different gff3 gene annotation files on the same reference 

from argparse import ArgumentParser
import csv

def parse_attrs(attrs_str, info_field_types):
	""" Parse in a string of attrs from a gff3 file. """
	attrs = {}
	attrs_list = attrs_str.split(';')
	for attr_pair in attrs_list:
		field_name = attr_pair.split('=')[0]
		field_val = attr_pair.split('=')[1]
		attrs[field_name] = field_val
		if field_name not in info_field_types:
			info_field_types.append(field_name)
	return attrs

def read_gff3_file(gff3_file_name):
	""" Read in a gff3 file. """

	all_feature_types = ['region', 'ncRNA_gene', 'lnc_RNA', 'exon', 'miRNA', 'pseudogene', 'pseudogenic_transcript', 'gene', 'mRNA', 'three_prime_UTR', 'CDS', 'five_prime_UTR', 'snRNA', 'ncRNA', 'scRNA', 'snoRNA', 'unconfirmed_transcript', 'V_gene_segment', 'C_gene_segment', 'J_gene_segment', 'rRNA', 'D_gene_segment', 'start_codon', 'stop_codon', 'intron', 'transcript']
	gene_types = ['ncRNA_gene', 'pseudogene', 'gene']
	# ebi tx types
	# tx_types = ['lnc_RNA', 'miRNA', 'pseudogenic_transcript', 'mRNA', 'snRNA', 'ncRNA', 'scRNA', 'snoRNA', 'unconfirmed_transcript', 'V_gene_segment', 'C_gene_segment', 'J_gene_segment', 'rRNA', 'D_gene_segment']
	tx_types = ['transcript', 'lnc_RNA', 'miRNA', 'pseudogenic_transcript', 'mRNA', 'snRNA', 'ncRNA', 'scRNA', 'snoRNA', 'unconfirmed_transcript', 'V_gene_segment', 'C_gene_segment', 'J_gene_segment', 'rRNA', 'D_gene_segment']
	# Other field types: region, exon, three_prime_UTR, five_prime_UTR, CDS
	all_info_field_types = []

	regions = []

	gene_ids = []
	parent_gene_ids = []
	tx_ids = []
	parent_tx_ids = []

	gene_locs = {}
	tx_locs = {}

	with open(gff3_file_name) as in_gff3:
		r = csv.reader(in_gff3, delimiter='\t')
		for row in r: 
			if not row[0].startswith('#'):
				feature_type = row[2]
				if feature_type not in all_feature_types:
					all_feature_types.append(feature_type)
					print("Warning: feature type: ", feature_type, "is not present in the code.")
				if feature_type == "region":
					regions.append((row[0], int(row[3]), int(row[4])))
				elif feature_type in gene_types:
					attrs_str = row[8]
					gene_attrs = parse_attrs(attrs_str, all_info_field_types)
					# print(gene_attrs)
					try:
						gene_ids.append(gene_attrs['gene_id'])
					except: 
						print("Error, no key 'gene_id'. Gene attrs: ", gene_attrs)
					try:
						parent_gene_id = gene_attrs['parent_gene']
					except: 
						parent_gene_id = gene_attrs['source_gene']
						# print("Error, no key 'parent_gene'. Gene attrs: ", gene_attrs)
					parent_gene_ids.append(parent_gene_id)
					gene_locs[parent_gene_id] = (row[0], int(row[3]), int(row[4]))
				elif feature_type in tx_types:
					attrs_str = row[8]
					tx_attrs = parse_attrs(attrs_str, all_info_field_types)
					# print(tx_attrs)
					try:
						tx_ids.append(tx_attrs['transcript_id'])
					except: 
						print("Error, no key 'transcript_id'. Tx attrs: ", tx_attrs)
					try:
						parent_tx_id = tx_attrs['parent_transcript']
					except: 
						parent_tx_id = tx_attrs['source_transcript']
						# print("Error, no key 'parent_transcript'. Tx attrs: ", tx_attrs)
					parent_tx_ids.append(parent_tx_id)
					tx_locs[parent_tx_id] = (row[0], int(row[3]), int(row[4]))
				else:
					# print(row)
					# break
					continue

	# print("regions: ", regions)
	# print("all info field types: ", all_info_field_types)
	print("\t# genes: ", len(gene_ids))
	print("\t# transcripts: ", len(tx_ids))
	# print(all_feature_types)
	return gene_locs, tx_locs

def compare_feature_locs(ref_locs, query_locs, ref_out_file_name, query_out_file_name, diff_locs_file_name):
	""" Compare the locations for two sets of features """
	num_overlapping_ids = 0
	num_ref_ids_only = 0
	num_query_ids_only = 0

	num_different_chroms = 0
	num_same_chrom_diff_loc = 0
	num_same_loc = 0 

	with open(ref_out_file_name, 'w') as ref_out, open(diff_locs_file_name, 'w') as diffs_out:
		for feature_id, loc in ref_locs.items():
			ref_chrom, ref_start, ref_stop = loc[0], loc[1], loc[2]
			if feature_id in query_locs:
				num_overlapping_ids += 1
				query_loc = query_locs[feature_id]
				query_chrom, query_start, query_stop = query_loc[0], query_loc[1], query_loc[2]
				if query_chrom == ref_chrom:
					if (query_start >= ref_start and query_start <= ref_stop) or (ref_start >= query_start and ref_start <= query_stop):
						# Is this sufficient to check for overlapping regions?
						num_same_loc += 1
					else:
						num_same_chrom_diff_loc += 1
						print(feature_id, ref_chrom, ref_start, ref_stop, query_chrom, query_start, query_stop, file=diffs_out, sep="\t")
				else:
					num_different_chroms += 1
					print(feature_id, ref_chrom, ref_start, ref_stop, query_chrom, query_start, query_stop, file=diffs_out, sep="\t")
			else:
				num_ref_ids_only += 1
				print(feature_id, ref_chrom, ref_start, ref_stop, file=ref_out, sep="\t")

	with open(query_out_file_name, 'w') as query_out:
		for feature_id, loc in query_locs.items():
			query_chrom, query_start, query_stop = loc[0], loc[1], loc[2]
			if feature_id not in ref_locs:
				num_query_ids_only += 1
				print(feature_id, query_chrom, query_start, query_stop, file=query_out, sep="\t")

	print("\tNumber of features in both files: ", num_overlapping_ids)
	print("\t\tSame locations: ", num_same_loc)
	print("\t\tSame chroms but different locations: ", num_same_chrom_diff_loc)
	print("\t\tDifferent chroms: ", num_different_chroms)
	print("\tNumber of features in reference only: ", num_ref_ids_only)
	print("\tNumber of features in query only: ", num_query_ids_only)



if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('ref_gff3_file')
	parser.add_argument('query_gff3_file')
	opts = parser.parse_args()
	print("> Reading in reference file...")
	ref_gene_locs, ref_tx_locs = read_gff3_file(opts.ref_gff3_file)
	print("> Reading in query file...")
	query_gene_locs, query_tx_locs = read_gff3_file(opts.query_gff3_file)
	print("> Comparing genes for both files...")
	compare_feature_locs(ref_gene_locs, query_gene_locs, "genes_only_ref.txt", "genes_only_query.txt", "gene_loc_diffs.txt")
	print("> Comparing transcripts for both files...")
	compare_feature_locs(ref_tx_locs, query_tx_locs, "txs_only_ref.txt", "txs_only_query.txt", "tx_loc_diffs.txt")


# Sample command: 
# python3 gff3_compare.py GCA_018472595.1_genes.gff3 HG00438.1.chm13.gff3 
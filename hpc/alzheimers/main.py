#!/bin/python
import truster
import os

# Setting up file paths
path_parent = os.path.dirname(os.getcwd())
config_path = os.path.join(path_parent, "src", "config")
software_config = os.path.join(config_path, "software.json")
lunarc_config = os.path.join(config_path, "lunarc.json")
cr_index = "/home/steinac/projects/alz/cellranger_index/premRNAREF_SingleCells/GRCh38_premRNA/"
raw = ["/home/steinac/projects/data/alz/"]
gene_gtf = "/home/steinac/projects/annotations/hg38/gencode/gencode.v38.chr_patch_hapl_scaff.annotation.gtf"
te_gtf = "/home/steinac/projects/annotations/hg38/repeatmask/hg38_rmsk_TEtranscripts.gtf"
star_index = "/home/steinac/projects/alz/STAR_index/" 

# Setting up Experiment object
alzheimers = truster.Experiment("alzheimers", lunarc_config, software_config)
alzheimers.register_samples_from_path(raw)

nuclei_dict = dict()
for sample in alzheimers.samples.keys():
	nuclei_dict[sample] = False

# Performing quantification
quantification_dir = os.path.join(path_parent, "1_counts")
print(quantification_dir)
print("Now quantifying...")
alzheimers.quantify(nuclei = nuclei_dict, cr_index = cr_index, outdir = quantification_dir, jobs = 7)

# Getting clusters
biggest_samples = [
'AD10-AK3431_S11_L004',
#'AD13-AK133_S21_L004',
'AD19-AK137_S25_L004',
'AD1-AK141_S13_L003',
'AD20-AK3570_S14_L004',
#'AD21-AK3474_S31_L003', this one had enough cells but a concerningly fat perc_mito tail
'AD2-AK4226_S1_L004',
'AD4-AK148_S32_L004',
'AD5-AK4955_S2_L004',
#'AD6-AK841_S11_L003',
#'AD8-AK129_S17_L004',
'AD9-AK3738_S6_L004',
#'NC11-AK836_S6_L003',
'NC12-AK3444_S4_L002',
'NC14-AK3433_S5_L002',
'NC15-AK3476_S25_L003',
#'NC16-AK4297_S11_L002', this one had enough cells but showed very low RNA counts per cell
'NC17-AK3566_S13_L002',
'NC18-AK3715_S12_L004',
'NC3-AK4232_S7_L004',
'NC7-AK831_S1_L003'
] # Samples with less than 2000 cells removed

for sampleId in list(alzheimers.samples.keys()):
	if sampleId in biggest_samples:
		alzheimers.set_quantification_outdir(sample_id = sampleId, cellranger_outdir = os.path.join(quantification_dir, sampleId))
	else:
		alzheimers.samples.pop(sampleId)

clusters_dir = os.path.join(path_parent, "2_getClusters")
alzheimers.get_clusters_all_samples(clusters_dir, perc_mitochondrial = 10, normalization_method = "CLR", max_size = 190000, max_genes = 3500, jobs = 10)
alzheimers.set_clusters_outdir(clusters_dir)

# Merging samples
merged_dir = os.path.join(path_parent, "3_mergeSamples")
alzheimers.merge_samples(merged_dir, "CLR", max_size = 2400000, res = 0.1, integrate_samples = True)
alzheimers.set_merge_samples_outdir(merged_dir)

indiv_dict = dict()
for i in biggest_samples:
	indiv_dict[i] = [i]

cond_dict = {"AD": [], "NC": []}
for i in biggest_samples:
	if i[0] == "A":
		cond_dict["AD"].append(i)
	else:
		cond_dict["NC"].append(i)

print(indiv_dict)
print(cond_dict)

# Running pipeline
merged_pipeline_dir_cond = os.path.join(path_parent, "3_mergeSamples", "clusterPipeline_cond_filterPCR")

alzheimers.process_clusters(
	mode = "merged", 
	outdir = merged_pipeline_dir_cond, 
	gene_gtf = gene_gtf, 
	te_gtf = te_gtf, 
	star_index = star_index, 
	RAM = 48725506423, 
	jobs = 10, 
	groups = cond_dict
)

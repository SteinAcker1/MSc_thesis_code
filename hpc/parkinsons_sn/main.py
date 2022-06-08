#!/bin/python
import truster
import os

# Heads up: This is the SUBSTANTIA NIGRA script, not the CORTEX script!

# Setting up file paths
path_parent = os.path.dirname(os.getcwd())
config_path = os.path.join(path_parent, "src", "config")
raw = ["/home/steinac/projects/data/ASAP_replacement/substantia_nigra/"]
gene_gtf = "/home/steinac/projects/annotations/hg38/gencode/gencode.v38.chr_patch_hapl_scaff.annotation.gtf"
te_gtf = "/home/steinac/projects/annotations/hg38/repeatmask/hg38_rmsk_TEtranscripts.gtf"
star_index = "/home/steinac/projects/alz/STAR_index/"
cr_index = "/home/steinac/projects/alz/cellranger_index/premRNAREF_SingleCells/GRCh38_premRNA/"
lunarc_config = os.path.join(config_path, "lunarc.json")
software_config = os.path.join(config_path, "software.json")

# Setting up Experiment object
parkinsons = truster.Experiment("parkinsons", lunarc_config, software_config)

parkinsons.register_samples_from_path(raw)
parkinsons.unregister_sample("substantia_nigra") #for some reason it keeps registering a sample called "substantia_nigra" even though no such file path exists

nuclei_dict = dict()
for sample in parkinsons.samples.keys():
    nuclei_dict[sample] = False
print(nuclei_dict)
print([i.sample_id for i in parkinsons.samples.values()])

# Performing quantification
quantification_dir = os.path.join(path_parent, "1_counts")
print("Now quantifying...")
parkinsons.quantify(nuclei = nuclei_dict, cr_index = cr_index, outdir = quantification_dir, jobs = 60)

# Getting clusters
for sampleId in list(parkinsons.samples.keys()):
        parkinsons.set_quantification_outdir(sample_id = sampleId, cellranger_outdir = os.path.join(quantification_dir, sampleId))

clusters_dir = os.path.join(path_parent, "2_getClusters")
print("Now clustering...")
parkinsons.get_clusters_all_samples(clusters_dir, perc_mitochondrial = 10, normalization_method = "CLR", max_size=2000, res = 0.7, jobs = 30)
parkinsons.set_clusters_outdir(clusters_dir)


# Merging samples
merged_dir = os.path.join(path_parent, "3_mergeSamples")
print("Now merging...")
parkinsons.merge_samples(merged_dir, "CLR", integrate_samples = True)
parkinsons.set_merge_samples_outdir(merged_dir)

# Creating group dicts

cond_dict = {"PD": [], "NC": []}
for i in parkinsons.samples.keys():
    if i.split("_")[1] == "PD":
        cond_dict["PD"].append(i)
    else:
        cond_dict["NC"].append(i)

print(cond_dict)

# Running pipeline
merged_pipeline_dir_cond = os.path.join(path_parent, "3_mergeSamples", "clusterPipeline_cond_filterPCR")

parkinsons.process_clusters(
	mode = "merged", 
	outdir = merged_pipeline_dir_cond, 
	gene_gtf = gene_gtf, 
	te_gtf = te_gtf, 
	star_index = star_index, 
	RAM = 48725506423, 
	jobs = 10, 
	groups = cond_dict
)


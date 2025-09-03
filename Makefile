TOP := $(shell git rev-parse --show-toplevel)
STUDY_DIR := $(TOP)/studies/MGYS00006491/MGYA00679207

all:
	python3 scripts/embeddings_clustering.py --study_dir "$(STUDY_DIR)" --n_sequences 10000 --use_mgnify_faa
	python3 scripts/functional_analysis_summary.py --cluster_results "$(STUDY_DIR)/report/mgnify/clustering_results.csv" --interpro_annotations "$(STUDY_DIR)/downloads/interPro.tsv"
	python3 scripts/deep_dive_cluster.py --cluster_results "$(STUDY_DIR)/report/mgnify/clustering_results.csv" --interpro_annotations "$(STUDY_DIR)/downloads/interPro.tsv" --cluster_id 85
	python3 scripts/find_biosensor_candidates.py --analyte "uranium" --cluster_results "$(STUDY_DIR)/report/mgnify/clustering_results.csv" --interpro_annotations "$(STUDY_DIR)/downloads/interPro.tsv"
	python3 scripts/analyze_genomic_context.py --candidates_report "$(STUDY_DIR)/report/mgnify/uranium_sensor_candidates.txt" --contigs_fasta "$(STUDY_DIR)/downloads/contigs.fasta" --interpro_annotations "$(STUDY_DIR)/downloads/interPro.tsv"
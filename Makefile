.PHONY: clean clean_env data lint environment serve_nb sync_data_to_s3 sync_data_from_s3

#################################################################################
# GLOBALS                                                                       #
#################################################################################
SHELL := /bin/bash

BUCKET = s3_bucket
CONDA_ENV_NAME = gmmgff
CONDA_ROOT = $(shell conda info --root)
CONDA_ENV_DIR = $(CONDA_ROOT)/envs/$(CONDA_ENV_NAME)
CONDA_ENV_PY = $(CONDA_ENV_DIR)/bin/python

#################################################################################
# COMMANDS                                                                      #
#################################################################################

serve_nb:
	source activate $(CONDA_ENV_NAME);\
	jupyter notebook --notebook-dir notebooks

environment:
	conda create -n $(CONDA_ENV_NAME) --file requirements.txt --yes
	source activate $(CONDA_ENV_NAME);\
	ipython kernel install --user --name $(CONDA_ENV_NAME) --display-name "$(CONDA_ENV_NAME)"


data:
	source activate $(CONDA_ENV_NAME);\
	python src/python/data/make_dataset.py

clean:
	find . -name "*.pyc" -exec rm {} \;

clean_env:
	source activate $(CONDA_ENV_NAME);\
	rm -rf $$(jupyter --data-dir)/kernels/$(CONDA_ENV_NAME)
	rm -rf $(CONDA_ENV_DIR)

lint:
	flake8 --exclude=lib/,bin/ .

sync_data_to_s3:
	aws s3 sync data/ s3://$(BUCKET)/data/

sync_data_from_s3:
	aws s3 sync s3://$(BUCKET)/data/ data/

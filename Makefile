.PHONY: clean data lint requirements

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROFILE = default
PROJECT_NAME = tcga_subtype_classification
PYTHON_INTERPRETER = python3
R_INTERPRETER = Rscript

ifeq (,$(shell which conda))
	HAS_CONDA=False
else
	HAS_CONDA=True
endif

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Install Python Dependencies
requirements: test_environment
	pip install -U pip setuptools wheel
	pip install -r requirements.txt

## Make Dataset
data: requirements
	$(R_INTERPRETER) src/data/download.R
	$(R_INTERPRETER) src/data/combine_datasets.R
	$(R_INTERPRETER) src/data/convert_to_eset.R
	$(R_INTERPRETER) src/data/convert_eset_to_feather.R
	$(R_INTERPRETER) src/data/generate_dataframe_prelim_subtypes_with_TCGA_data.R

## Train models and make predictions
models: data
	$(R_INTERPRETER) src/models/make_models.R
	$(R_INTERPRETER) src/models/extract_results.R

## Create visualizations
viz: models
	$(PYTHON_INTERPRETER) src/vizualization/make_viz.py

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Lint using flake8
lint:
	flake8 src

## Set up python interpreter environment
create_environment:
ifeq (True,$(HAS_CONDA))
	@echo ">>> Detected conda, creating conda environment."
	conda create --name $(PROJECT_NAME) python=3.6
	@echo ">>> New conda env created. Activating with:\nsource activate $(PROJECT_NAME)"
	@bash -c "source activate $(PROJECT_NAME)"
else
	@echo ">>> Please install conda via anaconda or miniconda."
	@echo ">>> See here: https://conda.io/miniconda.html for more details."
	$(error This package requires a conda installation.)
endif

## Test python environment is setup correctly
test_environment: install_conda_deps
	$(PYTHON_INTERPRETER) test_environment.py
	R CMD javareconf > /dev/null 2>&1 || true
	$(R_INTERPRETER) test_environment.R

## Install conda packages
install_conda_deps:
	conda config --add channels conda-forge defaults r bioconda
	conda install -n $(PROJECT_NAME) --file=requirements_conda.txt


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################



#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := show-help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')

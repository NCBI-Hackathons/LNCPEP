SHELL:=/bin/bash
MINICONDA_sh:=Miniconda3-4.4.10-Linux-x86_64.sh
MINICONDA_sh_url:=https://repo.continuum.io/miniconda/$(MINICONDA_sh)
# MINICONDA_sh_md5:=
CONDA_INSTALL_DIR:=$(shell echo '$$HOME')/conda
CONDA_ACTIVATE:=$(CONDA_INSTALL_DIR)/bin/activate
VAR:=
.PHONY: $(VAR)

# ~~~~~~ INSTALL ~~~~~~ #
$(MINICONDA_sh):
	wget "$(MINICONDA_sh_url)"


dl: $(MINICONDA_sh)

$(CONDA_INSTALL_DIR): dl
	@if [ ! -d "$(CONDA_INSTALL_DIR)" ]; then \
	bash "$(MINICONDA_sh)" -b -p "$(CONDA_INSTALL_DIR)"; \
	else \
	printf "Install dir already exists: %s\nExiting..." "$(CONDA_INSTALL_DIR)"; \
	exit 1; fi


# install conda in the current directory and install the conda-build package to it
install: $(CONDA_INSTALL_DIR)
	source "$(CONDA_ACTIVATE)" && \
	conda install -y conda-build && \
	conda config --add channels defaults && \
	conda config --add channels conda-forge && \
	conda config --add channels bioconda && \
	conda config --add channels anaconda


# ~~~~~ BUILD ~~~~~ #
check: 
	@if [ ! -d "$(VAR)" ]; then echo "ERROR: VAR is not a valid directory; VAR=$(VAR)"; exit 1; fi
	@if [ ! -f "$(VAR)/env.yml" ]; then echo "ERROR: file does not exist: $(VAR)/env.yml"; exit 1; fi

# standard env build methods
env: check
	source "$(CONDA_ACTIVATE)" && \
	conda env create --file "$(VAR)/env.yml" --name "$(VAR)"

remove:
	source "$(CONDA_ACTIVATE)" && \
	conda env remove -n "$(VAR)"

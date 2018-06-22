# `make build VAR=python3`
# `make test VAR=python3`
VAR=
REPO:=ncbi/tupac-lncur
.PHONY: $(VAR)
none:

check-Docker:
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

# ~~~~~~~ SETUP DOCKER CONTAINERS ~~~~~ #
build: check-Docker
	cd $(VAR) && \
	docker build -t $(REPO):$(VAR) .

pull:
	docker pull $(REPO):$(VAR)

# ~~~~~~ TEST CONTAINERS ~~~~~ #
test: build
	docker run --rm -ti $(REPO):$(VAR) bash

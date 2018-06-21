ZIP:=chickspress-20180621T194354Z-001.zip
RELABELDIR:=relabeled
ALL_CHICKSPRESS_FILE:=chickspress.tsv
ALL_CHICKSPRESS_BED:=chickspress.bed
RELABELEDFILES:=
FINDRELABELED:=
TSVS:=
FINDTSVS:=

ifneq ($(FINDTSVS),)
TSVS:=$(shell find chickspress -type f -name '*.tsv')
endif

ifneq ($(FINDRELABELED),)
RELABELEDFILES:=$(shell find $(RELABELDIR) -type f -name '*.tmp')
endif

all: $(ALL_CHICKSPRESS_BED)

chickspress:
	unzip "$(ZIP)"

relabel: chickspress
	@$(MAKE) relabel-recurse FINDTSVS=1

relabel-recurse: $(TSVS)
.PHONY: $(TSVS)

$(TSVS):
	@ { \
	mkdir -p "$(RELABELDIR)" ; \
	file_base="$$(basename "$@" | sed -e 's|.tsv||g')"; \
	tmp_file="$(RELABELDIR)/$${file_base}.tmp" ; \
	echo ">>> Writing udpated file to: $$tmp_file" ; \
	../scripts/paste-col.py -i "$@" -o "$$tmp_file" --header 'Tissue' --value "$${file_base}" ; \
	}

concat: relabel
	@$(MAKE) concat-recurse FINDRELABELED=1

concat-recurse:
	@echo ">>> Creating $(ALL_CHICKSPRESS_FILE)"
	@../scripts/concat-tables.py $(RELABELEDFILES) -o "$(ALL_CHICKSPRESS_FILE)"
.PHONY: $(RELABELEDFILES)

$(ALL_CHICKSPRESS_FILE): concat

$(ALL_CHICKSPRESS_BED): $(ALL_CHICKSPRESS_FILE)
	@echo ">>> Creating $(ALL_CHICKSPRESS_BED)"
	@cat "$(ALL_CHICKSPRESS_FILE)" | awk '{FS = "\t";OFS = "\t" ;print $$5,$$6,$$7,$$0}' | tail -n +2 > $(ALL_CHICKSPRESS_BED)

clean:
	rm -rf chickspress
	rm -rf "$(RELABELDIR)"
	rm -f "$(ALL_CHICKSPRESS_FILE)" "$(ALL_CHICKSPRESS_BED)"

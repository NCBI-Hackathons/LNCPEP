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

all: $(ALL_CHICKSPRESS_BED) split-tsv split-bed

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
	@../scripts/concat-tables.py $(RELABELEDFILES) > "$(ALL_CHICKSPRESS_FILE)"

$(ALL_CHICKSPRESS_FILE): concat

$(ALL_CHICKSPRESS_BED): $(ALL_CHICKSPRESS_FILE)
	@echo ">>> Creating $(ALL_CHICKSPRESS_BED)"
	@cat "$(ALL_CHICKSPRESS_FILE)" | awk '{FS = "\t";OFS = "\t" ;print $$5,$$6,$$7,$$0}' | tail -n +2 > $(ALL_CHICKSPRESS_BED).tmp ; \
	Rscript --vanilla -e 'prot_df <- read.delim(file = "$(ALL_CHICKSPRESS_BED).tmp", header = FALSE, sep = "\\t"); names(prot_df)[1] <- c("Chrom"); conv_df <- read.delim(file = "ncbi_ens_chrs_mapping.out", header = TRUE, sep = "\\t"); df <- merge(x = prot_df, y = conv_df, by.x = "Chrom", by.y = "Chr_Scaffold"); df[["Chrom"]] <- df[["NCBI_ACC"]]; write.table(x = df, file = "$(ALL_CHICKSPRESS_BED)", sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)' && \
	rm -f $(ALL_CHICKSPRESS_BED).tmp


SPLIT_DIR:=split
TSV_TISSUES:=
GET_TSV_TISSUES:=
ifneq ($(GET_TSV_TISSUES),)
TSV_TISSUES:=$(shell cut -f12 $(ALL_CHICKSPRESS_FILE) | sort -u | tr -d '\r' | sed -e 's|$$|.tsv|' -e 's|^|$(SPLIT_DIR)/|')
endif

BED_TISSUES:=
GET_BED_TISSUES:=
ifneq ($(GET_BED_TISSUES),)
BED_TISSUES:=$(shell cut -f15 $(ALL_CHICKSPRESS_BED) | sort -u | tr -d '\r' | sed -e 's|$$|.bed|' -e 's|^|$(SPLIT_DIR)/|')
endif

split-tsv: $(ALL_CHICKSPRESS_FILE)
	@$(MAKE) split-tsv-recurse GET_TSV_TISSUES=1
split-tsv-recurse: $(TSV_TISSUES)

$(TSV_TISSUES):
	@mkdir -p "$(SPLIT_DIR)" ; \
	output_file="$@" ; \
	pattern="$$(basename "$${output_file}" | sed -e 's|.tsv||')" ; \
	echo ">>> Making $${output_file} for pattern $${pattern}" ; \
	head -1 "$(ALL_CHICKSPRESS_FILE)" > "$${output_file}" ; \
	grep -w "$${pattern}" "$(ALL_CHICKSPRESS_FILE)" >> "$${output_file}"

split-bed: $(ALL_CHICKSPRESS_BED)
	@$(MAKE) split-bed-recurse GET_BED_TISSUES=1
split-bed-recurse: $(BED_TISSUES)
$(BED_TISSUES):
	@mkdir -p "$(SPLIT_DIR)" ; \
	output_file="$@" ; \
	pattern="$$(basename "$${output_file}" | sed -e 's|.bed||')" ; \
	echo ">>> Making $${output_file} for pattern $${pattern}" ; \
	grep -w "$${pattern}" "$(ALL_CHICKSPRESS_BED)" >> "$${output_file}"

clean:
	rm -rf chickspress
	rm -rf "$(RELABELDIR)"
	rm -f "$(ALL_CHICKSPRESS_FILE)" "$(ALL_CHICKSPRESS_BED)"
	rm -rf "$(SPLIT_DIR)"

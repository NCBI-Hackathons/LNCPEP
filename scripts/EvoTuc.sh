#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -g genome"
      echo ""

cat <<'EOF'

  -g </path/to/reference genome file>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":hg:" opt; do
  case $opt in
    h)
     usage
     exit 1
      ;;    
    g)
     referencegenome=$OPTARG
      ;;
    \?)
     echo "Invalid option: -$OPTARG" >&2
     exit 1
      ;;
    :)
     echo "Option -$OPTARG requires an argument." >&2
     exit 1
      ;;
  esac
done

for file in *.gtf

do
	new=$(basename $file ".combined.gtf")

	cut -f 9 $file | tr ";" "\n" | grep transcript | sed 's/^ //g' | cut -d " " -f 2 | grep -Ff - "$file" | gffread -w "${new}".fa -g "$referencegenome" - && \

	python /evolinc_docker/get_gene_length_filter.py "${new}".fa "${new}".putative_intergenic.genes.fa && sed -i 's/ .*//' "${new}".putative_intergenic.genes.fa

#	TransDecoder.LongOrfs -t "${new}".putative_intergenic.genes.fa -m 50

#	find . -type f -name longest_orfs.cds -exec cat '{}' \; | cat > "${new}".longest_orfs.fa

#	rm "${new}".fa "${new}".putative_intergenic.genes.fa

done

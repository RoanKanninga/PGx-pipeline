set -o pipefail
set -eu

awk '{print \\$1,\\$4,\\$4+1,\\$2,\\$5,\\$6,\\$2 "___" \\$5 "___" \\$6}' ${bim} > 'crossmap_input.bed'

CrossMap.py bed ${chain_file} 'crossmap_input.bed' 'crossmap_output.bed'
awk '{print \\$7}' 'crossmap_input.bed' | sort > 'input_ids.txt'
awk '{print \\$7}' 'crossmap_output.bed' | sort > 'output_ids.txt'
comm -23 'input_ids.txt' 'output_ids.txt' | awk '{split(\\$0,a,"___"); print a[1]}' > 'excluded_ids.txt'

plink2 \
--bfile "${bed.baseName}" \
--exclude 'excluded_ids.txt' \
--make-bed \
--output-chr MT 
--out 'crossmapped_plink' \
--keep-allele-order \
--update-map 'crossmap_output.bed' 2 4
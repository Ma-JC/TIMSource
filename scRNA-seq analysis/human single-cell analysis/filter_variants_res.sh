#!/bin/bash
mkdir -p sample
for file in ./sample_raw/*.SNP.hg19_multianno.xls; do
	        filename=$(basename "$file")
		        awk -F '\t' '{
				        if ($22 == "." && $23 == "." && ($45 == "." || $45 < 0.01) && ($55 == "." || $55 < 0.01) && ($65 == "." || $65 < 0.01)) {
						            print $0
							            }
							        }' "$file" > "./sample/$filename"
						done


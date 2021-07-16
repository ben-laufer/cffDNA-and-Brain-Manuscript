#!/bin/bash

# GNU parallel for DMRichR
# Ben Laufer

DMR(){

	contrast=$1
	echo ${contrast}

	mkdir ${contrast}
	cp cytosine_reports/*.gz ${contrast}
	mv sample_info_${contrast}.xlsx ${contrast}/sample_info.xlsx
	cd ${contrast}
	
	echo testCovariate is ${testCovariate}
	
	if [[ "${contrast}" == *"GD"* ]]
	then
	samplesource=cffDNA
	cutoff='0.10'
	else 
	samplesource=brain
	cutoff='0.05'
	fi
	
	echo Sample source is ${samplesource} and cutoff is ${cutoff}

	call="Rscript \
	--vanilla \
	/share/lasallelab/programs/DMRichR/DM.R \
	--genome rheMac10 \
	--coverage 1 \
	--perGroup '0.5' \
	--minCpGs 5 \
	--maxPerms 10 \
	--maxBlockPerms 10 \
	--cutoff ${cutoff} \
	--testCovariate Group \
	--GOfuncR FALSE \
	--cores 10"

	echo $call
	eval $call

	rm *.gz

}
export -f DMR

mkdir dmrLogs
cd /share/lasallelab/Ben/Obesity/

module load R/3.6.3

parallel --dry-run --will-cite --results dmrLogs -j 6 "DMR {}" :::: contrasts.txt
parallel --verbose --will-cite --results dmrLogs -j 6 "DMR {}" :::: contrasts.txt

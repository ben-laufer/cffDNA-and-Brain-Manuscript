# GAT Enrichments
# Ben Laufer

# Human ASD

cd /share/lasallelab/Ben/Obesity/meta/GAT/hg38

echo "Sorting bed files"

parallel --dry-run --will-cite "sort -k1,1 -k2,2n {} > {.}_sorted.bed" ::: *.bed
parallel --verbose --will-cite "sort -k1,1 -k2,2n {} > {.}_sorted.bed" ::: *.bed

echo "Done sorting bed files"

GAT(){

contrast=$1

echo "Testing ${contrast} for ASD enrichments"

call="gat-run.py \
--segments=${contrast}_sigRegions_consensus_hg38_sorted.bed \
--annotations=ASD_brain_male_sorted.bed \
--annotations=Rett_brain_female_sorted.bed \
--annotations=Dup15_brain_male_sorted.bed \
--annotations=DS_brain_male_sorted.bed \
--annotations=ASD_placenta_male_sorted.bed \
--annotations=ASD_placenta_female_sorted.bed \
--workspace=${contrast}_regions_consensus_hg38_sorted.bed \
--counter=nucleotide-overlap \
--num-samples=10000 \
--num-threads=25 \
--log=${contrast}.log \
> ${contrast}_results.tsv"

# --isochore-file=hg38isochores_sorted.bed \

echo $call
eval $call

echo "Done testing ${contrast} for ASD enrichments"
echo

}

export -f GAT

module load gat/1.3.4 

parallel --dry-run --will-cite "GAT {}" ::: brains cffDNA
parallel --verbose --will-cite "GAT {}" ::: brains cffDNA

echo "Removing temporary files"
rm *_sorted.bed
echo "Done removing temporary files"

# rheMac10

cd /share/lasallelab/Ben/Obesity/meta/GAT/rheMac10

echo "Sorting bed files"

parallel --dry-run --will-cite "sort -k1,1 -k2,2n {} > {.}_sorted.bed" ::: *.bed
parallel --verbose --will-cite "sort -k1,1 -k2,2n {} > {.}_sorted.bed" ::: *.bed

echo "Done sorting bed files"

echo "Testing for cffDNA enrichment in brain"

module load gat/1.3.4 

call="gat-run.py \
--segments=cffDNA_sigRegions_consensus_rheMac10_sorted.bed \
--annotations=brains_sigRegions_consensus_rheMac10_sorted.bed \
--workspace=cffDNA_regions_consensus_rheMac10_sorted.bed \
--counter=nucleotide-overlap \
--num-samples=10000 \
--num-threads=50 \
--log=GAT.log \
> GAT_results.tsv"

echo $call
eval $call

echo "Done testing for cffDNA enrichment in brain"

echo "Removing temporary files"
rm *_sorted.bed
echo "Done removing temporary files"



#!/bin/bash

# create precision-recall plots for for vg variants calling using
# gold standard vcf.  Directory structure
# of input is important, and explained in ../README.md and callVariants.py
# running is grouped by region to slightly facilitate debugging and resuming

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <out_dir>"
	 exit 1
fi

# vglr lrc_kir is a bad graph.  we can censor it as input in the wildcard arguments below
# to make it disappear from the analysis
GLOBIGNORE="*/lrc_kir/vglr*:*vglr-lrc_kir*"
# leave out debruijn mhc for now as it makes invalid augmented graph
GLOBIGNORE="*/mhc/debruijn*:debruijn*mhc*:${GLOBIGNORE}"
GLOBIGNORE="*/lrc_kir/debruijn*:debruijn*lrc_kir*:${GLOBIGNORE}"
#ditto camel -- vcf conversion doesn't work for crazy regions (probably due to caller making bizarre graph?)
GLOBIGNORE="*/mhc/camel*:camel*mhc*:${GLOBIGNORE}"
GLOBIGNORE="*/lrc_kir/camel*:camel*lrc_kir*:${GLOBIGNORE}"
# leave simons out for now as its hg19
GLOBIGNORE_COMP="*simons*:${GLOBIGNORE}"
GLOBIGNORE_CALL=$GLOBIGNORE

# command line args
GRAPHS=$1
ALIGNMENTS=$2
OUT_DIR=$3

TOIL_DIR=call_pr_toil5

# VCF comparisons : sompy happy vcf vcfeval
#COMPS=( "sompy" "happy" "vcf" )
# graph comparison : kmer corg
#COMPS=( "kmer" "corg")
COMPS=( "happy" )

# Parameters for indexing (used for corg / kmer comparison)
INDEX_OPTS="--kmer 20 --edge_max 5 --timeout 5000"

# Calling parameters
##REGIONS=( "brca2" "mhc" "brca1" "sma" "lrc_kir" )
REGIONS=( "lrc_kir" )
#OPTS="--maxCores 8 --vg_cores 4 --vg_only --depth 10 --skipBaseline"
OPTS="--maxCores 20 --vg_cores 4 --vg_only --depth 10 --skipBaseline"

# Comparison parameters


# Normalization (requires vt, recommended for som.py and vcf comparison)
#COMP_OPTS="--clip data/filters/platinum.bed --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 30"

# No normalization (only recommended for vcfeval and hap.py)
COMP_OPTS="--vg_cores 4 --maxCores 20 --gt"

# example how to use user-specified platypus and freebayes vcfs
#COMP_OPTS="--clip data/filters/platinum.bed --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 30 --platypus_path platinum_classic/platypus --freebayes_path platinum_classic/freebayes"

COMP_TAG=comp.gt

# call variants, compute and plot baseline comparison
# note how sample name (NA12878) is hardcoded.  Use wildcard to do all available samples
function run_pipeline {

	 VARIANTS_OUT_DIR=$1
	 CALL_OPTS=$2
	 PILEUP_OPTS=$3
	 FILTER_OPTS=$4
	 ROC=$5
	 	 
	 for i in "${REGIONS[@]}"
	 do
		  GLOBIGNORE=$GLOBIGNORE_CALL
		  mkdir ${VARIANTS_OUT_DIR}
		  #rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG}  ${ALIGNMENTS}/${i}/*/NA12878.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2>> ${VARIANTS_OUT_DIR}/call_log_${i}.txt

		  for j in "${COMPS[@]}"
		  do
				GLOBIGNORE=$GLOBIGNORE_COMP				
				# compute distances
				mkdir ${VARIANTS_OUT_DIR}.${COMP_TAG}
				rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/${i}/*/NA12878.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${VARIANTS_OUT_DIR}.${COMP_TAG} ${COMP_OPTS} ${INDEX_OPTS} ${ROC}  2>> ${VARIANTS_OUT_DIR}.${COMP_TAG}/comp_log_${i}_${j}.txt
				GLOBIGNORE=$GLOBIGNORE_CALL
		  done

	 done
}

# 1000 Genomes-like options fixed here for pileups
PILEUP_OPTS=" -w 40 -m 10 -q 10 "

# filter options (no secondary, primary must have > 90% identity and > 5% improvement over secondary)
mkdir ${OUT_DIR}
ident=0.90
delta=0.05
secscore=10000

# points on the precision-precall chart hard-coded here. 
depths=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 22 24 26 28 30 35 40 45)
supports=(01 02 03 04 05 06 07 08 09 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 15 15 15 15)
quals=(0.000 0.020 0.040 0.060 0.080 0.010 0.100 0.120 0.140 0.160 0.180 0.200 0.220 0.240 0.260 0.280 0.300 0.320 0.340 0.360 0.380 0.400 0.500 0.600 0.700 0.800 0.900 1.000)
#for i in {0..27}
#for i in 1 2 3 4 5 6 7 8 1 10 20
for i in 1 2 5 10 12 15 20 23
do
	 depth=${depths[${i}]}
	 qual=${quals[${i}]}
	 support=${supports[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore} "-r 0.0001 -b 0.4 -f 0.25 -s ${support} -d ${depth}" "$PILEUP_OPTS" "-r ${ident} -d ${delta} -e ${delta} -afu -s ${secscore} -o 10 " "$ROC_FLAG"
done
wait

# scrape together all results into one tsv / region / comp type
scripts/rocDistances.py ${OUT_DIR}/primary_*_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/pr_plots.${COMP_TAG} --best_comp happy

# finally, draw out the tables created above
scripts/plotVariantsDistances.py ${OUT_DIR}/pr_plots.${COMP_TAG}

exit 0


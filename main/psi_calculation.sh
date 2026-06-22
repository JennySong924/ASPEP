### calculate junction psi for a single cell type


## input -----------------------------------------

while getopts "f:c:o:n:j:" ARGS
do
        case $ARGS in
                f )
                    file=$OPTARG;;  ## junc count file
                c )
                    coveragecut=$OPTARG;;  ## junction score cutoff, keep junctions with score >= cutoff                    
                o )
                    outdir=$OPTARG;;  ## outdir folder
                    # cd $outdir;;
                n )
                    basename=$OPTARG;;  ## prefix
                ? )
                    echo "error"
                    exit 1;;
                esac
done


echo "outdir directory:" $outdir
mkdir -p $outdir

## filter junctions and create junction list file ----------------------------------------------
## junction file should contains: chrom   start   end     name    score   strand  splice_site     
## acceptors_skipped       exons_skipped   donors_skipped  anchor  known_donor     known_acceptor  
## known_junction  gene_names      gene_ids        transcripts (rmats outdir)

echo "Current directory: $(pwd)"
echo "File variable: $file"
ls -l "$file"

echo 'filter junctions and create junction list file...'
mkdir -p "$outdir"

awk -v s="$s" '$5 >= s' "$file" > "${outdir}/junction_score_filtered"
awk 'FNR > 1 {print $1, $2, $3, $6, $15, $16}' OFS="\t" "${outdir}/junction_score_filtered" \
  | sort -k1,1 -k2,2n -k3,3n -u \
  | awk -F"\t" '{print $0 "\t" $1 ":" $2 "-" $3 ":" $4}' \
  > "${outdir}/unique_junction_score_filtered"


## calculate bp depth ------------------------------------------------

SORTED_BED=${outdir}/${basename}_sorted_bed
BREAKPOINTS=${outdir}/${basename}_break_point
SEGMENTS=${outdir}/${basename}_segments

# 1. order bed file
cut -f1-6 "$file" |sort -k1,1 -k2,2n -k3,3n > "$SORTED_BED"

awk 'BEGIN{OFS="\t"} { $2=$2-1; print }' $SORTED_BED > ${SORTED_BED}_base0

# awk '$6=="+"' ${SORTED_BED}_base0 > ${SORTED_BED}_plus_base0
# awk '$6=="-"' ${SORTED_BED}_base0 > ${SORTED_BED}_minus_base0


PLUS_FILE="${SORTED_BED}_plus_base0"
MINUS_FILE="${SORTED_BED}_minus_base0"


> "$PLUS_FILE"
> "$MINUS_FILE"


awk -F"\t" '{
    if($6=="+") print >> "'"$PLUS_FILE"'";
    else if($6=="-") print >> "'"$MINUS_FILE"'";
}' ${SORTED_BED}_base0


# 2. collect starts and ends
awk '{print $1"\t"$2; print $1"\t"$3}' ${SORTED_BED}_base0 | sort -k1,1 -k2,2n | uniq > $BREAKPOINTS

# 3. divide into segments
awk '{
    chr=$1; pos=$2;
    if(chr != prev_chr){
        prev_chr = chr
        prev_pos = pos
        next
    }

    seg_start = prev_pos
    seg_end = pos
    if(seg_start <= seg_end){
        printf "%s\t%d\t%d\n", chr, seg_start, seg_end
    }
    prev_pos = pos
}' $BREAKPOINTS > $SEGMENTS


awk 'BEGIN{OFS="\t"} { $2=$2; print }' $SEGMENTS > ${SEGMENTS}_0base

# 4. calculate depth for segments
bedtools map -a ${SEGMENTS}_0base -b ${SORTED_BED}_plus_base0 -c 5 -o sum -null 0 > ${outdir}/${basename}.depth.plus0
bedtools map -a ${SEGMENTS}_0base -b ${SORTED_BED}_minus_base0 -c 5 -o sum -null 0 > ${outdir}/${basename}.depth.minus0
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1":"$2"-"$3, $4, "+"}' ${outdir}/${basename}.depth.plus0 > ${outdir}/${basename}.depth.plus
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1":"$2"-"$3, $4, "-"}' ${outdir}/${basename}.depth.minus0 > ${outdir}/${basename}.depth.minus
cat ${outdir}/${basename}.depth.plus ${outdir}/${basename}.depth.minus | sort -k1,1 -k2,2n -k3,3n > ${outdir}/${basename}.depth



# 5. remove temporal files
rm "$SORTED_BED" "$BREAKPOINTS" "$SEGMENTS" ${SORTED_BED}_base0 ${SEGMENTS}_0base
rm ${outdir}/${basename}.depth.plus ${outdir}/${basename}.depth.minus ${outdir}/${basename}.depth.plus0 ${outdir}/${basename}.depth.minus0
rm ${SORTED_BED}_plus_base0 ${SORTED_BED}_minus_base0



## junction ends depth preparation ------------------------------------------------

awk '{print $1, $2, $6; print $1, $3, $6}' $file | sort -k1,1 -k2,2n -u | awk '{print $1, $2-1, $2, ".", ".", $3}' OFS="\t" > $outdir/junction_ends

bedtools intersect -a $outdir/junction_ends -b ${outdir}/${basename}.depth -wa -wb -s > $outdir/${basename}_junction_ends_depth


## calculate psi ----------------------------------------------
Rscript psi_calculation.R $outdir/unique_junction_score_filtered $file $outdir $basename $coveragecut


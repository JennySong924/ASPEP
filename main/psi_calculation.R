if (!requireNamespace("data.table", quietly = TRUE)){
  install.packages("data.table")
}
if (!requireNamespace("dplyr", quietly = TRUE)){
  install.packages("dplyr")
}
suppressMessages(require(data.table))
suppressMessages(require(dplyr))
options(stringsAsFactors = F)
warnings('off')

## input -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide the path to your input and output file")
}

full_junc_list <- args[1] 
input_arg <- args[2]  
output_arg <- args[3] 
basename <- args[4]
coverage_cut <- as.integer(args[5])

## functions -------------------------------------------------------

weighted_psi <- function(data,coverage_cut){
  data[, start_psi := fifelse(start_coverage>=coverage_cut, score/start_coverage, NaN)]
  data[, end_psi   := fifelse(end_coverage>=coverage_cut,   score/end_coverage,   NaN)]
  data[, start_psi_weight := fifelse((start_coverage+end_coverage)>0,
                                     start_coverage/(start_coverage+end_coverage), NaN)]
  data[, start_psi_weight := fifelse(is.na(end_psi),1, start_psi_weight)]
  data[, end_psi_weight   := fifelse((start_coverage+end_coverage)>0,
                                     end_coverage/(start_coverage+end_coverage), NaN)]
  data[, end_psi_weight := fifelse(is.na(start_psi),1, end_psi_weight)]
  
  data[, start_psi_weighted := start_psi * start_psi_weight]
  data[, end_psi_weighted   := end_psi * end_psi_weight]
  
  data[, weighted_psi := rowSums(.SD, na.rm = TRUE), .SDcols = c("start_psi_weighted","end_psi_weighted")]
  # data$weighted_psi[data$start_coverage==0 & data$end_coverage==0] <- NaN
  data$weighted_psi[is.na(data$start_psi) & is.na(data$end_psi)] <- NaN
  
  return(data)
}

## main -------------------------------------------------------
cat(paste0('Now working on ',basename,'\n'))

cat('Loading input files...\n')

depth_path <- paste0(output_arg,'/depth/',basename)
full_junc_list_df <- as.data.frame(fread(full_junc_list))
junc.c <- as.data.frame(fread(input_arg))

cat('Calculating junction ends depth...\n')


chr    <- junc.c[[1]]
start  <- junc.c[[2]]
end    <- junc.c[[3]]
strand <- junc.c[[6]]


end1 <- data.frame(chr = chr, pos = start, strand = strand)
end2 <- data.frame(chr = chr, pos = end,   strand = strand)


junction_ends <- rbind(end1, end2)


junction_ends <- junction_ends[order(junction_ends$chr,
                                     junction_ends$pos,
                                     junction_ends$strand), ]
junction_ends <- unique(junction_ends)


bed_df <- data.frame(
  chr   = junction_ends$chr,
  start = pmax(junction_ends$pos - 1, 0),  # 防止出现负数
  end   = junction_ends$pos,
  name  = ".",
  score = ".",
  strand = junction_ends$strand,
  stringsAsFactors = FALSE
)


out_file <- paste0(depth_path, "_junction_ends")
write.table(bed_df, file = out_file,
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

system(paste0("bedtools intersect -a ",depth_path,"_junction_ends"," -b ",depth_path,'.depth'," -wa -wb -s > ",depth_path,"_junction_ends_depth"))


ends.depth <- as.data.frame(fread(paste0(output_arg,"/",basename,"_junction_ends_depth")))


cat('Calculating psi...\n')

full_junc_list_df_psi <- left_join(full_junc_list_df,ends.depth[,c(1,3,6,11)],by=c('V1','V2'='V3','V4'='V6'))
full_junc_list_df_psi <- left_join(full_junc_list_df_psi,ends.depth[,c(1,3,6,11)],by=c('V1','V3','V4'='V6'))
colnames(full_junc_list_df_psi) <- c('chr','start','end','strand','gene_names','gene_ids','junction_id','start_coverage','end_coverage')
full_junc_list_df_psi$start_coverage[is.na(full_junc_list_df_psi$start_coverage)] <- 0
full_junc_list_df_psi$end_coverage[is.na(full_junc_list_df_psi$end_coverage)] <- 0
full_junc_list_df_psi <- left_join(full_junc_list_df_psi,junc.c[,c(1,2,3,5,6)],by=c('chr'='chrom','start','end','strand'))
full_junc_list_df_psi$score[is.na(full_junc_list_df_psi$score)] <- 0
setDT(full_junc_list_df_psi)
output_df <- weighted_psi(full_junc_list_df_psi,coverage_cut)
output_df <- data.frame(output_df)

cat('Saving output...\n')  

psi_path <- paste0(output_arg,'/',basename,'.psi.gz') 
fwrite(output_df,psi_path,quote=F,col.names=T,row.names=F,sep='\t') 

cat('Done\n')

## junction_id

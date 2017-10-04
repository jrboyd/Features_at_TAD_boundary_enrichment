library(data.table)
source('functions.R')
if(exists("query_grs")) remove(query_grs)
#set file names
tad_bound_file = "data/MCF10A_tad_boundaries.txt"
bin_background_file = "data/40kb_bins.bed"
gencode_gtf_file = "data/gencode.v25.gene.annotation.gtf"
np_dir = "data"
IS_MALE = F
n_iterations = 10

###STEP 1 load GRanges objects to check vs TADs
###Two Strategies:
###1)loading peak from MACS2 narrowPeak
###2)loading annotation from gtf

#1) load peak files
#all are stored in query_grs
#you could easily point this at a different directory and load all the data
np_files = dir(np_dir, full.names = T, pattern = "narrowPeak$")
np_names = sapply(strsplit(basename(np_files), "_"), function(x)x[2])
for(i in 1:length(np_files)){
  add_query_gr_from_narrowPeak(np_file = np_files[i], gr_name = np_names[i], toss_chrY = !IS_MALE)
}

#2) load annotation from gencode gtf
source("parse_gtf.R")
gencode_dt = as.data.table(parse_gtf(gencode_gtf_file, additional_attrib = "gene_type"))
#GRanges compatible colname
colnames(gencode_dt)[colnames(gencode_dt) == "chrm"] = "seqnames"
#convert full transcript region to tss strand sensitively
gencode_dt[strand == "+", end := start]
gencode_dt[strand == "-", start := end]
#subset by gene_type
gr_pc = GRanges(gencode_dt[gene_type == "protein_coding"])
gr_lincs = GRanges(gencode_dt[gene_type %in% c("lincRNA", "macro_lncRNA", "bidirectional_promoter_lncRNA")])
#add to query_grs
query_grs = append(query_grs, list("tss:protein_coding" = gr_pc, "tss:lincs" = gr_lincs))

###STEP 2 load TAD boundaries and all background bins
#load tad boundaries and bins backgroud
tad_bounds_dt = fread(tad_bound_file)
tad_bounds_gr = GRanges(tad_bounds_dt)
#ASSUMPTION bin_size is equal to size of first tad boundary
bin_size = tad_bounds_dt[1, end - start]
all_bins_dt = fread(bin_background_file)
colnames(all_bins_dt) = c("seqnames", "start", "end", "index")
all_bins_gr = GRanges(all_bins_dt)

###STEP 3 cleanup chromosomes
all_chrm = paste0("chr", c(1:22, "X"))
if(IS_MALE) all_chrm = c(all_chrm, "chrY")
query_grs = lapply(query_grs, function(x){
  subset(x, seqnames %in% all_chrm)
})
tad_bounds_gr = subset(tad_bounds_gr, seqnames %in% all_chrm)
all_bins_gr = subset(all_bins_gr, seqnames %in% all_chrm)


###STEP 4 comparison test and plot
nc = 3
len = ceiling(length(query_grs) / nc) * nc
layout(matrix(1:len, ncol = nc))
calc_tad_overlaps_vs_random(query_grs = query_grs, 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = n_iterations)

###FOLLOW-UP derive new interval sets from intersection/exclusions
layout(1:3)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "RUNX1", "CTCF"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = n_iterations)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "RUNX1", "tss:protein_coding"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = n_iterations)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "RUNX1", "tss:lincs"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = n_iterations)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "H3K4ME3", "H3K4AC"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = n_iterations)


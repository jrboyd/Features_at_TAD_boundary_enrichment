library(data.table)
source('functions.R')
if(exists("query_grs")) remove(query_grs)
#load peak files
#all are stored in query_grs
#you could easily point this at a different directory and load all the data
np_files = dir("data", full.names = T, pattern = "narrowPeak$")
np_names = sapply(strsplit(basename(np_files), "_"), function(x)x[2])
for(i in 1:length(np_files)){
  add_query_gr_from_narrowPeak(np_file = np_files[i], gr_name = np_names[i], toss_chrY = T)
}

#load annotation from gencode gtf
source("parse_gtf.R")
gencode_dt = as.data.table(parse_gtf("data/gencode.v25.gene.annotation.gtf", additional_attrib = "gene_type"))
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

#load tad boundaries and bins backgroud
tad_bounds_dt = fread("data/MCF10A_tad_boundaries.txt")
tad_bounds_gr = GRanges(tad_bounds_dt)
#ASSUMPTION bin_size is equal to size of first tad boundary
bin_size = tad_bounds_dt[1, end - start]
all_bins_dt = fread("data/40kb_bins.bed")
colnames(all_bins_dt) = c("seqnames", "start", "end", "index")
all_bins_gr = GRanges(all_bins_dt)

#test and plot
nc = 3
len = ceiling(length(query_grs) / nc) * nc
layout(matrix(1:len, ncol = nc))
calc_tad_overlaps_vs_random(query_grs = query_grs, 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = 10)

#example intersection operations
layout(1:3)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "RUNX1", "CTCF"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = 10)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "RUNX1", "tss:protein_coding"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = 10)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "RUNX1", "tss:lincs"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = 10)
calc_tad_overlaps_vs_random(query_grs = compare_peak_sets(query_grs, "H3K4ME3", "H3K4AC"), 
                            tads_gr = tad_bounds_gr, 
                            all_hic_bins = all_bins_gr, n_iterations = 10)

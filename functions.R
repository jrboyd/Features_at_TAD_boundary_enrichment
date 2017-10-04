

#efficiently randomizes hic_1d (avoids creating 1 GRange per iteration)
#only valid for even bins across genomes - like with hic bins
#returns number of overlaps of qgr with each
#still effing slow
efficient_iterate_olap_random_gr = function(qgr, gr_src, gr_target, n_iter = 1000, seed = 0){
  require(data.table)
  require(GenomicRanges)
  set.seed(seed)
  print("randomize...")
  rand_i = unlist(pblapply(1:n_iter, function(x)sample(x = 1:length(gr_src), size = length(gr_target))))
  print("find overlaps (no progress bar)...")
  hits_dt = as.data.table(findOverlaps(qgr, gr_src[rand_i]))
  setkey(hits_dt, subjectHits)
  print("report...")
  pbsapply(1:n_iter, function(i){
    ir = i * length(gr_target) + c(1 - length(gr_target):0)
    hits_dt[.(ir)][!is.na(queryHits), .N]
  })
}

#calculate enrichment of query in subject compared to supplied randomly generated equivalent regions
calc_tad_overlaps_vs_random = function(query_grs, tads_gr, all_hic_bins, n_iterations = 1000){
  require(GenomicRanges)
  require(pbapply)
  bin_w = table(width(all_hic_bins))
  all_hic_mode = names(bin_w)[bin_w == max(bin_w)]
  bin_w = table(width(tads_gr))
  tads_mode = names(bin_w)[bin_w == max(bin_w)]
  if(!all_hic_mode == tads_mode) stop("problem with bin widths mismatch")
  
  subject_hit_counts = lapply(query_grs, function(qgr){
    length(unique(queryHits(findOverlaps(qgr, tads_gr))))
  })
  random_hit_counts = lapply(query_grs, function(qgr){
    efficient_iterate_olap_random_gr(qgr, gr_src = all_hic_bins, gr_target = tads_gr, n_iter = n_iterations)
  })
  
  cbind(sapply(1:length(query_grs), function(i){
    grp_name = names(query_grs)[i]
    subj_hit = subject_hit_counts[[i]]
    rand_hit = random_hit_counts[[i]]
    
    u = format(round(mean(rand_hit), 2), nsmall = 2)
    fe = format(round(subj_hit / mean(rand_hit), 2), nsmall = 2)
    txt_summary = paste(subj_hit, "/", u, "=", fe, "fold-enrichment :", grp_name)
    
    hist(rand_hit, main = paste(grp_name, ":", fe, "fold-enrichment"),
         xlab = "Hit Count",
         xlim = c(min(rand_hit, 0), max(rand_hit, subj_hit)*1.2))
    ymax = par('usr')[4]
    lines(rep(subj_hit,2), c(0, ymax), col = "red")
    par(xpd = NA)
    text(mean(rand_hit), ymax, "Expected", adj = c(.5,-.2))
    text(subj_hit, ymax, "Observed", adj = c(.5,-.2))
    return(txt_summary)
  }))
}

#create gr from input narrowPeak file, np_file. 
#add gr to globally environment query_grs list.
#create list if necessary
#can discard chrY if desired
add_query_gr_from_narrowPeak = function(np_file, gr_name = NULL, toss_chrY = F){
  require(GenomicRanges)
  #load and convert narrowPeak file
  np = read.table(np_file)
  cn = c("seqnames", "start", "end", "peak_id", "score", "dot", "fe", "p", "q", "width")
  colnames(np) = cn
  if(toss_chrY) np = subset(np, seqnames != "chrY")
  gr = GRanges(np)
  #initialize query_grs if needed
  if(!exists("query_grs")) query_grs <<- list()
  #create gr_name if needed
  if(is.null(gr_name)){
    gr_name = paste0("gr_", length(query_grs))
  }
  #add to query_grs
  toadd = list(gr)
  names(toadd) = gr_name
  query_grs <<- append(query_grs, toadd)
  return()
}

compare_peak_sets = function(gr_list, a_name, b_name){
  a = gr_list[[a_name]]
  elementMetadata(a) = NULL
  b = gr_list[[b_name]]
  elementMetadata(b) = NULL
  olaps = findOverlaps(a, b)
  
  a_and_b = reduce(c(a[queryHits(olaps)], b[subjectHits(olaps)]))
  a_not_b = a[-queryHits(olaps)]
  b_not_a = b[-subjectHits(olaps)]
  int_grs = list(a_and_b, a_not_b, b_not_a)
  names(int_grs) = c(paste(a_name, "and", b_name),
                     paste(a_name, "not", b_name),
                     paste(b_name, "not", a_name))
  return(int_grs)
}

parse_gtf = function(gtf_name, rownames_attrib = 'gene_id', feature_type = 'gene', additional_attrib = c()) {
  print('loading gtf contents...')
  raw_lines = read.table(gtf_name, sep = "\n", stringsAsFactors = F)[, 1]
  get_col = function(i) {
    unlist(lapply(strsplit(raw_lines, "\t"), function(x) x[[i]]))
  }
  print('filtering gtf contents...')
  types = get_col(3)
  keep = types == feature_type
  raw_lines = raw_lines[keep]
  
  attribs = get_col(9)
  all_attribs = strsplit(attribs, "; +")
  
  
  get_attrib = function(key) {
    out = lapply(all_attribs, function(x) {
      keep = grepl(key, x)
      str = "NA"
      if (sum(keep) > 0) 
        str = strsplit(x[keep], " ")[[1]][2]
      return(str)
    })
    return(unlist(out))
  }
  
  print('parsing gtf attributes')
  gene_id = get_attrib("gene_id")
  gene_name = get_attrib("gene_name")
  chrm = get_col(1)
  strand = get_col(7)
  start = as.numeric(get_col(4))
  end = as.numeric(get_col(5))
  rnames = get_attrib(rownames_attrib)
  if(sum(duplicated(rnames)) > 0){
    warnings("the rowname_attrib was not unique, using arbitraty number instead. rowname_attrib will be included as column.")
    additional_attrib = c(additional_attrib, rownames_attrib)
    rnames = 1:length(gene_name)
  }
  ref_dict = data.frame(gene_id, gene_name, chrm, start, end, strand, row.names = rnames, 
                        stringsAsFactors = F)
  for(attrib in additional_attrib){
    ref_dict = cbind(ref_dict, get_attrib(attrib))
    colnames(ref_dict)[ncol(ref_dict)] = attrib
  }
  
  return(ref_dict)
  # plot.pie = function(dat){ uniq = unique(dat) tmp = sapply(uniq, function(x)return(sum(x == dat)))
  # pie(tmp) } plot.pie(support) enst2support = support names(enst2support) = enst_id enst2ensg = ensg_id
  # names(enst2ensg) = enst_id
}

find_gene_in_gtf = function(gtf_name, rownames_attrib = 'gene_id', target_gene = 'CASC7', additional_attrib = c()) {
  print('loading gtf contents...')
  raw_lines = read.table(gtf_name, sep = "\n", stringsAsFactors = F)[, 1]
  get_col = function(i) {
    unlist(lapply(strsplit(raw_lines, "\t"), function(x) x[[i]]))
  }
  print('filtering gtf contents...')
  
  #   keep = types == feature_type
  #   raw_lines = raw_lines[keep]
  
  attribs = get_col(9)
  all_attribs = strsplit(attribs, "; ")
  
  
  get_attrib = function(key) {
    out = lapply(all_attribs, function(x) {
      keep = grepl(key, x)
      str = "NA"
      if (sum(keep) > 0) 
        str = strsplit(x[keep], " ")[[1]][2]
      return(str)
    })
    return(unlist(out))
  }
  
  print('parsing gtf attributes')
  types = get_col(3)
  gene_id = get_attrib("gene_id")
  gene_name = get_attrib("gene_name")
  chrm = get_col(1)
  strand = get_col(7)
  start = as.numeric(get_col(4))
  end = as.numeric(get_col(5))
  
  
  
  ref_dict = data.frame(gene_id, types, gene_name, chrm, start, strand, end,  
                        stringsAsFactors = F)
  for(attrib in additional_attrib){
    ref_dict = cbind(ref_dict, get_attrib(attrib))
    colnames(ref_dict)[ncol(ref_dict)] = attrib
  }
  
  keep = gene_name 
  ref_dict = ref_dict[keep,]
  
  return(ref_dict)
  # plot.pie = function(dat){ uniq = unique(dat) tmp = sapply(uniq, function(x)return(sum(x == dat)))
  # pie(tmp) } plot.pie(support) enst2support = support names(enst2support) = enst_id enst2ensg = ensg_id
  # names(enst2ensg) = enst_id
}

if (F) {
  parse_gtf("mycounts_annotation.gff")
  save(ref_dict, file = "ref_dicts.save")
} 
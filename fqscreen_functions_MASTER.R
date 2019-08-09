#
#
#
# ------------- FUNCTION | Helper for naming read ID ----------- #
find_read_id = function(sample_name) {
  lkup_id = grep("(R\\d+)", sample_name)
  if(length(lkup_id) == 0) {
    cat("Nothing to grep here.", lkup_id, sep = "\n")
    read_id = "Single-End"
  } else {
    cat("You grepped it!", lkup_id, sep = "\n")
    read_id = gsub("(R\\d+)(*SKIP)(*FAIL)|.", "", sample_name, perl = TRUE)
  }
  return(read_id)
}

#
#
#
# ----------- FUNCTION | Get complete plotting df from all files ------------ #
df_from_all_files = function(all_files) {
  #Read in all files
  files = lapply(all_files$datapath, read.delim,
                 sep = "\t", header = FALSE,
                 stringsAsFactors = FALSE)
  
  # ----------- OPERATIONS ON EACH INDIVIDUAL FILE ----------- #
  # Remove dirty first row
  files = lapply(files, function(x) x[-c(1),])
  # Rename the df columns appropriately
  data_cols = c("genome", "reads_processed",
                "num_unmapped", "percent_unmapped",
                "num_one_hit_one_genome", "percent_one_hit_one_genome",
                "num_multiple_hits_one_genome", "percent_multiple_hits_one_genome", 
                "num_one_hit_multiple_genomes", "percent_one_hit_multiple_genomes",
                "num_multiple_hits_multiple_genomes", "percent_multiple_hits_multiple_genomes")
  files = lapply(files, setNames, data_cols)
  # Remove old header
  files = lapply(files, function(x) x[-c(1),])
  # Clean No Hits to fit into df as column
  files = lapply(files, function(x) rbind(head(cbind(x, percent_no_hits = 0), 14),
                                          rep(c("No hits", 0, as.numeric(gsub("[^.0-9]+", "", x[15,1]))), c(1, 11, 1))))
  # Make columns numeric (except genomes)
  files = lapply(files, function(x) {
    x[,-c(1)] = lapply(x[,-c(1)], as.numeric)
    x
  })
  # Lists for additional info needed in df
  #sample_name = sapply(gsub("(_screen.txt)", "", all_files$name), `[`, USE.NAMES = FALSE)
  sample_name = sapply(gsub("(_screen.txt)", "", all_files$name), `[`, USE.NAMES = FALSE)
  #sample = sapply(gsub("_R\\d+", "", sample_name), `[`, USE.NAMES = FALSE)
  sample = sapply(gsub("(_R\\d+)|(_S\\d+)", "", sample_name, perl = TRUE), `[`, USE.NAMES = FALSE)
  run_id = sapply(gsub("(S\\d+)(*SKIP)(*FAIL)|.", "", sample_name, perl = TRUE), `[`, USE.NAMES = FALSE)
  #read = sapply(gsub(".+?[?=_R\\d+]", "", sample_name), `[`, USE.NAMES = FALSE)
  #read = sapply(gsub("(R\\d+)(*SKIP)(*FAIL)|.", "", sample_name, perl = TRUE), `[`, USE.NAMES = FALSE)
  read = sapply(sample_name, find_read_id)
  
  genome_ord = c("Human", "Mouse", "rRNA", "MT", "No hits", "Adapters", "Rat", "Drosophila", "Ecoli", 
                 "Worm", "Yeast", "Arabidopsis", "Lambda", "PhiX", "Vectors")
  
  hits = c("One Hit/One Genome", "Multiple Hits/One Genome",
           "One Hit/Multiple Genomes", "Multiple Hits/Multiple Genomes")
  
  per_mapping = c(hits, "No Hits")
  
  num_mapping = c("Genome", "Reads Processed", "Unmapped", hits)
  
  # Create df for plotting
  df_per = lapply(files, function(in_dat) {
    idx_per = grep("percent", colnames(in_dat))
    per_dat = in_dat[c(1, idx_per[-c(1)])]
    df_per = melt(per_dat, id.vars = c("genome"))
    df_per = cbind(df_per, genome_ord = factor(df_per$genome, levels = genome_ord))
    df_per$variable = factor(df_per$variable, labels = per_mapping)
    df_per
  })
  # Add IDs
  df_per = Map(cbind, df_per, sample_name = sample_name, sample = sample, run_id = run_id, read = read)
  
  # Create df for number mapping
  df_num = lapply(files, function(in_dat) {
    idx_num = grep("^num", colnames(in_dat))
    num_dat = in_dat[c(1:2, idx_num)]
    df_num = setNames(num_dat, num_mapping)
    df_num = cbind(df_num, genome_ord = factor(df_num$Genome, levels = genome_ord))
    df_num
  })
  # Add IDs
  df_num = Map(cbind, df_num, sample_name = sample_name, sample = sample, run_id = run_id, read = read)
  
  # Change factors to character
  char_cols = c("sample_name", "sample", "run_id", "read")
  
  # df_per
  df_per = lapply(df_per, function(df){
    df[char_cols] = lapply(df[char_cols], as.character)
    df
  })
  
  # df_num
  df_num = lapply(df_num, function(df){
    df[char_cols] = lapply(df[char_cols], as.character)
    df
  })
  
  
  # ------------- BIND ALL FILES TO MASTER DFS ------------- #
  df_per = bind_rows(df_per)
  df_num = bind_rows(df_num)
  # Order samples numerically
  sample_levels = unique(df_per$sample[order(nchar(df_per$sample), df_per$sample)])
  cat("Sample names:", sample_levels, sep = "\n")
  df_per$sample_ord = factor(df_per$sample, levels = sample_levels)
  df_num$sample_ord = factor(df_num$sample, levels = sample_levels)
  
  # Return complete df
  return(list(df_per, df_num))
}

#
#
#
# -------------- FUNCTION | Plot data per sample ------------- #
bar_plot_per_sample = function(data, sample_name) {
  ggplot(data=data, aes(x = interaction(read, genome_ord, sep = " | "), y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    geom_text(data=data[data$value > 3,], aes(label = round(value, 1)),
              size = 3.5, position = position_stack(vjust = 0.5),
              color = "#ffffff", fontface = "bold") +
    scale_y_continuous(breaks = c(0, 50, 100)) +
    coord_flip() +
    #facet_wrap(~sample_ord, dir = "v") +
    labs(title = paste0("Sample: ", sample_name),
         fill = "Percent Mapping") +
    scale_fill_manual(values=c("#8AC4DE", "#006FB0", "#FBA683", "#D40822", "#999999")) +
    theme_classic() +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
          plot.title = element_text(size = 24, face = "bold", margin = margin(t=0, r=0, b=12, l=0, unit = "pt")),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.ticks.y = element_blank(),
          #legend.background = element_rect(fill="#D9D9D9",
          #                                 size=0.5, linetype="solid", 
          #                                 colour ="#3B3B3B"),
          legend.position = "top",
          legend.justification = "left",
          #legend.title = element_blank(),
          strip.background = element_rect(fill = "#ffffff", color = "#000000"),
          strip.text = element_text(size = 8, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}

# ---------- END SCRIPT ---------- #


### Basic multispecies data analysis
# The goal of this file is conduct all of the basic analysis that underlies the multispecies data analysis.

### Get and format xenoGI information

#Read in the data on the locus family for each gene from xenoGI


locus_family_column_names <-
  c(
    "gene_name",
    "family_origin",
    "gene_history",
    "locus_island",
    "initial_family",
    "origin_family",
    "locus_family",
    "loc_fam_MRCA",
    "gene_description"
  )

paths <- list.files("xenoGI_outputs/analysis", pattern = "genes", full.names = TRUE)

all_locus_fam_genes <-
  read_tsv(
    paths,
    skip = 2,
    col_names = locus_family_column_names,
    show_col_types = FALSE,
    comment = "#"
  ) |>
  filter(!str_detect(gene_name, "contig")) %>%
  extract(
    col = gene_name,
    into = c("xenoGI_num", "species", "abrv_ID"),
    regex = "^([0-9]+)_([A-Za-z0-9_]+)-(.+)$",
    remove = TRUE,
    convert = TRUE
  )

rm(locus_family_column_names)
rm(paths)


# Importing geneinfo.txt as geneinfo

# The file `geneinfo.txt` contains the most comprehensive set of gene identifiers for all genes in this study. 
#I will then add the locus families to this data set.

geneinfo <-
  read_tsv(
    "xenoGI_outputs/geneInfo.txt",
    comment = "#",
    show_col_types = FALSE,
    col_names = c(
      "xenoGI_num",
      "species",
      "gene",
      "locus_tag",
      "WP_num",
      "gene_description",
      "genbank_file",
      "start",
      "stop",
      "strand"
    )
  ) |>
  extract(species , "species", "^[0-9]+_(.+)-", convert = TRUE) 


geneinfo_with_locus_fam <- geneinfo |>
  left_join(all_locus_fam_genes) |>
  dplyr::select(locus_family, loc_fam_MRCA, species, everything()) |>
  arrange(locus_family) 



## Cleaning up the geneinfo tibble to remove duplicates

#The entries in `geneinfo` *should* be unique, but they aren't. 
#In the table below, we can see a set of 6 genes that have either alternative start codons or alternative stop codons. 
#This is going to cause problems later, where we should have unique identifiers for each gene.

# geneinfo_with_locus_fam |>
#   group_by(locus_tag, species) |>
#   filter(n() > 1) |>
#   arrange(locus_tag) |>
#   dplyr::select(locus_tag, start, stop, everything())


#In addition, there are some cases where there are paralogues in genomes that have the same locus_family. 
#This is also going to cause problems with the analysis, because I want to use locus_family to identify orthologues.


# geneinfo_with_locus_fam |>
#   group_by(locus_family, species)  |>
# filter(n() > 1) |>
#   arrange(locus_tag) |>
#   dplyr::select(locus_tag, start, stop, everything())


#In both cases, I'm going to only take the first member of each of these duplicated families.


geneinfo_with_locus_fam_no_dups <- geneinfo_with_locus_fam |>
  group_by(locus_tag, species) |>
  slice_head(n = 1) |> 
  ungroup() |> 
  group_by(locus_family, species) |> 
  slice_head(n = 1)



# Importing and formatting counts

## Import counts

#These are the counts that I created via BWA and and htseq-count. 
#First I will get the list of all of the files and then remove the ones that didn't pass intial QC.


paths <- list.files("RNAseq_data/counts", full.names = TRUE)

#Don't read in files that didn't pass quality control.
#See the file `quality control.Rmd` for details

bad_quality <- c("E_JH04_3", "B_JH05_0", "E_JH05_0","E_JH07_0")

path_positions_to_keep <- bad_quality |> 
  map(\(x) str_detect(paths, x, negate = TRUE)) |> 
  as.data.frame(col.names = paste0("X", 1:4)) |> 
 rowSums() == 4 

paths_to_use <- paths[path_positions_to_keep]
rm(paths)
  

#The command below reads in all of the files and creates one giant (tidy) tibble, with each count from a single replicate as an experiment.
#Lots of data about the experiment was encoded in the name of the file, so the file name is saved to the new column `id`

read_counts <-
  read_tsv(
    paths_to_use,
    col_names = c("locus_tag", "count"),
    id = "file",
    show_col_types = FALSE
  )

#Now I can start converting all that information in the file name into useful information



counts_with_sample_info <- read_counts |> 
  mutate(replicate = case_when(
    str_detect(file, "counts/A") ~ "A",
    str_detect(file, "counts/B") ~ "B",
    str_detect(file, "counts/D") ~ "D",
    str_detect(file, "counts/E") ~ "E",
    TRUE ~ NA_character_
  ),
  species = case_when(
    str_detect(file, "JH03") ~ "S_typhimurium_14028s",
    str_detect(file, "JH04") ~ "C_rodentium_ICC168",
    str_detect(file, "JH05") ~ "K_pneumoniae_KPNIH1",
    str_detect(file, "JH06") ~ "E_cloacae_ATCC13047",
    str_detect(file, "JH07") ~ "S_marcescens_Db11",
    str_detect(file, "JH08") ~ "E_coli_K12",
    str_detect(file, "JH09") ~ "E_coli_K12",
    str_detect(file, "JH20") ~ "S_typhimurium_14028s",
    TRUE ~ NA_character_
  ),
  genotype = case_when(
    str_detect(file, "JH03") ~ "wt",
    str_detect(file, "JH04") ~ "wt",
    str_detect(file, "JH05") ~ "wt",
    str_detect(file, "JH06") ~ "wt",
    str_detect(file, "JH07") ~ "wt",
    str_detect(file, "JH08") ~ "wt",
    str_detect(file, "JH09") ~ "rpoS",
    str_detect(file, "JH20") ~ "rpoS",
    TRUE ~ NA_character_
  ),
  temp = case_when(
    str_detect(file, "_0_") ~ "37C",
    str_detect(file, "_3_") ~ "15C",
    TRUE ~ NA_character_
  )) 

counts_by_xenoGI_num <- counts_with_sample_info |> 
  left_join(geneinfo_with_locus_fam_no_dups, by = c("locus_tag", "species")) |> 
  filter(!is.na(xenoGI_num))


### An aside

# Am I justified in removing the genes that don't have a xenoGI number?
#   (There are 9210 *rows* that don't have xenoGI info.) What is going on
# here?
# 
# 
# library(skimr)
# 
# counts_with_sample_info |> 
#   left_join(geneinfo_with_locus_fam_no_dups, by = c("locus_tag", "species")) |> 
#   skim()
# 
# 
# Turns out that these entries are all either RNA genes (rRNA, tRNA,
# etc...) or there is no entries associated with them. They are not CDS
# genes. So I can discard them for this analysis.


# #Get the list of all the locus_tags without a xenoGI_num
# 
# locus_tag_without_xenoGI_num <- counts_by_xenoGI_num |>
#   filter(is.na(xenoGI_num)) |>
#   dplyr::select(locus_tag) |>
#   unique()
# 
# #Get all the gff files
# 
# gff_paths <-
#   list.files(
#     path = "RNAseq_data/annotation files",
#     pattern = ".gff$",
#     recursive = TRUE,
#     full.names = TRUE
#   )
# 
# combined_gffs <-
#   read_tsv(
#     gff_paths,
#     comment = "#",
#     col_names = c(
#       "genome",
#       "source",
#       "kind_of_element",
#       "start",
#       "end",
#       "X6",
#       "strand",
#       "X8",
#       "the_good_stuff"
#     )
#   )
# 
# # Extract the locus_tags and the gene_biotypes from the gffs.
# 
# combined_gffs_locus_tags <- combined_gffs |> 
#   mutate(locus_tag = str_match(the_good_stuff, ";locus_tag=(.+?);")[,2],
#          gene_biotype = str_match(the_good_stuff, ";gene_biotype=(.+?);")[,2]) |> 
#   filter(kind_of_element == "gene")
# 
# # Look at the gene_biotypes associated with the locus_tags with no xenoGI info
# locus_tag_without_xenoGI_num |> 
#   left_join(combined_gffs_locus_tags, by = "locus_tag") |> 
#   dplyr::select(gene_biotype) |> 
#   unique()
# 
# #What about the locus_tags with NA gene_biotype?
# 
# locus_tag_without_xenoGI_num |> 
#   left_join(combined_gffs_locus_tags, by = "locus_tag") |> 
#   filter(is.na(gene_biotype)) |>
#   group_by(is.na(the_good_stuff)) |> 
#   tally()
# 
# 
# 
# 
# If a gene with no xenoGI_num isn't annotated as an RNA, it has no
#    annotation at all. I think I'm safe dropping these genes from the
# analysis.


# Create a table with sample information that has the information on each sample. 
# The `sample_name` column created here will match the column names for each sample created when I `pivot_wider` below. 
# It will need to be moved to rownames in order to work with DESeq2


sample_table <- counts_by_xenoGI_num |>
  dplyr::select(replicate, species, genotype, temp) |>
  distinct() |>
  mutate(
    condition = paste(species, genotype, temp, sep = "_"),
    sample_name = paste(replicate, species, genotype, temp, sep = "_"),
    species = factor(species),
    genotype = factor(genotype),
    temp = factor(temp),
    condition = factor(condition)
  ) |>
  column_to_rownames(var = "sample_name")


# I think that `counts_by_xenoGI_with_NA` doesn't get used anywhere, but I'm not deleting it just yet in case I'm wrong. 
# counts_by_xenoGI_with_NA <- counts_by_xenoGI_num |> 
#      pivot_wider(id_cols = locus_family, names_from =c(replicate, species, genotype, temp), values_from = count) |> 
#      filter(is.na(count))
   
   
counts_by_xenoGI_num_wide <- counts_by_xenoGI_num |> 
     pivot_wider(id_cols = locus_family, names_from =c(replicate, species, genotype, temp), values_from = count, values_fill = 0) 


   
# Find significant orthologs for each species.

sig_alpha <- 0.01 # Our significance threshold for this analysis

   
# General function to find signficiant genes
  
   
run_DESeq_single_contrast <- function(species_to_use, genotype_to_use, temp_to_use, variable_to_test){
     
     args <- list(species_to_use, genotype_to_use, temp_to_use)
     single_level <- map(args, length) == 1
     
     output_object <- paste(paste(args[single_level], collapse = "_"), "significance_of", variable_to_test, sep="_")
     output_object
     output_file <- paste0("outputs/sig_tables/",output_object,".tsv")
     design_formula <- formula(paste0("~ ", variable_to_test))
     
     
     samples_to_test <- sample_table |>
       filter(genotype %in% genotype_to_use,
              species %in% species_to_use,
              temp %in% temp_to_use)
     
     counts_to_test <- counts_by_xenoGI_num_wide |>
       dplyr::select(c("locus_family", rownames(samples_to_test))) |>
       # rowwise() |>
       # filter(sum(c_across(rownames(samples_to_test))) != 0) |>
       column_to_rownames("locus_family")
     
     
     dds <- DESeqDataSetFromMatrix(countData = counts_to_test,
                                   colData = samples_to_test,
                                   design = design_formula)
     keep <- rowSums(counts(dds)) >= 10
     dds <- dds[keep,]
     
     dds$temp <- relevel(dds$temp, ref = "37C")
     dds$genotype <- relevel(dds$genotype, ref = "wt")
     dds <- DESeq(dds)
     
    
     
     output_object <- dds |>
       results() |>
       as.data.frame() |>
       as_tibble(rownames = "locus_family") |>
       mutate(locus_family = as.numeric(locus_family)) |>
       #Defines which genes are significant for the rest of the analysis. Change this if you want a different adjusted p-value of lfc
       mutate(significant = padj < sig_alpha &
                abs(log2FoldChange) > 1) |>
       mutate(significant = case_when(is.na(padj) ~ NA,
                                      TRUE ~ significant)) |>
       left_join(filter(geneinfo_with_locus_fam_no_dups, species == species_to_use)) |>
       dplyr::select(
         locus_tag,
         gene,
         baseMean,
         log2FoldChange,
         padj,
         significant,
         species,
         locus_family,
         WP_num,
         gene_description,
         start,
         stop,
         strand
       )
     
     write_tsv(output_object, output_file)
     
     output_object
     
   }

#Calculate the significant genes for each species in the wt case
   
species <- sample_table |> 
     pull(species) |> 
     unique()

dir.create("outputs/sig_tables/")


all_wt_significance <- species |>
  as.character() |> 
  set_names() |>
  purrr::map(
    \(x) run_DESeq_single_contrast(
      species_to_use = x,
      genotype_to_use = "wt",
      temp_to_use = c("37C", "15C"),
      variable_to_test = "temp"
    )
  )

#I can use species labels both as factor levels and with the labeller function for plots.

species_labels <- c(E_coli_K12 = "E. coli", S_typhimurium_14028s = "S. enterica", C_rodentium_ICC168 = "C. rodentium", E_cloacae_ATCC13047 = "E. cloacae", K_pneumoniae_KPNIH1 = "K. pneumoniae", S_marcescens_Db11 ="S. marcescens")



#The flat version (rather than list) will be useful for many downstream applications. 

all_wt_sig_flat <-bind_rows(all_wt_significance) |> 
  mutate(species = fct_relevel(species, names(species_labels)))


   
#Calculate the significant genes for each species in the ∆rpoS case
   

   mutant_species <- sample_table |> 
     filter(genotype == "rpoS") |> 
     pull(species) |> 
     unique()
   
   all_rpoS_significance <- mutant_species |>
     as.character() |> 
     set_names() |>
     map(\(x) run_DESeq_single_contrast(species_to_use = x, genotype_to_use = "rpoS", temp_to_use = c("37C", "15C"),  variable_to_test = "temp"))
   
   wt_vs_rpoS_37_significance <- mutant_species |>
     as.character() |> 
     set_names() |>
     map(\(x) run_DESeq_single_contrast(species_to_use = x, genotype_to_use = c("wt", "rpoS"), temp_to_use = "37C",  variable_to_test = "genotype"))
   
   wt_vs_rpoS_15_significance <- mutant_species |>
     as.character() |> 
     set_names() |>
     map(\(x) run_DESeq_single_contrast(species_to_use = x, genotype_to_use = c("wt", "rpoS"), temp_to_use = "15C",  variable_to_test = "genotype"))

   
all_dds <-
  DESeqDataSetFromMatrix(
    countData = column_to_rownames(counts_by_xenoGI_num_wide, "locus_family"),
    colData = sample_table,
    design = ~ condition
  ) |>
  DESeq() 

normed_counts <- counts(all_dds, normalized = TRUE) |>
  as.data.frame() |>
  rownames_to_column(var = "locus_family") |>
  pivot_longer(
    cols = !contains("locus_family"),
    names_to = "sample",
    values_to = "normed_count"
  )


sample_table_tidy <- sample_table |> 
  rownames_to_column(var = "sample")

locus_fam_ecoli_names <- geneinfo_with_locus_fam_no_dups |> 
  filter(species == "E_coli_K12") |> 
  ungroup() |> 
  dplyr::select(locus_family, E_coli_K12_gene = gene) 

locus_tag <- geneinfo_with_locus_fam_no_dups |> 
  dplyr::select(locus_family, species, locus_tag)



merged_data <-
  left_join(normed_counts, sample_table_tidy, by = "sample") |>
  mutate(locus_family = as.numeric(locus_family)) |> 
  left_join(locus_fam_ecoli_names, by = join_by(locus_family)) |> 
  left_join(locus_tag, by = join_by(species, locus_family)) |> 
  mutate(normed_count = case_when(
    is.na(locus_tag) ~ NA_real_,
    TRUE ~ normed_count),
    species = fct_relevel(
      species,
      names(species_labels)
    )
  )








plot_gene <- function(count_data = merged_data, ecoli_gene_name = NA, locus_family_num = NA, extra_title_text = "") {
  
   

    if(is.na(ecoli_gene_name) && is.na(locus_family_num)){
      stop("Either an E. coli gene name or a locus family number are required")
      
    }

      
    

  
   if (!is.na(ecoli_gene_name)) {
    count_data %>%
      filter(E_coli_K12_gene == ecoli_gene_name) %>%
      ggplot(aes(x = temp, y = normed_count, color = genotype)) +
      geom_jitter(width = .1,
                  alpha = .1,
                  size = 5) +
      stat_summary(
        data = subset(
          count_data,
          genotype == "wt" &
            E_coli_K12_gene == ecoli_gene_name
        ),
        fun.data = mean_se,
        geom = "errorbar",
        width = .1
      ) +
      stat_summary(
        data = subset(
          count_data,
          genotype == "wt" & E_coli_K12_gene == ecoli_gene_name
        ),
        fun = mean
      ) +
      facet_grid(~ species,
                 labeller = labeller(species = species_labels)) +
      scale_color_brewer(palette = "Set1") +
      theme(text = element_text(size = 10),
            strip.text.x = element_text(face = "italic")) +
      labs(title = paste(ecoli_gene_name,extra_title_text), 
             x = "Temperature",
           y = "Normalized count") +
       scale_x_discrete(labels = c("15°C", "37°C")) +
      NULL
  }
  
  else{
    count_data %>%
      filter(locus_family_num == locus_family) %>%
      ggplot(aes(x = temp, y = normed_count, color = genotype)) +
      geom_jitter(width = .1,
                  alpha = .1,
                  size = 5) +
      stat_summary(
        data = subset(count_data, genotype == "wt" &
                        locus_family_num == locus_family),
        fun.data = mean_se,
        geom = "errorbar",
        width = .1
      ) +
      stat_summary(
        data = subset(count_data, genotype == "wt" &
                        locus_family_num == locus_family),
        fun = mean
      ) +
      facet_grid( ~ species,
                  labeller = labeller(species = species_labels)) +
      scale_color_brewer(palette = "Set1") +
      theme(text = element_text(size = 10),
            strip.text.x = element_text(face = "italic")) +
      labs(title = paste("locus family",locus_family_num,extra_title_text),
           x = "Temperature",
           y = "Normalized count") +
      scale_x_discrete(labels = c("15°C", "37°C")) +
      NULL
  }
}

plot_wt_gene <- function(count_data = merged_data, ecoli_gene_name = NA, locus_family_num = NA, extra_title_text = "") {
  
  
  if(is.na(ecoli_gene_name) && is.na(locus_family_num)){
    stop("Either an E. coli gene name or a locus family number are required")
    
  }
  

  if (!is.na(ecoli_gene_name)) {
    count_data %>%
      filter(E_coli_K12_gene == ecoli_gene_name,
             genotype == "wt") %>%
      ggplot(aes(x = temp, y = normed_count)) +
      geom_jitter(width = .1,
                  alpha = .1) +
      stat_summary(
        data = subset(
          count_data,
          genotype == "wt" &
            E_coli_K12_gene == ecoli_gene_name
        ),
        fun.data = mean_se,
        geom = "errorbar",
        width = .1
      ) +
      stat_summary(
        data = subset(
          count_data,
          genotype == "wt" & E_coli_K12_gene == ecoli_gene_name
        ),
        fun = mean
      ) +
      facet_grid(~ species,
                 labeller = labeller(species = species_labels)) +
      scale_color_brewer(palette = "Set1") +
      theme(text = element_text(size = 10),
            strip.text.x = element_text(face = "italic")) +
      labs(title = paste(ecoli_gene_name,extra_title_text),
           x = "Temperature",
           y = "Normalized count") +
      scale_x_discrete(labels = c("15°C", "37°C")) +
      NULL
  }
  
  else{
    count_data %>%
      filter(locus_family_num == locus_family,
             genotype == "wt") %>%
      ggplot(aes(x = temp, y = normed_count, color = genotype)) +
      geom_jitter(width = .1,
                  alpha = .1,
                  size = 1.25) +
      stat_summary(
        data = subset(count_data, genotype == "wt" &
                        locus_family_num == locus_family),
        fun.data = mean_se,
        geom = "errorbar",
        width = .1
      ) +
      stat_summary(
        data = subset(count_data, genotype == "wt" &
                        locus_family_num == locus_family),
        fun = mean,
        size = 1.5
      ) +
      facet_grid( ~ species,
                  labeller = labeller(species = species_labels)) +
      scale_color_brewer(palette = "Set1") +
      theme(text = element_text(size = 10),
            strip.text.x = element_text(face = "italic")) +
      labs(title = paste("locus family",locus_family_num,extra_title_text),
           x = "Temperature",
           y = "Normalized count") +
      scale_x_discrete(labels = c("15°C", "37°C")) +
      NULL
  }
}


### The below function is used in several .Rmd files
format_species_names <- function(name){
  case_when(name == "E_coli_K12" ~ "E. coli",
            name == "S_typhimurium_14028s" ~ "S. enterica",
            name == "C_rodentium_ICC168" ~ "C. rodentium",
            name == "E_cloacae_ATCC13047" ~ "E. cloacae",
            name == "K_pneumoniae_KPNIH1" ~ "K. pneumoniae",
            name == "S_marcescens_Db11" ~ "S. marcescens",
            TRUE ~ name)
}


save(list = ls(), file = "basic_analysis.RData")

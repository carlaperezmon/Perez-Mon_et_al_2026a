### SET WORKING DIRECTORY, LOAD PACKAGES, METADATA AND DATA --------

path_work="M:/postdoc_2021/sequencing_results/data_analyses/2026-01_trnl_who"

setwd(path_work)

load("trnl_analyses_hist_seqs2026.RData")


# load packages, working paths and metadata ------
ipak <- function(pkg){
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  } #install BiocManager if needed 
  lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        install.packages(pkg) #install packages from CRAN
        library(pkg, character.only = TRUE) 
      }, error = function(e) {
        BiocManager::install(pkg, force = TRUE) # if package is not found in CRAN it installs it from Biocmanager. Force=TRUE install the package even if incompatibility of version of R. This needs to be treated with Caution though, as it is not warranty that the package will run smoothly if created in a more recent R version than the one used by the user
        library(pkg, character.only = TRUE)  
      })
    }
  })
  sapply(pkg, require, character.only = TRUE) #tells you if packages have been loaded
}


packages <- c("R.utils","Biostrings","phyloseq","Biostrings","ggplot2","openxlsx","ape","reshape2","qiime2R",
              "stringr","dplyr","Polychrome","Rmisc","rgl","EcolUtils","vegan","plot3D","gridExtra","scatterplot3d",
              "indicspecies", "DESeq2","gplots","Rmisc","rworldmap","countrycode",'rWCVP','BIEN',"kewr","viridis",
              "ggpubr","gridExtra","rWCVPdata","maps","rnaturalearth","patchwork")

ipak(packages)


#create subdirectories

subDir <- "tables_general_stats"
dir.create(file.path(path_work, subDir), showWarnings = TRUE)

subDir1 <- "plots_counts_taxa"
dir.create(file.path(path_work, subDir1), showWarnings = TRUE)

subDir2 <- "pow_ind_asvs_maps"
dir.create(file.path(path_work, subDir2), showWarnings = TRUE)

subDir3 <- "pow_aggregated_maps"
dir.create(file.path(path_work, subDir3), showWarnings = TRUE)


# import metadata and abundance tables #

metadata=read.xlsx("metadata_who_trnl.xlsx", sheet= "Sheet1")
colnames(metadata)[colnames(metadata)=="sample-id"] <-"sample.id"
metadata$lot_id2_rep=paste(metadata$lot_id2,metadata$replicate,sep = "_")

rownames(metadata)=metadata$sample.id

shortcodes=read.xlsx("who_shortcodes.xlsx", sheet= "Sheet1")
metadata$lot_id3=shortcodes$Short.code[match(metadata$lot_id,shortcodes$PACS.code)]
metadata$lot_id3_rep=paste(metadata$lot_id3,metadata$replicate,sep = "_")
metadata=metadata[order(metadata$lot_id3_rep),]

table_trnl=read_qza('140_120/table-dada2_trnl.qza')
table_trnl=table_trnl$data

seqs_trnl=read_qza('140_120/rep-seqs-dada2_trnl.qza')
seqs_trnl=seqs_trnl$data


#filter control from abundance table and filter ASVs for overall abundance #

number_filtering=9

names_neg_control=na.omit(metadata$sample.id[metadata$material=="control"]) #remove also the reads in the sept 2024 controls from the other A5J3 runs
  
neg_control=table_trnl[,colnames(table_trnl) %in% names_neg_control]
neg_control=neg_control[rowSums(neg_control) > number_filtering,]

table_trnl_filt=table_trnl[!(rownames(table_trnl) %in% rownames(neg_control)),]
table_trnl_filt=table_trnl[,!(colnames(table_trnl) %in% names_neg_control)]
table_trnl_filt=table_trnl_filt[rowSums(table_trnl_filt) > 9,]
  
  
### LOAD trnl-blast SPECIES ANNOTATIONS -----

all_blast=read.table("140_120/blast_trnl_pe0.9_qc0.9_NCBI_export/blast6.tsv",sep = "\t")
colnames(all_blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
all_blast=all_blast[all_blast$pident > 0,]

taxonomy_identities=read.csv("K:/Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/databases/trnl_ncbi/taxonomy_r/taxa_tableforqiime.txt",sep='\t', header=FALSE)
taxonomy_identities[c("domain", "phylum", "class", "order", "family", "genus","species")] <- str_split_fixed(taxonomy_identities$V2,";",7) 

all_blast$species=taxonomy_identities$species[match(all_blast$sseqid,taxonomy_identities$V1)]
all_blast$species=gsub("s__","", all_blast$species)

all_blast$genus=taxonomy_identities$genus[match(all_blast$sseqid,taxonomy_identities$V1)]
all_blast$genus=gsub("g__","", all_blast$genus)


### FILTER OUT ZEA (CORN) READS ------
# REMOVE ASVs that  annotated to zea at percentage identity > 0.9 

blast_zea=all_blast[all_blast$genus=='Zea',] 

blast_zea2 <- blast_zea  %>% 
  dplyr::select(qseqid,pident,species) %>%
  dplyr::group_by(qseqid,species) %>%
  dplyr::summarise(mean_pident=round(mean(pident),2))%>%
  arrange(desc(mean_pident)) %>%
  mutate(qseqid_pident=paste(qseqid,' (',mean_pident,')'))

all_zea_sps_by_featureid <- blast_zea2 %>%
  group_by(qseqid) %>%
  summarise(values_collapsed = paste(species, collapse = ","), .groups = "drop")

table_trnl_filt_zea=table_trnl_filt[rownames(table_trnl_filt) %in% blast_zea2$qseqid,]
table_trnl_filt_zea=melt(table_trnl_filt_zea)

table_trnl_filt_zea$qseqid_pident=blast_zea2$qseqid_pident[match(table_trnl_filt_zea$Var1, blast_zea2$qseqid)]
table_trnl_filt_zea$mean_pident=blast_zea2$mean_pident[match(table_trnl_filt_zea$Var1, blast_zea2$qseqid)]
table_trnl_filt_zea$lot_id3_rep=metadata$lot_id3_rep[match(table_trnl_filt_zea$Var2,metadata$sample.id)]

#generate abundance table to evaluate how much Zea interference
png(paste(path_work,'/', subDir1, '/asvs_zea_heatmap.png', sep=""), width = 1000, height = 480)
ggplot(na.omit(table_trnl_filt_zea), aes(lot_id3_rep, factor(qseqid_pident, levels=unique(rev(blast_zea2$qseqid_pident))), fill=log10(value))) + 
  geom_tile() +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "red", 
    midpoint = 0
  ) + 
  ylab('Log10 read counts') + 
  xlab('')+
  theme_classic2()
dev.off()

table_trnl_filt_zea2=table_trnl_filt_zea %>% 
  group_by(lot_id3_rep) %>%
  dplyr::summarise(values_aggregated = sum(value, na.rm = TRUE), .groups = "drop")

  
ggplot(table_trnl_filt_zea2, aes(lot_id3_rep,values_aggregated)) + 
  geom_col() + xlab('') + ylab('read counts annotated to Zea ( > 0.90 pident)') + 
  theme_classic2()


# retention of non-zea ASVs

all_blastnozea=all_blast[!(all_blast$qseqid %in% blast_zea$qseqid),] 
length(unique(all_blastnozea$qseqid))

threshold=98 


### SELECTION OF BEST MATCHES PER ASVS FOR FURTHER ANALYSES AND WITH SPECIES ASSIGNMENTS > 98% PERCENTAGE IDENTITY-------

all_blastnozea=all_blastnozea[all_blastnozea$pident >= threshold,] 
all_blast_ev=all_blastnozea %>%
  dplyr::count(qseqid, sort = TRUE)

feature_ids=unique(all_blastnozea$qseqid)
best_hits_selection=list()
  
for (i in 1:length(feature_ids)) {
  print(feature_ids[i])
  
  max_hits <- all_blastnozea %>%
    filter(qseqid==feature_ids[i]) %>%
    filter(pident==max(pident))
  
  best_hits_selection [[i]] <- max_hits
  
}

  
all_best_hits_selection=dplyr::bind_rows(best_hits_selection)
  
all_best_hits_selection2 <- all_best_hits_selection  %>% 
  dplyr::select(qseqid,pident,species) %>%
  dplyr::group_by(qseqid,species) %>%
  dplyr::summarise(mean_pident=mean(pident))

all_blast_ev=all_best_hits_selection2 %>%
  select(qseqid,species)%>%
  na.omit()%>%
  unique()%>%
  dplyr::count(qseqid, sort = TRUE)

all_blast_ev_10=all_blast_ev[all_blast_ev$n < 6,]

all_best_hits_selection2=all_best_hits_selection2[all_best_hits_selection2$qseqid %in% all_blast_ev_10$qseqid,]  

# merge blast information with count information
  
trnl_blast_map=merge(table_trnl_filt,all_best_hits_selection2,by.x="row.names",by.y="qseqid")


### CREATION OF TABLE FOR SELECTED ASVS PER SPECIES AND THEIR ABUNDANCE ACROSS SAMPLES (FigS19) -------

table_trnl_filt_subset=as.data.frame(t(table_trnl_filt[rownames(table_trnl_filt) %in% trnl_blast_map$Row.names,]))  
table_trnl_filt_subset$lot_id3=metadata$lot_id3[match(rownames(table_trnl_filt_subset),metadata$sample.id)]  

table_trnl_filt_subset=melt(table_trnl_filt_subset)
table_trnl_filt_subset$country=metadata$collection_country[match(table_trnl_filt_subset$lot_id3,metadata$lot_id3)]

sps_by_featureid <- all_best_hits_selection2 %>%
  group_by(qseqid) %>%
  dplyr::summarise(values_collapsed = paste(species, collapse = ","), .groups = "drop")

sps_by_featureid$short_name=paste('trnl-ASV',seq(1,length(sps_by_featureid$qseqid)),sep='-')

table_trnl_filt_subset$short_name=sps_by_featureid$short_name[match(table_trnl_filt_subset$variable,sps_by_featureid$qseqid)]

asvs_annotated_abundance_plot=ggplot(table_trnl_filt_subset,aes(x=short_name,y=value,fill = interaction(lot_id3, country))) + 
  labs(y='Total reads',x='') +
  geom_col(width=.5, position = "stack")+
  scale_fill_viridis_d(option = "magma",name='Sample & Country')+
  theme_classic()+
  coord_flip()

write.xlsx(sps_by_featureid[,!colnames(sps_by_featureid) %in% 'qseqid'], 'asv_trnl_select_species.xlsx')
ggsave('asv_trnl_select_abunplot_species.png',asvs_annotated_abundance_plot)


### LOAD POWO DATA FOR MAP REPRESENTATIONS (FigS20) ------
names <- rWCVPdata::wcvp_names
distributions <- rWCVPdata::wcvp_distributions

world_map <- ne_countries(scale = "medium", returnclass = "sf")

base_world <- ggplot() +
  geom_sf(data = world_map, fill = "grey90", color = "grey80")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="left",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())


### REPRESENTATION POW for each asv and species separately ####  

#arrange tables
trnl_blast_map2=trnl_blast_map
colnames(trnl_blast_map2)[colnames(trnl_blast_map2)=='Row.names'] <-'Feature.ID'

add_dipsacus=c('ASV-18S-dipsacus',rep(as.numeric(1),35),'Dipsacus asper',100)   #dipsacus asper was detected with 18S-marker, blast approach, in multiple samples

trnl_blast_map2=rbind(trnl_blast_map2,add_dipsacus)

asv_id=unique(trnl_blast_map2$Feature.ID)
 
#download POWO distributions   
all_distributions=data.frame()
    
for (j in (1:length(asv_id))) {
      print(asv_id[j])
      
      asv_id2=asv_id[j]
      
      species_to_query <- unique(trnl_blast_map2$species[trnl_blast_map2$Feature.ID==asv_id[j]])
      
      plot_analysis <- list()
      
      for (k in (1:length(species_to_query))) {
        print(species_to_query[k])
        
        powo_results=purrr::map_dfr(search_powo(species_to_query[k])$results, unlist)
        
        if(nrow(powo_results)==0) {
          print(species_to_query[k])
          
          print ("not found")} else {
            
            species_to_query2=powo_results$name[powo_results$accepted==TRUE]
            
            try({
              
              distribution <- wcvp_distribution(species_to_query2[1], taxon_rank = "species", wcvp_names = names, 
                                                wcvp_distributions = distributions, extinct = FALSE,
                                                location_doubtful = FALSE)
              
              plot_analysis[[paste0(species_to_query[k],'_plot')]] <-wcvp_distribution_map(distribution) + ggtitle(paste(species_to_query[k], "=",species_to_query2))
              
              distribution$species=rep(species_to_query[k], nrow(distribution))
              distribution$asv_id=rep(asv_id[j], nrow(distribution))
              distribution$country_seq=asv_id2
              
              all_distributions=rbind(all_distributions,distribution)
              
            },silent = TRUE) }
        
      }
      
      if (length(plot_analysis) ==0) { print ("nothing to plot")} else {
        if (length(plot_analysis) > 4) {
          range_loop=seq(1,length(plot_analysis),by=4) 
          
          pdf(file = paste(path_work,'/', subDir2, '/',asv_id2,"_noZea_species.pdf",sep=""), paper="a4r",width=10, height=7)
          for(o in 1:length(range_loop)) {
            l=range_loop[o]+3
            
            if (l < length(plot_analysis)) {
              plot_analysis2=plot_analysis[range_loop[o]:l]
              
            } else { 
              plot_analysis2=plot_analysis[range_loop[o]:length(plot_analysis)]
            }
            do.call('grid.arrange',c(plot_analysis2, ncol= 2,nrow=2))
          }
          
          dev.off()      
          
          
          
        } else { 
          pdf(file = paste(path_work,'/', subDir2, '/',asv_id2,"_noZea_species.pdf",sep=""), paper="a4r",width=10, height=7)
          do.call('grid.arrange',c(plot_analysis, ncol= 2,nrow=2))  
          dev.off()
        }      
      }
    
    assign("all_distributions_asvs",all_distributions)
    
}
  


### REPRESENTATION POW for asv separately (aggregation of all possible sps per ASV) ####  

all_distributions_asvs=all_distributions
    
asv_id=unique(all_distributions$asv_id)
    
for (j in (1:length(asv_id))) {
      print(asv_id[j])
      
      asv_id2=asv_id[j]
      
      possible_sps=paste(unique(all_distributions$species[all_distributions$asv_id==asv_id[j]]), collapse = ",")
      
      all_distributions_subset=all_distributions[all_distributions$asv_id==asv_id[j],]
      all_distributions_subset$occurrence=rep(1,nrow(all_distributions_subset))
      all_distributions_subset$rel_occurrence=all_distributions_subset$occurrence/nrow(all_distributions_subset)
      
      all_distributions_subset <- all_distributions_subset %>%
        group_by(LEVEL3_COD) %>%
        dplyr::summarise(sum_occurrence = sum(occurrence), sum_rel_occurrence=sum(rel_occurrence)) %>%
        arrange(sum_occurrence)
      
      
      color_palette <- viridis(10)
      
      # Map numeric values to colors
      all_distributions_subset$numeric_colors <- color_palette[cut(all_distributions_subset$sum_occurrence, 10)]
      
      bot_regions=unique(all_distributions_subset$LEVEL3_COD)
      
      p <- base_world
      
      
      for (l in (1:length(bot_regions))) {
        print(unique(bot_regions[l]))
        
        all_distributions_subset2=all_distributions_subset[all_distributions_subset$LEVEL3_COD==bot_regions[l],]
        
        p=p + geom_sf(data = all_distributions_subset2$geometry[1], fill=all_distributions_subset2$numeric_colors)
        
      }
      
      p2=p+ggtitle(paste("asv_id:",asv_id2)) + labs(subtitle=str_wrap(paste('Aggregated distribution of species:', possible_sps), 100)) + theme(plot.subtitle=element_text(face = "italic"))
      
      color_bar=cbind(all_distributions_subset$sum_occurrence,all_distributions_subset$numeric_colors)
      color_bar=data.frame(unique(color_bar))
      #color_bar=(color_bar[!duplicated(color_bar),])
      color_bar$y=rep(1,nrow(color_bar))
      levels=as.character(color_bar$X1)
      
      color_bar_plot=ggplot(color_bar, aes(x=factor(X1, levels=levels), y=y)) + 
        geom_bar(stat = "identity",  fill=color_bar$X2) + 
        scale_y_continuous(expand = c(0,0)) + 
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.x=element_blank())
      
      
      all_plot=ggpubr::ggarrange(p2, color_bar_plot, heights = c(10, 0.7),widths = c(10, 0.7),
                                 ncol = 1, nrow = 2, align = "none")
      
      
      ggsave(paste(path_work,"/",subDir2,'/',asv_id2,"_noZea_sps_aggregated.pdf",sep=""),all_plot,device = cairo_pdf,
             width = 30, height = 15, units = "cm")
    }
    


### Representation of all species distributions for all ASVs aggregated ####  
#split between including natives and introduced, or only native#

all_distributions_asvs$asv_id_short=sps_by_featureid$short_name[match(all_distributions_asvs$asv_id,sps_by_featureid$qseqid)]

occurrence_type=c('all','native')
titles=c('including native & introduce ranges','only native ranges')

all_distributions_asvs_all=all_distributions_asvs
all_distributions_asvs_native=all_distributions_asvs[all_distributions_asvs$occurrence_type=='native',]

for (i in 1:length(occurrence_type)) {
  print(occurrence_type[i])
  
  get(paste("all_distributions_asvs",occurrence_type[i],sep ="_")) -> all_distributions_subset
  
  #all_distributions_subset=all_distributions_subset %>%
    #select(-asv_id,-asv_id_short,-country_seq)%>%
    #unique()
  
  all_distributions_subset$occurrence=rep(1,nrow(all_distributions_subset))
  all_distributions_subset$rel_occurrence=all_distributions_subset$occurrence/nrow(all_distributions_subset)

  sps=paste(unique(all_distributions_subset$species), collapse = ", ")
    
  subtitle_rep=paste("Distribution of",length(unique(all_distributions_subset$species)),"possible species, represented by",length(unique(all_distributions_subset$asv_id)), 'ASVs')
  subtitle_rep2=paste(subtitle_rep,sps,sep="\n")

 
  all_distributions_subset <- all_distributions_subset %>%
    group_by(LEVEL3_COD) %>%
    dplyr::summarise(sum_occurrence=sum(occurrence), sum_rel_occurrence=sum(rel_occurrence)) %>%
    arrange(sum_occurrence)
    
    
  color_palette <- viridis(10)
    
    # Map numeric values to colors

  all_distributions_subset$numeric_colors <- color_palette[cut(all_distributions_subset$sum_occurrence, 10)]
    
  bot_regions=unique(all_distributions_subset$LEVEL3_COD)
    
  p <- base_world
    
    for (l in (1:length(bot_regions))) {
      print(unique(bot_regions[l]))
      
      all_distributions_subset2=all_distributions_subset[all_distributions_subset$LEVEL3_COD==bot_regions[l],]
      
      p=p + geom_sf(data = all_distributions_subset2$geometry[1], fill=all_distributions_subset2$numeric_colors, color="grey20")
      
    }
    
    
  p2=p+ggtitle(subtitle_rep2) + 
      labs(caption=cat(subtitle_rep2)) + 
      theme(plot.caption = element_text(hjust = 0))+
      plot_annotation(caption = titles[i],
                        theme = theme(
                          plot.caption = element_text(size = 16, face = "bold")))
    

  color_bar=data.frame(x=seq(min(all_distributions_subset$sum_occurrence),max(all_distributions_subset$sum_occurrence),1),
                         y=0.5,
                         color=seq(min(all_distributions_subset$sum_occurrence),max(all_distributions_subset$sum_occurrence),1))
    
 
  color_bar_plot=ggplot(color_bar, aes(x=x, y = y, fill=color)) +
      geom_tile() +
      scale_y_continuous(breaks = c(0,0.5))+
      scale_fill_viridis(option = "viridis",name = "Total species\noverlapping") +  # Use the viridis color palette
      theme_void() #+  # Remove axis, grid, and background
      #theme(legend.position = "none",
            #axis.text.y=element_text(color="black", size=14),plot.margin = unit(c(5, 9.5, 5, 10), "cm"))
    
    legend <- get_legend(color_bar_plot)
    
    all_plot_agreggated=ggarrange(p2,legend, heights = c(10, 0.7),widths = c(10, 0.7),
                       ncol = 2, nrow = 1, align = "none")
    
    ggsave(paste(path_work,"/",subDir3,'/',occurrence_type[i],"_noZea_all_aggregated.pdf",sep=""),all_plot_agregatted,device = cairo_pdf,
           width = 30, height = 15, units = "cm")
   
}   
 





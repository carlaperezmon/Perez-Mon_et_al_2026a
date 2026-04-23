###### BLAST 18s GBIF differential abundant, excipients removed ######

##### set working directory and paths for data and figures ------

path_work="M:/postdoc_2021/sequencing_results/data_analyses/2024-2025_who/clustered_runs/18S"

setwd(path_work)

load("gbif_hist.RData")


###### load packages ------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


packages <- c("R.utils","Biostrings","ggplot2","openxlsx","reshape2","qiime2R",
              "stringr","dplyr","Polychrome","Rmisc","rgl","vegan","gridExtra",
              "gplots","Rmisc","rgbif","rworldmap","countrycode","ggmap","hexbin","viridis",
              "grattantheme",'RColorBrewer')

ipak(packages)


##### import WHO metadata and create representation vectors for later ########

metadata=read.xlsx("../../metadata_who_falsified_sel.xlsx", sheet= "Sheet1")
colnames(metadata)[colnames(metadata)=="sample-id"] <-"sample.id"
metadata$lot_id2_rep=paste(metadata$lot_id2,metadata$replicate,sep = "_")

rownames(metadata)=metadata$sample.id

shortcodes=read.xlsx("../../who_shortcodes.xlsx", sheet= "Sheet1")
metadata$lot_id3=shortcodes$Short.code[match(metadata$lot_id,shortcodes$PACS.code)]
metadata$lot_id3_rep=paste(metadata$lot_id3,metadata$replicate,sep = "_")
metadata=metadata[order(metadata$lot_id3_rep),]


##color vectors

order_countries=c("Niger","Chad","Cameroon","Nigeria","Ghana","Burkina Faso","Liberia")
color_countries<-brewer.pal(7,"Set3")
names(color_countries) <-order_countries



#import filtered count table

count_filt=read.table("tables_general_stats/count_filt_18S.csv",sep = ",", header=TRUE,check.names = FALSE,row.names = 1)
count_filt$total_reads=rowSums(count_filt)


# Import PR2 lineage with metadata

lineage_with_metadata=read.xlsx('K:/ogden_grp/Conservation_Science/Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/databases/pr2/pr2_version_5.0.0_merged.xlsx',sheet=1) 
lineage_with_metadata <- lineage_with_metadata[grep("Eukaryota",lineage_with_metadata$domain),] #keep only the eukaryotic sequences
lineage_with_metadata$species=gsub (".*_sp.$", NA, lineage_with_metadata$species)

# Load blast results and process -> result is a tsv table in blastn format 6: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

all_blast=read.table("blast_18S_pe0.9_qc0.9_PR2_export/blast6.tsv",sep = "\t")
colnames(all_blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

#select ASVs in counts and filter by identity threshold#

all_blast=all_blast[all_blast$qseqid %in% rownames(count_filt),]

threshold=98
all_blast=all_blast[all_blast$pident > threshold,]

all_blast=all_blast[order(all_blast$qseqid,-all_blast$pident),]


#### OPTIONAL ###### 

#best hits selection###

feature_ids=unique(all_blast$qseqid)

best_hits_selection=list()

for (i in 1:length(feature_ids)) {
  print(feature_ids[i])
  
  max_hits <- all_blast %>%
    filter(qseqid==feature_ids[i]) %>%
    filter(pident==max(pident) | bitscore==max(bitscore))
  
  best_hits_selection [[i]] <- max_hits
  
}

all_best_hits_selection=dplyr::bind_rows(best_hits_selection)

all_best_hits_selection=all_best_hits_selection[order(all_best_hits_selection$qseqid,-all_best_hits_selection$pident),]
all_best_hits_selection$species=lineage_with_metadata$species[match(all_best_hits_selection$sseqid, lineage_with_metadata$pr2_accession)]

all_blast_ev=all_best_hits_selection %>%
  select(qseqid,species)%>%
  na.omit()%>%
  unique()%>%
  dplyr::count(qseqid, sort = TRUE) #28 ASVs across all dataset 

all_blast_ev_10=all_blast_ev[all_blast_ev$n < 11,] #24 with less than 10 assigned possible species at the selected identity threshold


#link with metadata

all_blast_select=all_blast[all_blast$qseqid %in% all_blast_ev_10$qseqid,]


all_blast_pr2_metadata=merge(all_blast_select,lineage_with_metadata, by.x='sseqid',by.y='pr2_accession')%>%
  filter(genus!='Zea',
         !qseqid %in% qseqid[species=='Homo_sapiens'])%>% #remove human
  mutate(species=gsub('_',' ',species))


#Load metadata and import filtered count table

count_filt_select=count_filt[rownames(count_filt) %in% all_blast_select$qseqid,]
count_filt_select=count_filt_select[count_filt_select$total_reads>99,] #focuss on samples with 100 reads or more

all_blast_pr2_metadata_count=merge(all_blast_pr2_metadata,count_filt_select,by.x='qseqid',by.y='row.names')

species_to_search_fungi=unique(na.omit(all_blast_pr2_metadata_count$species[all_blast_pr2_metadata_count$subdivision=='Fungi']))
species_to_search_notfungi=unique(na.omit(all_blast_pr2_metadata_count$species[all_blast_pr2_metadata_count$subdivision!='Fungi']))

asv_ids_species=all_blast_pr2_metadata_count %>% 
  select(qseqid,species)%>%
  na.omit()%>% 
  unique()%>%
  mutate(species=gsub('Dipsacus asperoides','Dipsacus asper',species)) #change, otherwise problem later in representing ASVs


taxa_group=c('species_to_search_fungi','species_to_search_notfungi')
group_name=c('Fungi','Others')

all_gbif_results2=data.frame() # do it with a dataframe because couldnt manage to store element within element in list

for (i in (1:length(taxa_group))) {
  
  species_to_query <- get(taxa_group[i])

for (j in (1:length(species_to_query))) {
      print(species_to_query[j])
      
      name=name_backbone(species_to_query[j])$species
      
      gbif_results<-occ_data(
        scientificName=name, #taxon key was giving a mistake for a query so I searched by name instead
        hasGeospatialIssue = FALSE,
        occurrenceStatus = "PRESENT")
      
      gbif_results_df=gbif_results[['data']]
      
      if(is.null(gbif_results_df)) {
        print ("not found") } else { 
          
          gbif_results_df <- gbif_results_df[!(gbif_results_df$basisOfRecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN","PRESERVED_SPECIMEN")),]
          gbif_results_df <- gbif_results_df %>% mutate(group=group_name[i])
          
          all_gbif_results2<-dplyr::bind_rows(all_gbif_results2, gbif_results_df)
        }
}
}  

length(unique(all_gbif_results2$species))
length(unique(all_blast_pr2_metadata_count$species))

setdiff(unique(all_gbif_results2$species),unique(all_blast_pr2_metadata_count$species))



#### Representation of results #########

### selected ASVS across samples
asv_ids_species2=unique(asv_ids_species$qseqid[asv_ids_species$species %in% all_gbif_results2$species])
asv_ids_species2=asv_ids_species2[order(asv_ids_species2)]

asv_id_species_collapse=asv_ids_species%>%
  group_by(qseqid)%>%
  dplyr::summarise(all_species = paste(species, collapse = ","), .groups = "drop")


asv_short_names=cbind.data.frame(Feature.ID=asv_ids_species2, 
                                 short_name=paste('18S-ASV',seq(01,length(asv_ids_species2)),sep='-'),
                                 all_species=asv_id_species_collapse$all_species[match(asv_ids_species2,asv_id_species_collapse$qseqid)])


count_filt_select2=count_filt_select[rownames(count_filt_select) %in% asv_ids_species2,]
count_filt_select2=count_filt_select2[,!colnames(count_filt_select2) %in% 'total_reads']

count_filt_select2$asv_id=asv_short_names$short_name[match(rownames(count_filt_select2),asv_short_names$Feature.ID)]
count_filt_select2=melt(count_filt_select2)
count_filt_select2$lot_id3=metadata$lot_id3[match(count_filt_select2$variable,metadata$sample.id)]
count_filt_select2=count_filt_select2[!count_filt_select2$lot_id3 %in% c('EXT-NC','PCR-NC','NA'),]
count_filt_select2$country=metadata$collection_country[match(count_filt_select2$variable, metadata$sample.id)]

asvs_annotated_abundance_plot=ggplot(count_filt_select2,aes(x=asv_id,y=value,fill = interaction(lot_id3, country))) + 
  labs(y='Total reads',x='') +
  geom_col(width=.5, position = "stack")+
  scale_fill_viridis_d(option = "magma",name='Sample & Country')+
  theme_classic()+
  coord_flip()


table_grob <- tableGrob(asv_short_names[,!colnames(asv_short_names) %in% 'Feature.ID'], cols=c('ASV-ID', 'Species matched'), rows=rep('',nrow(asv_short_names)))
composed_plot1=grid.arrange(asvs_annotated_abundance_plot, table_grob, ncol = 1)

ggsave('ASV_abundance&species_matched.png',composed_plot1,width = 30, height = 20, units = "cm")

write.xlsx(asv_short_names[,!colnames(asv_short_names) %in% 'Feature.ID'], 'asv_shortname_species.xlsx')


count_filt_select2$all_species=asv_short_names$all_species[match(count_filt_select2$asv_id, asv_short_names$short_name)]

count_filt_select2 = count_filt_select2 %>%
  mutate(lot_id3 = case_when(
    lot_id3 %in% c('W04','W05','W06','W07','W08','W09') ~ "W4-9",
    lot_id3 %in% c("W13", "W14",'W15') ~ "W13-15",
    lot_id3 %in% c("W17-W", "W17-Y") ~ "W17",
    TRUE ~ lot_id3))%>%
    group_by(lot_id3, all_species) %>%
    dplyr::summarise(total = sum(value), .groups = "drop")

library(writexl)
df_list <- count_filt_select2 %>%
  filter(total != 0)%>%
  arrange(desc(total)) %>%
  group_split(lot_id3)

write_xlsx(df_list, "output2.xlsx")  




## barplots per continent ### 

#to get the continent with a customized function
get_continent_simple <- function(lat, lon) {
  if (lat >= -35 && lat <= 37 && lon >= -17 && lon <= 51) return("AFRICA")
  if (lat >= 7 && lat <= 83 && lon >= -168 && lon <= -52) return("NORTH_AMERICA")
  if (lat >= -56 && lat <= 12 && lon >= -81 && lon <= -34) return("SOUTH_AMERICA")
  if (lat >= 34 && lat <= 72 && lon >= -25 && lon <= 60) return("EUROPE")
  if (lat >= 8 && lat <= 80 && lon >= 60 && lon <= 180) return("ASIA")
  if (lat >= -50 && lat <= -10 && lon >= 110 && lon <= 180) return("OCEANIA")
  if (lat <= -60) return("ANTARCTICA")
  return(NA)
}


all_gbif_results2$continent[is.na(all_gbif_results2$continent) & !is.na(all_gbif_results2$decimalLatitude) & !is.na(all_gbif_results2$decimalLongitude)] <-
  mapply(get_continent_simple,all_gbif_results2$decimalLatitude[is.na(all_gbif_results2$continent) & !is.na(all_gbif_results2$decimalLatitude) & !is.na(all_gbif_results2$decimalLongitude)],
         all_gbif_results2$decimalLongitude[is.na(all_gbif_results2$continent) & !is.na(all_gbif_results2$decimalLatitude) & !is.na(all_gbif_results2$decimalLongitude)])


##### REPRESENTATION OF RESULTS --------------------

all_gbif_results_fungi=all_gbif_results2[all_gbif_results2$group=='Fungi',]
all_gbif_results_others=all_gbif_results2[all_gbif_results2$group=='Others',]

all_gbif_continent_fungi=all_gbif_results_fungi%>%
  group_by(species,continent)%>%
  dplyr::summarise(total_observations=n()) %>%
  na.omit() %>%
  ggplot(aes(x=continent,y=total_observations)) + 
  geom_col(aes(fill=species), width=.5, position = "stack")+
  labs(y='Number observations for all possible fungal species',x='',title = 'Fungal species GBIF DB')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


all_gbif_continent_others=all_gbif_results_others%>%
  group_by(species,continent)%>%
  dplyr::summarise(total_observations=n()) %>%
  na.omit() %>%
  ggplot(aes(x=continent,y=total_observations)) + 
  geom_col(aes(fill=species), width=.5, position = "stack")+
  labs(y='Number observations for all possible non fungal species',x='', title ='Not fungal species GBIF DB')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


#using global fungi
## Downloaded manually, all found except Alternaria passiflorae, Cochliobolus eleusines and Pleospora herbarum

write.csv2(unique(all_gbif_results_fungi$species),'selected_fungi.txt',row.names = FALSE, col.names = FALSE)

globalfung_results_list <- list.files(path = paste(path_work,"/global_fungi_results/",sep = ""), pattern = "*.txt", full.names = TRUE)
globalfung_results <- do.call(rbind, lapply(globalfung_results_list, function(file) {
  data <- read.delim(file, header = TRUE, sep = "\t")  # adjust sep as needed
  data$species <- gsub(".txt","",basename(file))
  return(data)}))


all_gbif_continent_fungi_globalfun=globalfung_results%>%
  group_by(species,Continent)%>%
  dplyr::summarise(total_observations=n()) %>%
  na.omit() %>%
  ggplot(aes(x=Continent,y=total_observations)) + 
  geom_col(aes(fill=species), width=.5, position = "stack")+
  labs(y='Number observations for all possible fungal species',x='',title = 'Fungal species Global Fungi DB')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


composite_plots2=ggarrange(all_gbif_continent_others,all_gbif_continent_fungi,all_gbif_continent_fungi_globalfun, ncol=3)  



####### representation of maps #####

world_map <- map_data("world")

base_world <- ggplot() + coord_fixed() +
  xlab("") + ylab("") +
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="grey80", fill="grey90") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="left",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())


gbif_results_geoinfo <- all_gbif_results2 %>%
    select(species,decimalLatitude,decimalLongitude) %>%
    na.omit()%>%
    unique()

level_species=na.omit(unique(gbif_results_geoinfo$species))
level_species=level_species[order(level_species)]
  
viridis_colors <- viridis(length(level_species), option = "F")  

### using GBIF ------    

  map_data <- 
    base_world +
    geom_point(data=gbif_results_geoinfo, 
               aes(x=decimalLongitude, y=decimalLatitude, fill=factor(species, levels = level_species)), colour="black", 
               pch=21, size=5, alpha=I(0.7)) +
    scale_fill_manual(values = setNames(viridis_colors, level_species), name="Possible species identity of ASVs") +
    #ggtitle(paste(names_plot[i])) + 
    guides(fill = guide_legend(ncol = 2, title.position="top"))
  
  
  ggsave(paste(path_work,"/blast_manual/manual",countries_separation[i],'_map.pdf',sep=""),map_data,device = cairo_pdf,width = 30, height = 15, units = "cm")
  

  map_data2 <- 
    base_world + 
    stat_bin_hex(data = gbif_results_geoinfo, aes(x = decimalLongitude, y = decimalLatitude), bins = 20, color = "black") +
    stat_bin_hex(data = gbif_results_geoinfo, aes(x = decimalLongitude, y = decimalLatitude, label=after_stat(count)), bins = 20, color = "black", geom="text") +
    scale_fill_gradient(low="white",high="yellow") +
    #ggtitle(paste(names_plot[i])) +
    #labs(subtitle=subtitle_rep) +
    theme(legend.position="right")
  
  ggsave(paste(path_work,"/blast_manual/manual",countries_separation[i],'_map2.pdf',sep=""),map_data2,device = cairo_pdf,
         width = 15, height = 10, units = "cm")

  
### using globalFungi ------    
    
  globalfung_map <- 
    base_world +
    geom_point(data=globalfung_results, aes(x = Longitude, y = Latitude, fill=factor(species, levels = level_species)), colour="black",pch=21,size=5, alpha=I(0.7)) +  
    #scale_fill_viridis_d(option = "H", name="Possible species identity of ASVs") +
    scale_fill_manual(values = setNames(viridis_colors, level_species), name="Possible species identity of ASVs")+
    guides(fill = guide_legend(ncol = 2, title.position="top"))
  
  globalfung_map2 <- 
    base_world + 
    stat_bin_hex(data = globalfung_results, aes(x = Longitude, y = Latitude), bins = 20, color = "black") +
    stat_bin_hex(data = globalfung_results, aes(x = Longitude, y = Latitude, label=after_stat(count)), bins = 20, color = "black", geom="text") +
    scale_fill_gradient(low="white",high="yellow") +
    theme(legend.position="right") 
  

##### write fasta and check #####

rep_seqs=read_qza('rep-seqs_clustered_18S.qza')
rep_seqs=rep_seqs$data  
  
rep_seqs_select=rep_seqs[names(rep_seqs) %in% asv_ids_species2]  
writeXStringSet(rep_seqs_select, filepath = "selected_18S_ASVs.fa")

### notes 1e1210c1c5b32ec835be5aa305aff648 match to spiders from Africa in NCBI 100%
## cd033c98ec17fa22ec9e65fe6e0ddaeb to cucumis sativus and other plants native to Asia/China region
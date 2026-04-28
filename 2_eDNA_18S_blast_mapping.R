
### SET WORKING DIRECTORY, LOAD PACKAGES AND METADATA ------

path_work="M:/postdoc_2021/sequencing_results/data_analyses/2024-2025_who/clustered_runs/18S"

setwd(path_work)

load("gbif_hist.RData")


### LOAD PACKAGES, WORKING PATHS, METADATA, 18S ASSIGNMENTS AND ABUNDANCE TABLES ------

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


# load metadata

metadata=read.xlsx("../../metadata_who_falsified_sel.xlsx", sheet= "Sheet1")
colnames(metadata)[colnames(metadata)=="sample-id"] <-"sample.id"
metadata$lot_id2_rep=paste(metadata$lot_id2,metadata$replicate,sep = "_")

rownames(metadata)=metadata$sample.id

shortcodes=read.xlsx("../../who_shortcodes.xlsx", sheet= "Sheet1")
metadata$lot_id3=shortcodes$Short.code[match(metadata$lot_id,shortcodes$PACS.code)]
metadata$lot_id3_rep=paste(metadata$lot_id3,metadata$replicate,sep = "_")
metadata=metadata[order(metadata$lot_id3_rep),]


# import filtered 18S ASV abundance table

count_filt=read.table("tables_general_stats/count_filt_18S.csv",sep = ",", header=TRUE,check.names = FALSE,row.names = 1)
count_filt$total_reads=rowSums(count_filt)


# Import PR2 lineage with metadata

lineage_with_metadata=read.xlsx('K:/ogden_grp/Conservation_Science/Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/databases/pr2/pr2_version_5.0.0_merged.xlsx',sheet=1) 
lineage_with_metadata <- lineage_with_metadata[grep("Eukaryota",lineage_with_metadata$domain),] #keep only the eukaryotic sequences
lineage_with_metadata$species=gsub (".*_sp.$", NA, lineage_with_metadata$species)


# Load blast results and process -> result is a tsv table in blastn format 6: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

all_blast=read.table("blast_18S_pe0.9_qc0.9_PR2_export/blast6.tsv",sep = "\t")
colnames(all_blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


#select ASVs/ filter counts by taxa assingment identity threshold

all_blast=all_blast[all_blast$qseqid %in% rownames(count_filt),]

threshold=98
all_blast=all_blast[all_blast$pident > threshold,]

all_blast=all_blast[order(all_blast$qseqid,-all_blast$pident),]


### SELECT BEST TAXA ASSIGNMENT HITS ------  

#best hits selection

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
  dplyr::count(qseqid, sort = TRUE) 

all_blast_ev_10=all_blast_ev[all_blast_ev$n < 6,] #select ASVs that only equally match to a maximum of 6 species


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
  unique()


### CREATION OF TABLE FOR SELECTED ASVS PER SPECIES AND THEIR ABUNDANCE ACROSS SAMPLES (FigS19) -------

asv_ids_species2=unique(asv_ids_species$qseqid)
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

write.xlsx(asv_short_names[,!colnames(asv_short_names) %in% 'Feature.ID'], 'asv_18S_select_species.xlsx')
ggsave('asv_18S_select_abunplot_species.png',asvs_annotated_abundance_plot)







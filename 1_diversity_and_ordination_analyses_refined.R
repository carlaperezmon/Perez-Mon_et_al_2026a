##### DIVERSITY AND ORDINATION ANALYSES OF WHO SAMPLES SEPT 2024 and April 2025 ========

##### set working directory, load packages and metadata ------

#path_work="M:/postdoc_2021/sequencing_results/data_analyses/2024-2025_who/clustered_runs"
path_work="M:/postdoc_2021/sequencing_results/data_analyses/2024-2025_who/clustered_runs"

setwd(path_work)

load("../1_diversity_and_ordination_analyses_refined_hist.RData")


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


packages <- c("ggplot2","openxlsx","reshape2","qiime2R","Rmisc",
              "stringr","dplyr","EcolUtils","vegan","gridExtra",
              "ggpubr","adegenet","DA","RColorBrewer","ggrepel",
              "scatterplot3d","plot3D","animation","dendextend",
              "cluster","factoextra","tidyr","sf",
              'rnaturalearth','rnaturalearthdata','ggspatial','pROC',
              'gridGraphics','ggordiplots','bios2mds','ggrepel',
              'RColorBrewer','lubridate','paletteer','MASS','ggord','scales',
              'gridExtra','patchwork','lubridate','colorspace')

ipak(packages)



# load metadata ------
metadata=read.xlsx("../metadata_who_falsified_sel.xlsx", sheet= "Sheet1")
colnames(metadata)[colnames(metadata)=="sample-id"] <-"sample.id"
metadata$lot_id2_rep=paste(metadata$lot_id2,metadata$replicate,sep = "_")

rownames(metadata)=metadata$sample.id

shortcodes=read.xlsx("../who_shortcodes.xlsx", sheet= "Sheet1")
metadata$lot_id3=shortcodes$Short.code[match(metadata$lot_id,shortcodes$PACS.code)]
metadata$lot_id3_rep=paste(metadata$lot_id3,metadata$replicate,sep = "_")
metadata=metadata[order(metadata$lot_id3_rep),]
metadata$collection_country[metadata$collection_country=='Liberia?'] <-'Liberia' #keep it like this for now, although we will need to check this info
metadata$collection_country[metadata$collection_country=='?'] <-'Unknown'
metadata$repeated[is.na(metadata$repeated)] <-'n'
metadata=metadata[metadata$seq_batch!="feb23",]
metadata=metadata[!(metadata$extract_controls=='batch24_rep' & metadata$material=='control'),]
metadata=metadata[metadata$lot_id3!='W24',]


metadata_16S=metadata[metadata$marker %in% c('16S+18S','16S'),]
metadata_18S=metadata[metadata$marker %in% c('16S+18S','18S'),]

dates=read.xlsx('../who_dates.xlsx')
dates$manuf_mont_year=paste(dates$Mfg_date_month,dates$Mfg_date_year,sep='-')
dates$manuf_mont_year=as.Date(paste0("01-", dates$manuf_mont_year), format = "%d-%b-%Y")
dates$lot_id3=metadata$lot_id3[match(dates$lot_id,metadata$lot_id)]


#parameters for representation

#lot_ids=c("NYLO","2CU1","TLZP","2T3Q","Q9RW","9UJH","X3D7","TS1C","X136","9AKZ","OEPW","1ACP","F2A4","A5J3")
lot_ids=c("W15","W02","W14","W13","W07","W06","W04","W08","W09","W10","W11","W03","W05","W12",
          "W01","W16","W17-W","W17-Y","W18","W19","W20","W21","W22","W23")

shapes_vector=c(22,21,24,3,4,23,25,7,8,9,10,11,12,13,
                0,1,2,5,6,14,15,16,17,18)

names(shapes_vector)=lot_ids

# country colors
order_countries=c("Niger","Chad","Cameroon","Nigeria","Ghana","Burkina Faso","Liberia")
color_countries<-brewer.pal(7,"Set3")
names(color_countries) <-order_countries

color_countries[2] <- "#DAA520"
color_countries[6] <- '#FF00FF'

color_countries=color_countries[order(names(color_countries))]

#hca group colors
color_function <- function (x) {if(x=="1") "blue4"
  else if(x=="2") "burlywood3"
  else if (x=="3") "azure4" 
  else if (x=="4") "chartreuse3" 
  else if (x=="5")"blueviolet"
  else if (x=="6") 'brown'
  else if (x=="7") 'yellow2'
  else if (x=="8") 'red'
  else if (x=="9") 'lightgreen'
  else if (x=="10") 'black'
}

#alternative to yellow is turquoise3

colors_hca_groups=c("blue4", "burlywood3", "azure4",'chartreuse3',"blueviolet",'brown','yellow2','red','lightgreen','black')


### date colors ####

order_dates=unique(dates$Mfg_date_year[order(dates$Mfg_date_year)])
#color_dates=hue_pal()(9)
color_dates=c('black','grey','grey','grey','grey','seashell2','seashell2','seashell2','lightpink')
names(color_dates) <-order_dates


#stylistic choices for ggplot representation  ####
size_text_y=10
size_text_x=10
size_text_label=14
size_axis_labs=10

rep_theme=theme(panel.background= element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #strip.background =element_rect(fill="white"),
                axis.text.y = element_text(size = size_text_y, color="black"),
                axis.text.x = element_text(size = size_text_x, color="black", angle=90),
                axis.title.y = element_text(size = size_text_label),
                strip.text.x = element_text(size = size_text_label, color="black"),
                legend.text = element_text(size = size_text_label),
                legend.title = element_text(size = size_text_label),
                legend.position = "rigth",
                axis.line = element_line(colour="black", linewidth= 0.5, linetype = 1), 
                axis.ticks = element_line(colour="black", linewidth= 0.5),
                axis.ticks.length = unit(0.1, "cm"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))






### REPRESENT AFRICA MAP ######
# Load African countries 
africa <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")

# Subset a few countries (e.g., Kenya, Nigeria, and South Africa)

selected_countries <- africa %>% 
  filter(admin %in% order_countries)

png(file="map_countries_collection.png")
ggplot(data = africa) +
  geom_sf(fill = "gray90", color = "black") + # Full African map
  geom_sf(data = selected_countries, fill = color_countries[sort(selected_countries$sovereignt, decreasing = TRUE)], color = "darkblue") + # Highlighted countries
  #geom_sf(data = selected_countries, fill = color_countries, color = "darkblue") + # Highlighted countries
  annotation_scale(location = "bl", width_hint = 0.5) + # Scale bar (bottom-left)
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) + # North arrow (top-left)
  coord_sf(expand = FALSE) + # Keeps coordinate limits tight
  theme_minimal() +
  labs(title = "Countries of tablet collection",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid.major = element_line(color = "gray80", size = 0.3), 
        panel.grid.minor = element_blank())+
  geom_text(data = selected_countries %>%
              st_centroid() %>%  # Get centroids of countries
              cbind(st_coordinates(.)),  # Extract coordinates
            aes(x = X, y = Y, label = name), 
            size = 5, color = "black", fontface = "bold")
dev.off()



#create objects -----
marker=c("16S","18S")
directories=c("16S","18S") #use the file with the merged sequence for 16S
objects=c("table_clustered","rep-seqs_clustered") 


#load asvs tables and repseqs (needed for later) -------

for (i in 1:length(marker)) {
  print(marker[i]) 
  
  for (j in 1:length(objects)) {
    
    data_qza=read_qza(paste(path_work,"/",directories[i],"/",objects[j],"_",marker[i],".qza", sep=""))
    data=data_qza$data
    
    assign(paste(objects[j],"_",marker[i],sep=""),data)
    rm(data_qza,data)
  }
}


hist(width(`rep-seqs_clustered_16S`))
summary(width(`rep-seqs_clustered_16S`))
`rep-seqs_clustered_16S_filt`=`rep-seqs_clustered_16S`[width(`rep-seqs_clustered_16S`) > 250] #correspond to 1st quartil
hist(width(`rep-seqs_clustered_16S_filt`))

hist(width(`rep-seqs_clustered_18S`))
summary(width(`rep-seqs_clustered_18S`))
`rep-seqs_clustered_18S_filt`=`rep-seqs_clustered_18S`[width(`rep-seqs_clustered_18S`) > 190] #correspond to 1st quartil
hist(width(`rep-seqs_clustered_18S_filt`))


writeXStringSet(`rep-seqs_clustered_18S`, 'rep-seqs_clustered_18S.fa')

#load taxonomies Silva to filter for ordination analyses -------


for (i in 1:length(marker)) {
  print(marker[i]) 
  
    data_qza=read_qza(paste(path_work,"/",directories[i],"/taxonomy_dada2_clustered_",marker[i],"_silva.qza", sep=""))
    data=data_qza$data
    
    assign(paste("taxonomy_silva_",marker[i],sep=""),data)
    rm(data_qza,data)
}


# WORK THROUGH SILVA TAXONOMY =======


for (i in 1:length(marker)) {
  get(paste("taxonomy_silva",marker[i],sep = "_")) -> taxonomy
  
  taxonomy[c("domain", "phylum", "class", "order", "family", "genus","species")] <-str_split_fixed(taxonomy$Taxon,";",7) 
 
  taxonomy=taxonomy[,-which(names(taxonomy) %in% c("Taxon"))]
  
  taxonomy=as.data.frame(apply(taxonomy, 2, function(x) gsub(" [a-z]__", "", x)))
  taxonomy=as.data.frame(apply(taxonomy, 2, function(x) gsub("[a-z]__", "", x)))
  taxonomy[taxonomy == ""] <- "Unassigned"
  
  
  taxonomy=apply(taxonomy, 2, function (x) gsub ("unclassified.*", "Unassigned", x, ignore.case=T)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub ("unidentified.*", "Unassigned", x, ignore.case=T)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub ("uncultured.*", "Unassigned", x, ignore.case=T)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub ("unknown.*", "Unassigned", x, ignore.case=T)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub ("Candidatus_", "Ca. ", x, ignore.case=T)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub (".*_or$", "Unassigned", x)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub (".*_fa$", "Unassigned", x)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub (".*_ge$", "Unassigned", x))
  taxonomy=apply(taxonomy, 2, function (x) gsub (".*_cl$", "Unassigned", x))
  taxonomy=apply(taxonomy, 2, function (x) gsub (".*_sp$", "Unassigned", x))
  taxonomy=apply(taxonomy, 2, function (x) gsub (".*_sp.$", "Unassigned", x))
  #taxonomy=apply(taxonomy, 2, function (x) gsub ("_$", "", x)) 
  taxonomy=apply(taxonomy, 2, function (x) gsub ("metagenome", "Unassigned", x, ignore.case=T))
    
    
  taxonomy=as.data.frame(taxonomy)
    
  taxonomy$genus=gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia",taxonomy$genus)
  taxonomy$genus=gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Rhizobium",taxonomy$genus)
    
 
  taxonomy_prok=taxonomy[taxonomy$domain=="Bacteria" | taxonomy$domain=="Archaea",]
  taxonomy_prok=taxonomy_prok[taxonomy_prok$family!="Mitochondria",]
  taxonomy_prok=taxonomy_prok[taxonomy_prok$family!="Chloroplast",]
    
    
  taxonomy_euk=taxonomy[taxonomy$domain=="Eukaryota" | taxonomy$domain=="Bacteria" & taxonomy$family=="Mitochondria" | taxonomy$domain=="Bacteria" & taxonomy$family=="Chloroplast",] #I did not included the reads that match to Eukaryota because there was not other classification beyond domain
  taxonomy_euk_sps=taxonomy_euk[taxonomy_euk$species!="Unassigned",]
    
  assign(paste("taxonomy_prok",marker[i],sep="_"),taxonomy_prok)
  assign(paste("taxonomy_euk",marker[i],sep="_"),taxonomy_euk)
  assign(paste("taxonomy_euk_sps",marker[i],sep="_"),taxonomy_euk_sps)
    
}



# REMOVING READS FROM CONTROLS AND CREATE RAREFIED COUNTS =======


#process negative controls ----

number_filtering=9

subDir <- "tables_general_stats"

for (i in 1:length(marker)){
  get(paste(objects[1],marker[i],sep = "_")) -> count
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  
  count=count[,colnames(count) %in% metadata_subset$sample.id]
  
  mainDir <- paste(path_work,directories[i],sep="/")
  dir.create(file.path(mainDir, subDir), showWarnings = TRUE)
  
  names_neg_control=na.omit(metadata_subset$sample.id[metadata_subset$tablet_id=="PCR_neg_control"]) #remove also the reads in the sept 2024 controls from the other A5J3 runs
  
  neg_control=count[,names_neg_control]
  neg_control=neg_control[rowSums(neg_control) > number_filtering,]
  
  counts_filt1=count[!(rownames(count) %in% rownames(neg_control)),]
  
  names_ext_control=na.omit(metadata_subset$sample.id[metadata_subset$tablet_id=="NC"])
  
  ext_control=count[,names_ext_control]
  ext_control=ext_control[rowSums(ext_control) > number_filtering,]
  
  counts_filt2=counts_filt1[!(rownames(counts_filt1) %in% rownames(ext_control)),]
  
  counts_filt2=counts_filt2[rowSums(counts_filt2) > 9,]
  
  write.csv(counts_filt2, paste(mainDir,"/", subDir,"/count_filt_",marker[i],".csv",sep="")) 
  
  
  ##### save the negative controls for later
  
  neg_ext_controls=count[,c(names_ext_control,names_neg_control)]
  neg_ext_controls=neg_ext_controls[rowSums(neg_ext_controls) > number_filtering,]
  
  colnames(neg_ext_controls)=metadata_subset$DNA_aliquot_id[match(colnames(neg_ext_controls),metadata_subset$sample.id)]
  
  write.csv(neg_ext_controls, paste(mainDir,"/",subDir,"/neg_ext_control_",marker[i],".csv",sep="")) 
  
  
  #### create filt stats table
  
  controls_filt_stats=rbind(initial_counts=colSums(count), 
                            filt_pcr_control_counts=colSums(counts_filt1), 
                            percent_remaining_counts1=round(colSums(counts_filt1)/colSums(count)*100),
                            filt_ext_control_counts=colSums(counts_filt2),
                            percent_remaining_counts2=round(colSums(counts_filt2)/colSums(count)*100))
  
  
  counts_ap=count
  counts_ap[counts_ap > 0] <- 1 
  colSums(counts_ap)
  
  asv_filt1=counts_filt1
  asv_filt1[asv_filt1 > 0] <- 1 
  colSums(asv_filt1)
  
  
  asv_filt2=counts_filt2
  asv_filt2[asv_filt2 > 0] <- 1 
  colSums(asv_filt2)
  
  
  controls_filt_stats=rbind(controls_filt_stats, 
                            inital_ASV=colSums(counts_ap),
                            filt_pcr_control_ASV=colSums(asv_filt1), 
                            percent_remaining_ASV1=round(colSums(asv_filt1)/colSums(counts_ap)*100),
                            filt_ext_control_ASV=colSums(asv_filt2),
                            percent_remaining_ASV2=round(colSums(asv_filt2)/colSums(counts_ap)*100))
  
  colnames(controls_filt_stats)=metadata_subset$lot_id3_rep[match(colnames(controls_filt_stats),metadata_subset$sample.id)]
  
  write.csv(t(controls_filt_stats), paste(mainDir,"/",subDir,"/controls_filt_stats_",marker[i],".csv",sep="")) #no much more differences between discarding all ASV reads from controls or the ones with more than 10 reads so lets be strict and remove all ASV reads from controls
  
  assign(paste("count_filt",marker[i],sep="_"),counts_filt2)
  assign(paste("neg_ext_controls",marker[i],sep="_"),neg_ext_controls)

} 



# RAREFACTION ----  

abundance_threshold=1000 #will try distinct abundance thresholds, probably also 5000, but keep in mind that some samples will get kick out
rarDir=paste(abundance_threshold,"_rarefaction", sep="")

for (i in 1:length(marker)){
  get(paste("count_filt",marker[i],sep="_")) -> count

  mainDir <- paste(path_work,directories[i],sep="/")
  dir.create(file.path(mainDir, rarDir), showWarnings = TRUE)
  
  count=count[,colSums(count) > abundance_threshold]
  
  t_count=t(count)
  print(min(rowSums(t_count))) #1129 for 16S, 2041 for 18S
 
  rel_count=decostand(t_count, method = "total")
  
  rar_count= rrarefy.perm(t_count, sample=min(rowSums(t_count)), n=100, round.out=T)
  print(rowSums(rar_count))
  
  assign(paste("rar_count",marker[i],sep="_"),rar_count)
  assign(paste("rel_count",marker[i],sep="_"),rel_count)
}


##### FIRST HIERARCHICAL CLUSTERING FOR VISUALIZATION AND CUTTING PARAMETERS SELECTION =======

#other distances -> Jensen-Shannon Divergence https://search.r-project.org/CRAN/refmans/philentropy/html/JSD.html
#Aitchison distance https://search.r-project.org/CRAN/refmans/robCompositions/html/aDist.html

#separate old WHO and newer bottles for now

all=metadata$sample.id
new=metadata$sample.id[metadata$source!='WHO-old']
old=metadata$sample.id[metadata$source=='WHO-old']

names=c('all','new','old')

subDir3 <- paste(rarDir,"/hca",sep="")

for (i in 1:length(marker)){
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  get(paste('rep-seqs_clustered',marker[i],'filt',sep = "_")) -> rep_seqs_filt
  
  for (j in 1:length(names)) {
  
  metadata_subset1=metadata_subset[metadata_subset$sample.id %in% get(names[j]),]

  count1=count[rownames(count) %in% metadata_subset1$sample.id,]
  count1=count1[,colnames(count1) %in% names(rep_seqs_filt)]
  
  if(marker[i]=='18S'){count1=count1[,colnames(count1) %in% taxonomy_euk_18S$Feature.ID]}
  
  count1<- count1[, colSums(count1 > 0) > 1]
  count1<- count1[rowSums(count1) > 0,]
  
  mainDir <- paste(path_work,directories[i],sep="/")
  dir.create(file.path(mainDir, subDir3), showWarnings = TRUE)
  
  dist.matrix=vegdist(count1,method="bray")
  hclust_avg <- hclust(dist.matrix, method = 'average')
  #hclust_avg <- hclust(vegdist(count1,method="robust.aitchison"), method = 'average')
  
  png(paste(mainDir,"/",subDir3,"/hca_",names[j],'_',marker[i],"_1.png", sep=""), width = 800, height = 480)
  plot(hclust_avg)  
  dev.off()
 
  dend <- as.dendrogram(hclust_avg)
  
  labels(dend)=metadata$lot_id3_rep[match(labels(dend),metadata$sample.id)]
  labels_colors(dend)=color_countries[match(metadata$collection_country[match(labels(dend),metadata$lot_id3_rep)],names(color_countries))]
  labels_cex(dend) <- 0.7
  dend= dend %>% set("leaves_pch", shapes_vector[match(metadata$lot_id3[match(labels(dend),metadata$lot_id3_rep)],names(shapes_vector))])
  
  dend <- sort(dend)

  plot(dend)
  
  png(paste(mainDir,"/",subDir3,"/hca_",names[j],'_',marker[i],"_2.png", sep=""), width = 1000, height = 480)
  plot(dend, main = paste(marker[i],"HCA, average distances"), pch=19)#,xlab="",sub="")
  dev.off()
  
  assign(paste('dend',names[j],marker[i],sep='_'),dend)
}

}


##### HIERARCHICAL CLUSTERING ANALYSES TO INFER SAMPLES RELATIONSHIPS (removing outliers and drawing groups) =======

subDir7 <- paste(rarDir,"/hca_lau",sep="")

#outliers from seq sep24 and april 25 
outliers_16S=c('36060OR0008L01','32800OR0065L01','36060OR0060L01','36060OR0079L01')# minimum outlier removal, new batch: '36060OR0008L01','36060OR0061L01','32800OR0065L01', old batch 
outliers_18S=c("32800OR0041L01",'36060OR0003L01','36060OR0113L01') #maybe'36060OR0139L01'

#for all group
#outliers_16S=c('36060OR0008L01','32800OR0065L01','36060OR0060L01','36060OR0061L01','36060OR0083L01','36060OR0084L01')# '36060OR0078L01','36060OR0079L01','36060OR0080L01'
outliers_18S=c("32800OR0041L01",'36060OR0003L01','36060OR0113L01',
               '36060OR0112L01','36060OR0115L01','36060OR0116L01','36060OR0139L01') 


group_W24=na.omit(metadata$sample.id[metadata$lot_id3=='W24'])
group_W22=na.omit(metadata$sample.id[metadata$lot_id3=='W22'])
group_W20=na.omit(metadata$sample.id[metadata$lot_id3=='W20'])
group_W19=na.omit(metadata$sample.id[metadata$lot_id3=='W19'])
group_W16=na.omit(metadata$sample.id[metadata$lot_id3=='W16'])

#cut values
cut_vals_new=c(0.87,0.93) #0.87 #0.83 paritions better in the distance distribution
cut_vals_old=c(0.95,0.8)
cut_vals_all=c(0.865,0.92)

### stick to all samples rep #### 
for (i in 1:length(marker)) {
  print(marker[i])
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  get(paste('rep-seqs_clustered',marker[i],'filt',sep = "_")) -> rep_seqs_filt
  get(paste('outliers',marker[i],sep = "_")) -> outliers
  
  j=1
  names[j]
  get(paste('cut_vals',names[j],sep = "_")) -> cut_vals
  metadata_subset1=metadata_subset[metadata_subset$sample.id %in% get(names[j]),]
    
    count1=count[rownames(count) %in% metadata_subset1$sample.id,]
    count1=count1[!rownames(count1) %in% outliers,]
    count1=count1[!rownames(count1) %in% group_W24,]
    count1=count1[!rownames(count1) %in% group_W20,]
    count1=count1[!rownames(count1) %in% group_W19,]
    
    if (marker[i]=='16S' & names[j]=='all') {
    count1=count1[!rownames(count1) %in% group_W16,]} #otherwise too disperse
    
    if (marker[i]=='16S') {
      count1=count1[!rownames(count1) %in% group_W22,] 
    }
    
    count1=count1[,colnames(count1) %in% names(rep_seqs_filt)]
    
    if(marker[i]=='18S'){count1=count1[,colnames(count1) %in% taxonomy_euk_18S$Feature.ID]} #filtering or not the 16S reads by taxonomy gives same patterns
    
    count1<- count1[, colSums(count1 > 0) > 1]
    count1<- count1[rowSums(count1) > 0,]
  
    mainDir <- paste(path_work,directories[i],sep="/")
    dir.create(file.path(mainDir, subDir7), showWarnings = TRUE)
  
    dist.matrix=vegdist(count1,method="bray")
    #dist.matrix=dist(count1,method="euclidean")

    #draw the dendrogram to choose a threshold
  
    hclust_avg <- hclust(dist.matrix, method = 'average')
    #hclust_avg <- hclust(dist.matrix, method = 'ward.D2')
    plot(hclust_avg)  
    
    dend <- as.dendrogram(hclust_avg)
    labels(dend)=metadata$lot_id3[match(labels(dend),metadata$sample.id)]
    labels_cex(dend) <- 0.65
    #dend=sort(dend)  
    
    plot(dend)
    
    cut_val=cut_vals[i]
    cut_avg <- cutree(hclust_avg, h = cut_val)
  
    #labels color for representation
    labels_cutavg=cut_avg
    names(labels_cutavg)=metadata_subset$lot_id3[match(hclust_avg$labels,metadata_subset$sample.id)]
    labels_cutavg=labels_cutavg[order(match(names(labels_cutavg),labels(dend)))]
    labels_cutavg=unique(as.character(labels_cutavg))
    
    color_vector=sapply(labels_cutavg, color_function)
    
    plot(dend)
    rect_positions=rect.dendrogram(dend, k=max(cut_avg), border = 'grey', xpd=0,cluster = cut_avg) #h=cut_val
    labels_position <- cumsum(c(1, lengths(rect_positions)))
    labels_position2 <- round((labels_position[-length(labels_position)] + labels_position[-1]) / 2, 0)
    
    dates_bar=color_dates[match(dates$Mfg_date_year[match(labels(dend),dates$lot_id3)],names(color_dates))]
    countries_bar=color_countries[match(metadata$collection_country[match(labels(dend),metadata$lot_id3)],names(color_countries))]
    dates_contries_bars=cbind(dates_bar,countries_bar)
    
    
    par(mar=c(5.1, 4.1, 2, 10))
    plot(dend,  ylab='Distance') #xaxt="n" 
    abline(h = cut_val, col = "blue", lty = 2, lwd = 1.5,xpd = FALSE) #plot the abline of the height you think should be correct
    rect.dendrogram(dend, k=max(cut_avg), border = 'grey', xpd=0,cluster = cut_avg) #k=max(cut_avg)
    text(x=labels_position2, y=cut_val, col=color_vector, labels=labels_cutavg, font=2, cex=2)
    colored_bars(colors = dates_contries_bars, dend = dend, sort_by_labels_order = FALSE,rowLabels=c('',''),y_shift=-0.11)
    par(xpd=TRUE)
    legend("topright", title=expression(bold("Country col.")),
           legend = names(color_countries), 
           fill= color_countries, bty = "n",  pt.cex = 0.5, cex = 1.2, x.intersp = 0.5, 
           text.col = "black", horiz = FALSE, inset = c(-0.3, 0.2),
           title.adj=c(0.1, 0.1))
    legend("topright",
           title=expression(bold("Period mfg.")),
           legend = c('2009','2011-2014','2017-2019','2021'), 
           fill= c('black','grey','seashell2','lightpink'), bty = "n",  pt.cex = 0.5, cex = 1.2, x.intersp = 0.5,
           text.col = "black", horiz = FALSE, inset = c(-0.275, 0.62),
           title.adj=c(0.1, 0.1))
    p0 <- recordPlot()
    dev.off()
    
    png(paste(mainDir,"/",subDir7,"/hca_",names[j],'_',marker[i],"_2.png", sep=""),width = 1000, height = 600)
    p0
    dev.off()
    
  #use cutree for your groupings using value between selected threshold and previous
  
  # create pairwise distance matrix and plot
    pairwise_dist <- as.data.frame(as.table(as.matrix(dist.matrix)))%>%
      dplyr::rename(distance=Freq) %>%
      filter(Var1 != Var2) %>%
      mutate(group1=cut_avg[match(Var1,names(cut_avg))],
             group2=cut_avg[match(Var2,names(cut_avg))],
             variation=ifelse(group1==group2,'Intra-Cluster','Inter-Cluster'),
             freq=1)  
    
    closest_smaller <- max(hclust_avg$height[hclust_avg$height <= cut_val])
    closest_larger <- min(hclust_avg$height[hclust_avg$height >= cut_val])
    
    
    p1=ggdensity(pairwise_dist, x = 'distance', fill = 'variation', color = 'variation',rug = TRUE)+
      labs(x = "Euclidean distance", y = "Density") +
      geom_vline(xintercept = cut_val, linetype = "dashed", color = "blue") +
      annotate("rect", xmin = closest_smaller, xmax = closest_larger, ymin = 0, ymax = Inf, alpha = 0.2, fill = "grey")+
      theme_pubr()+
      theme(legend.title = element_blank(), legend.position = c(0.2, 0.8),
            plot.title=element_text(size=10))
    
       
    roc_curve <- roc(pairwise_dist$variation, pairwise_dist$distance)  
    
    p2=ggplot(data = data.frame(tpr = roc_curve$sensitivities, fpr = 1 - roc_curve$specificities), aes(x = fpr, y = tpr)) +
      geom_line(color = "black", linewidth = 1) +
      geom_abline(linetype = "dashed", color = "red") +  # Diagonal reference line
      labs(x = "False Positive Rate", y = "True Positive Rate") +
      theme_pubr()+
      geom_text(aes(x = 0.25, y = 0.7, label = paste('AUC:',round(roc_curve$auc,2))),size=6)
    
    
    #uncostrained PCoA
    
    hca_clusters=cut_avg
    
    res.pca<- prcomp(dist.matrix,  scale = FALSE) #PCoA using distances rather than raw data, labelling the groups
    
    p3=fviz_pca_ind(res.pca, axes = c(1, 2), geom="point", habillage=hca_clusters,addEllipses=TRUE,ellipse.type='convex',ellipse.level=0.90, 
                    palette = colors_hca_groups, ggtheme =theme_pubr())+ggtitle(NULL) #ellipse.type='confidence'
    
    
    plots_below=ggarrange(p3,p1,p2,ncol = 3, labels = c('B',"C", "D"), widths = c(1, 0.6, 0.6))
    
    
    all_plots=ggarrange(p0,                                                
                        plots_below +
                          theme(plot.margin = margin(0.1, 2, 0, 0.8, "cm")), 
                        nrow = 2,
                        labels='A', heights = c(1,0.6)) + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))     
    
    ggsave(all_plots,filename=paste(mainDir,"/",subDir7,"/hca_lau_",names[j],'_',marker[i],".png", sep=""), width = 18, height = 10.5)
  
    assign(paste('hca_clusters',names[j],marker[i],sep='_'),cut_avg)
    assign(paste('hca_dend',names[j],marker[i],sep='_'),dend)
}  


### save tables as xlsx for summary of groupings #######

for (i in 1:length(marker)) {
  print(marker[i])
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  
  mainDir <- paste(path_work,directories[i],sep="/")
  dir.create(file.path(mainDir, subDir7), showWarnings = TRUE)
  
  for (j in 1:length(names)) {
    print(names[j])
    get(paste('hca_clusters',names[j],marker[i],sep = "_")) -> hca_clusters
    
    metadata_subset1=metadata_subset %>%
      select(sample.id,lot_id3,lot_id) %>%
      filter(sample.id %in% names(hca_clusters)) %>%
      mutate(group=hca_clusters[match(sample.id, names(hca_clusters))])

    metadata_subset2=metadata_subset1 %>%
      select(lot_id3,lot_id, group) %>%
      distinct()
    
    write.xlsx(metadata_subset1,paste(mainDir,"/",subDir7,"/summary_groupings_",names[j],'_',marker[i],".xlsx", sep=""))
    write.xlsx(metadata_subset2,paste(mainDir,"/",subDir7,"/summary_groupings_",names[j],'_',marker[i],"2.xlsx", sep=""))
    }
}


### NMDS for classification of too dispersed groups ######

missing_groups_16S=c('W16','W19','W20','W22')
missing_groups_18S=c('W19','W20')

subDir8=paste(subDir7,'mnds_groups_add',sep='/')

for (i in 1:length(marker)) {
  print(marker[i])
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  get(paste('rep-seqs_clustered',marker[i],'filt',sep = "_")) -> rep_seqs_filt
  get(paste('hca_clusters_all',marker[i],sep = "_")) -> hca_clusters
  get(paste('missing_groups',marker[i],sep = "_")) -> missing_groups
  
  count1=count[,colnames(count) %in% names(rep_seqs_filt)]
  
  if(marker[i]=='18S'){count1=count1[,colnames(count1) %in% taxonomy_euk_18S$Feature.ID]}
  
  count1<- count1[, colSums(count1 > 0) > 1]
  count1<- count1[rowSums(count1) > 0,]
  
  mainDir <- paste(path_work,directories[i],sep="/")
  dir.create(file.path(mainDir, subDir8), showWarnings = TRUE)
  
  count_main=count1[rownames(count1) %in% names(hca_clusters),]
  dist.matrix.main=as.matrix(vegdist(count_main,method="bray"))
  mmds_main=mmds(dist.matrix.main,pc=3)
  
  for (j in 1:length(missing_groups)) {
  sample_ids=na.omit(metadata_subset$sample.id[metadata_subset$lot_id3==missing_groups[j]])  
  count_sup=count1[rownames(count1) %in% sample_ids,]
  
  dist.matrix.sup=as.matrix(vegdist(count1,method="bray"))
  dist.matrix.sup=dist.matrix.sup[rownames(dist.matrix.sup) %in% rownames(count_sup),]
  dist.matrix.sup=dist.matrix.sup[,colnames(dist.matrix.sup) %in% rownames(count_main)]
  
  mmds_sup=mmds.project(mmds_main, dist.matrix.sup, pc = 3, group.file = NULL)
  
  #for colors and shape in plot
  
  active_vars=metadata$lot_id3[match(rownames(mmds_main$coord),metadata$sample.id)]
  shapes_active_vars=shapes_vector[match(active_vars, names(shapes_vector))]
  
  active_vars_color=hca_clusters[match(rownames(mmds_main$coord),names(hca_clusters))]
  color_vector=sapply(active_vars_color, color_function)
  
  color_vector2=color_vector
  names(color_vector2)=hca_clusters[match(names(color_vector2),names(hca_clusters))]
    
  #active_vars_country=metadata$collection_country[match(rownames(mmds_main$coord),metadata$sample.id)]
  #color_active_vars=color_vector[match(active_vars, names(shapes_vector))]
  sup_vars=metadata$lot_id3[match(rownames(mmds_sup$coord),metadata$sample.id)]
  shapes_sup_vars=shapes_vector[match(sup_vars, names(shapes_vector))]
  
  png(paste(mainDir, '/',subDir8,'/mmds_add_groups_',missing_groups[j],'_',marker[i],'.png', sep = ''),width = 1200, height = 500)
  par(mfrow = c(1, 3))
  mmds.2D.plot(mmds_main, mmds_sup, active.pch=shapes_active_vars, sup.pch=shapes_sup_vars, 
               active.col=color_vector, sup.legend.name = missing_groups[j],new.plot = FALSE,
               axis = c(1, 2))
  legend("topright", legend = unique(names(color_vector2)), col = unique(color_vector2[order(names(color_vector2))]), pch = 19)
  
  mmds.2D.plot(mmds_main, mmds_sup, active.pch=shapes_active_vars, sup.pch=shapes_sup_vars, 
               active.col=color_vector, sup.legend.name = missing_groups[j],new.plot = FALSE,
               axis = c(1, 3))
  legend("topright", legend = unique(names(color_vector2)), col = unique(color_vector2[order(names(color_vector2))]), pch = 19)
  
  
  mmds.2D.plot(mmds_main, mmds_sup, active.pch=shapes_active_vars, sup.pch=shapes_sup_vars, 
               active.col=color_vector, sup.legend.name = missing_groups[j],new.plot = FALSE,
               axis = c(2, 3))
  legend("topright", legend = unique(names(color_vector2)), col = unique(color_vector2[order(names(color_vector2))]), pch = 19)
  
  dev.off()
  
  #look directly at pairwise distances for better decision
  
  all_names=c(names(hca_clusters),sample_ids)
  count2=count1[rownames(count1) %in% all_names,]
  
  dist.matrix=vegdist(count2,method="bray")
  
  pairwise_dist <- as.data.frame(as.table(as.matrix(dist.matrix)))%>%
    rename(distance=Freq) %>%
    filter(Var1 != Var2) %>%
    mutate(Var1=metadata_subset$lot_id3[match(Var1,metadata_subset$sample.id)],
           Group2=hca_clusters[match(Var2,names(hca_clusters))],
           Var2=metadata_subset$lot_id3[match(Var2,metadata_subset$sample.id)]) %>%
    filter(Var1 == missing_groups[j]) 
  
  pairwise_dist$Group2[is.na(pairwise_dist$Group2)] <- missing_groups[j]
  
  pairwise_dist$unique_combination <- apply(pairwise_dist, 1, function(x) {
    paste(sort(c(x['Var1'], x['Group2'])), collapse = "_")
  })
  
  pairwise_dist$unique_combination <- fct_reorder(pairwise_dist$unique_combination, pairwise_dist$distance, .fun = mean)
  
  pairwise_dist2=pairwise_dist %>% group_by(unique_combination) %>%
    summarize(mean_distance = mean(distance, na.rm = TRUE),
              sd_distance = sd(distance, na.rm = TRUE))
  
  p0=ggplot(pairwise_dist, aes(x=unique_combination, y=distance))+
    geom_boxplot(color='black')+
    geom_point()
  
  hclust_avg <- hclust(dist.matrix, method = 'average')
  dend <- as.dendrogram(hclust_avg)
  labels(dend)=metadata$lot_id3[match(labels(dend),metadata$sample.id)]
  plot(dend)
  p1 <-recordPlot()
  
  png(paste(mainDir, '/',subDir8,'/pairwise_dist_hca_',missing_groups[j],'_',marker[i],'.png', sep = ''),width = 1200, height = 500)
  ggarrange(p1,p0)
  dev.off()
  
  write.xlsx(pairwise_dist2, paste(mainDir, '/',subDir8,'/pairwise_dist_means',missing_groups[j],'_',marker[i],'.xlsx', sep = ''))
  
  }

}




### EXTRA NMDS for classification of 16 groups to 4 or 6 ######

missing_groups2_16S=c('W16','W19','W22','W20')


i=1
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  get(paste('rep-seqs_clustered',marker[i],'filt',sep = "_")) -> rep_seqs_filt
  get(paste('hca_clusters_all',marker[i],sep = "_")) -> hca_clusters
  get(paste('missing_groups2',marker[i],sep = "_")) -> missing_groups
  
  count1=count[,colnames(count) %in% names(rep_seqs_filt)]
  
  if(marker[i]=='18S'){count1=count1[,colnames(count1) %in% taxonomy_euk_18S$Feature.ID]}
  
  count1<- count1[, colSums(count1 > 0) > 1]
  count1<- count1[rowSums(count1) > 0,]
  
  for (j in 1:length(missing_groups)) {
  hca_clusters_subset=hca_clusters[hca_clusters %in% c(4,6)] #groups of interest for main nmds later
  sample_ids2=na.omit(metadata_subset$sample.id[metadata_subset$lot_id3==missing_groups[j]]) #groups to interrogate as supplementary  
  count_main2=count1[rownames(count1) %in% names(hca_clusters_subset),]
  
  dist.matrix.main2=as.matrix(vegdist(count_main2,method="bray"))
  mmds_main2=mmds(dist.matrix.main2,pc=3)
  
  sample_ids2=na.omit(metadata_subset$sample.id[metadata_subset$lot_id3==missing_groups[j]])  
  count_sup2=count1[rownames(count1) %in% sample_ids2,]
  
  dist.matrix.sup2=as.matrix(vegdist(count1,method="bray"))
  dist.matrix.sup2=dist.matrix.sup2[rownames(dist.matrix.sup2) %in% rownames(count_sup2),]
  dist.matrix.sup2=dist.matrix.sup2[,colnames(dist.matrix.sup2) %in% rownames(count_main2)]
  
  mmds_sup2=mmds.project(mmds_main2, dist.matrix.sup2, pc = 3, group.file = NULL)
  
  #for colors and shape in plot
  
  active_vars2=metadata$lot_id3[match(rownames(mmds_main2$coord),metadata$sample.id)]
  shapes_active_vars2=shapes_vector[match(active_vars2, names(shapes_vector))]
  
  active_vars_color2=hca_clusters[match(rownames(mmds_main2$coord),names(hca_clusters))]
  color_vector3=sapply(active_vars_color2, color_function)
  
  color_vector4=color_vector3
  names(color_vector4)=hca_clusters[match(names(color_vector4),names(hca_clusters))]
  
  active_vars_country=metadata$collection_country[match(rownames(mmds_main$coord),metadata$sample.id)]
  color_active_vars=color_vector[match(active_vars, names(shapes_vector))]
  sup_vars2=metadata$lot_id3[match(rownames(mmds_sup2$coord),metadata$sample.id)]
  shapes_sup_vars2=shapes_vector[match(sup_vars2, names(shapes_vector))]
  
  png(paste(mainDir, '/',subDir8,'/mmds_add_groups_',missing_groups[j],'_',marker[i],'_subset.png', sep = ''),width = 1200, height = 500)
  par(mfrow = c(1, 3))
  mmds.2D.plot(mmds_main2, mmds_sup2, active.pch=shapes_active_vars2, sup.pch=shapes_sup_vars2, 
               active.col=color_vector3, sup.legend.name = missing_groups[j],new.plot = FALSE,
               axis = c(1, 2))
  legend("topright", legend = unique(names(color_vector2)), col = unique(color_vector2[order(names(color_vector2))]), pch = 19)
  
  mmds.2D.plot(mmds_main2, mmds_sup2, active.pch=shapes_active_vars2, sup.pch=shapes_sup_vars2, 
               active.col=color_vector3, sup.legend.name = missing_groups[j],new.plot = FALSE,
               axis = c(1, 3))
  legend("topright", legend = unique(names(color_vector2)), col = unique(color_vector2[order(names(color_vector2))]), pch = 19)
  
  
  mmds.2D.plot(mmds_main2, mmds_sup2, active.pch=shapes_active_vars2, sup.pch=shapes_sup_vars2, 
               active.col=color_vector3, sup.legend.name = missing_groups[j],new.plot = FALSE,
               axis = c(2, 3))
  legend("topright", legend = unique(names(color_vector2)), col = unique(color_vector2[order(names(color_vector2))]), pch = 19)
  dev.off()
  }  





### EXPLORING RELATIONSHIP WITH DATES AND PCAS BASED ON TIME vs PLACE ####### 

#function for color
darken_sequence <- function(base_color, n = 7) {
  # Convert to HCL
  hcl_col <- as(hex2RGB(base_color), "polarLUV")
  hue <- coords(hcl_col)[, "H"]
  chroma <- coords(hcl_col)[, "C"]
  
  # Lightness sequence from light to dark
  lightness_seq <- seq(90, 30, length.out = n)
  
  # Create colors
  hcl(h = hue, c = chroma, l = lightness_seq)
}



for (i in 1:length(marker)) {
  print(marker[i])
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste('metadata',marker[i],sep = "_")) -> metadata_subset
  get(paste('rep-seqs_clustered',marker[i],'filt',sep = "_")) -> rep_seqs_filt
  
  count1=count[,colnames(count) %in% names(rep_seqs_filt)]
  count1=count[rownames(count) %in% metadata_subset$sample.id,]
  
  if(marker[i]=='18S'){count1=count1[,colnames(count1) %in% taxonomy_euk_18S$Feature.ID]}
  
  count1<- count1[, colSums(count1 > 0) > 1]
  count1<- count1[rowSums(count1) > 0,]
  
  #count1=count1[rownames(count1) %in% names(hca_clusters),]
  dist.matrix=as.matrix(vegdist(count1,method="bray"))
  
  pairwise_dist <- as.data.frame(as.table(as.matrix(dist.matrix)))%>%
    dplyr::rename(distance=Freq) %>%
    filter(Var1 != Var2) %>%
    mutate(Var1=metadata_subset$lot_id3[match(Var1,metadata_subset$sample.id)],
           Var2=metadata_subset$lot_id3[match(Var2,metadata_subset$sample.id)],
           lot1=metadata$lot_id[match(Var1,metadata$lot_id3)],
           lot2=metadata$lot_id[match(Var2,metadata$lot_id3)],
           manuf_year1=dates$Mfg_date_year[match(lot1,dates$lot_id)],
           manuf_year2=dates$Mfg_date_year[match(lot2,dates$lot_id)],
           manuf_date1=dates$manuf_mont_year[match(lot1,dates$lot_id)],
           manuf_date2=dates$manuf_mont_year[match(lot2,dates$lot_id)],
           col_place1=metadata$collection_country[match(lot1,metadata$lot_id)],
           col_place2=metadata$collection_country[match(lot2,metadata$lot_id)],
           diff_manuf_date=abs(time_length(interval(manuf_date1, manuf_date2), unit = "months")))
  
  
  pairwise_dist$pairs_time <- apply(pairwise_dist, 1, function(x) {
    paste(sort(c(x['manuf_year1'], x['manuf_year2'])), collapse = "_")
  })
  
  pairwise_dist$pairs_countries <- apply(pairwise_dist, 1, function(x) {
    paste(sort(c(x['col_place1'], x['col_place2'])), collapse = " - ")
  })  
  
 
  
  pairwise_dist2 <- pairwise_dist%>%
    #filter(col_place1 != col_place2)%>%
    dplyr::select(distance,pairs_time,pairs_countries)%>%
    unique()
  
  
  pairwise_dist2$year_block=sapply(strsplit(as.character(pairwise_dist2$pairs_time), "_"), `[`, 1)
  pairwise_dist2$country_block=sapply(strsplit(pairwise_dist2$pairs_countries, " - "), `[`, 1)
  
  
  my_colors=c()
  for(k in 1:length(color_countries)){
    countries_vector=unique(pairwise_dist2$pairs_countries[pairwise_dist2$country_block==names(color_countries[k])])
    countries_vector=countries_vector[order(countries_vector)]
    
    base_col <- color_countries[k]
    
    palette_color <- darken_sequence(base_col, length(countries_vector))
    
    names(palette_color)=countries_vector
    
    my_colors=c(my_colors,palette_color)
  }
  
  
  rep_time=ggplot(pairwise_dist2, aes(x=pairs_time, y=distance))+
    geom_boxplot(color='black')+
    geom_point(aes(color=pairs_countries))+
    scale_color_manual(values = my_colors, name='Country pairs')+
    labs(x='',y='Distance') + #'Pairwise distances between manufacturing years, grouped by year of manufacture'
    facet_grid(~year_block,scales='free_x',space = "free_x")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 10))
  
  pairwise_dist3=pairwise_dist%>%
    filter(col_place1 != col_place2)%>%
    dplyr::select(distance,pairs_countries)%>%
    group_by(pairs_countries)%>%
    mutate(median_dist=median(distance))%>%
    arrange(median_dist)%>%
    unique()
  
  order_pairs_countries=unique(pairwise_dist3$pairs_countries)
  
  pairwise_dist3$country_block=sapply(strsplit(pairwise_dist3$pairs_countries, " - "), `[`, 1)
  
  rep_place=ggplot(pairwise_dist3, aes(x=factor(pairs_countries, levels=order_pairs_countries), y=distance))+
    facet_grid(cols=vars(country_block),scales='free_x',space = "free_x")+
    geom_boxplot(color='black')+
    geom_point(aes(color=pairs_countries))+
    labs(x='',y='Distance')+
    scale_color_manual(values = my_colors, name='Country pairs',guide = "none")+
    labs(x='')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")+
    theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 10))
  
  ## PCAs comparing place of collection or time ###
  
  pcoa=cmdscale(dist.matrix,eig=TRUE, k=3)
  
  pcoa_df=as.data.frame(pcoa$points[, 1:3])
  
  Eigenvalues <- eigenvals(pcoa) 
  Variance <- Eigenvalues / sum(Eigenvalues) 
  Variance1 <- 100 * signif(Variance[1], 2)
  Variance2 <- 100 * signif(Variance[2], 2)
  Variance3 <- 100 * signif(Variance[3], 2)
  
  pcoa_df$lot_id3=metadata_subset$lot_id3[match(rownames(pcoa_df),metadata_subset$sample.id)]
  pcoa_df$source=metadata_subset$source[match(pcoa_df$lot_id3,metadata_subset$lot_id3)]
  pcoa_df$collection_country=metadata_subset$collection_country[match(pcoa_df$lot_id3,metadata_subset$lot_id3)]
  pcoa_df$mfg_year=as.factor(dates$Mfg_date_year[match(pcoa_df$lot_id3,dates$lot_id3)])
  pcoa_df=pcoa_df[order(pcoa_df$lot_id3),]
  
  pcoa_time=ggplot(pcoa_df, aes(V1,V2))+ 
    geom_point(size=5, aes(shape=lot_id3, 
                           color=mfg_year,
                           fill=mfg_year)) +
    scale_color_manual(values=color_dates, name="Year of manufacture") +
    scale_fill_manual(values=color_dates, guide="none") +
    scale_shape_manual(values=shapes_vector, name="Samples")+
    xlab(paste("Axis 1 (", Variance1, "% )")) + 
    ylab(paste("Axis 2 (", Variance2, "% )")) +  
    theme_classic()+
    guides(shape = guide_legend(ncol = 4))
  
  
  pcoa_place=ggplot(pcoa_df, aes(V1,V2))+ 
    geom_point(size=5, aes(shape=lot_id3, 
                           color=collection_country,
                           fill=collection_country)) +
    scale_color_manual(values=color_countries, name="Country of collection") +
    scale_fill_manual(values=color_countries, guide="none") +
    scale_shape_manual(values=shapes_vector, name="Samples")+
    xlab(paste("Axis 1 (", Variance1, "% )")) + 
    ylab(paste("Axis 2 (", Variance2, "% )")) +  
    theme_classic()+
    guides(shape = guide_legend(ncol = 4))
  
  mainDir <- paste(path_work,directories[i],sep="/")
  ggsave(  
    free(rep_time,type = "label") + free(pcoa_time,type = "label") + free(rep_place,type = "label") + free(pcoa_place,type = "label") +
      plot_layout(ncol = 2) +
      plot_layout(widths = c(3, 1,3,1))+
      plot_layout(guides = "collect")  & guides(color = guide_legend(ncol = 3)), 
    filename=paste(mainDir,"/",subDir7,"/pcoa_time_place_",names[j],'_',marker[i],".png", sep=""), width = 20, height = 10.5)
  
}
  
  
### general taxonomy tables ######
  
t_rar_count_16S=as.data.frame(t(rar_count_16S))
colnames(t_rar_count_16S)=metadata$lot_id3_rep[match(colnames(t_rar_count_16S),metadata$sample.id)]
t_rar_count_16S=t_rar_count_16S[, order(colnames(t_rar_count_16S))]
t_rar_count_16S$Feature.ID=rownames(t_rar_count_16S)
taxonomy_counts_16S=merge(taxonomy_silva_filt_16S, t_rar_count_16S, by='Feature.ID')
write.xlsx(taxonomy_counts_16S, 'taxonomy_counts_16S.xlsx')  

t_rar_count_18S=as.data.frame(t(rar_count_18S))
colnames(t_rar_count_18S)=metadata$lot_id3_rep[match(colnames(t_rar_count_18S),metadata$sample.id)]
t_rar_count_18S=t_rar_count_18S[, order(colnames(t_rar_count_18S))]
t_rar_count_18S$Feature.ID=rownames(t_rar_count_18S)
taxonomy_counts_18S=merge(taxonomy_silva_filt_18S, t_rar_count_18S, by='Feature.ID')
write.xlsx(taxonomy_counts_18S, 'taxonomy_counts_18S.xlsx')  

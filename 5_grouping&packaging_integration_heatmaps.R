
### SET WORKING DIRECTORY, LOAD PACKAGES, METADATA AND GROUPING DATA --------

path_work="M:/postdoc_2021/who_extra_analyses/all_groupings_comparison"

setwd(path_work)

load("all_groupings_comparison_inclpackaging.RData")

## load packages ###

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


packages <- c("ggplot2","openxlsx","reshape2","dplyr",'gplots','RColorBrewer','circlize','dendextend',
              'stringr','rnaturalearth','rnaturalearthdata','ggspatial')

ipak(packages)

#load metadata

metadata=read.xlsx('../who_metadata.xlsx')

dates=read.xlsx('../who_dates.xlsx')
dates$manuf_mont_year=paste(dates$Mfg_date_month,dates$Mfg_date_year,sep='-')
dates$manuf_mont_year=as.Date(paste0("01-", dates$manuf_mont_year), format = "%d-%b-%Y")
dates$lot_id3=metadata$lot_id3[match(dates$lot_id,metadata$lot_id)]
dates$collection_country=metadata$collection_country[match(dates$lot_id,metadata$lot_id)]
dates=na.omit(dates)

morph=read.xlsx('morphometry_summary_hca_groupings.xlsx') %>%
  dplyr::select(lot_id3,hca_grouping)%>%
  rename(morph_groups=hca_grouping)%>%
  distinct()

nir=read.xlsx('NIR_Data_Clusters.xlsx')  %>%
  dplyr::select(lot_id3,group) %>%
  dplyr::rename(nir_groups=group)%>%
  distinct()

ftir= read.xlsx('ftir_summary_hca_groupings.xlsx') %>%
  dplyr::select(lot_id3,hca_grouping) %>%
  rename(ftir_groups=hca_grouping)%>%
  distinct()

xrd=read.xlsx('xrd_summary_hca_groupings.xlsx') %>%
  dplyr::select(lot_id3,hca_grouping) %>%
  rename(xrd_groups=hca_grouping)%>%
  distinct()%>%
  na.omit()

isotopes=read.xlsx('isotopes_summary_hca_groupings.xlsx') %>%
  dplyr::select(lot_id3,hca_grouping)%>%
  rename(iso_groups=hca_grouping)%>%
  distinct() #in the file I outgrouped W17

hydro=read.xlsx('hydro_summary_hca_groupings.xlsx')%>%
  dplyr::select(lot_id3,hca_grouping)%>%
  rename(hydro_groups=hca_grouping)%>%
  distinct()%>%
  na.omit()

caco3=read.xlsx('caco3_summary_hca_groupings.xlsx')  %>%
  dplyr::select(lot_id3,hca_grouping) %>%
  dplyr::rename(caco3_groups=hca_grouping)%>%
  distinct()

icpms=read.xlsx('summary_hca_groupings_icpms_edited.xlsx')  %>%
  dplyr::select(lot_id3,hca_grouping) %>%
  dplyr::rename(icpms_groups=hca_grouping)%>%
  distinct()

eDNA_16S=read.xlsx('summary_groupings_all_16S_final.xlsx')  %>%
  dplyr::select(lot_id3,group) %>%
  tidyr::separate_rows(group, sep = ",")%>%
  mutate(group=as.numeric(group))%>%
  mutate(lot_id3=make.unique(lot_id3, sep = ".")) %>% #make a note in figure that W19 can be 4 or 6
  dplyr::rename(all_16S_groups=group)%>%
  distinct()

api=read.xlsx('summary_groupings_api.xlsx') %>% #groups created manually by identifying presence/absence of compounds (Orbitrap LC-MS data)
  dplyr::select(lot_id3,group) %>%
  tidyr::separate_rows(group, sep = ",")%>%
  mutate(group=as.numeric(group))%>%
  mutate(lot_id3=make.unique(lot_id3, sep = ".")) %>% #make a note in figure that W19 can be 4 or 6
  dplyr::rename(all_ipas_groups=group)%>%
  distinct()

exci=read.xlsx('excip_summary_groupings.xlsx') #groups created manually by identifying presence/absence of compounds (XRD data)


# create a unique table with all groups per analyses across all samples

table_all <- merge(morph,nir,by=c('lot_id3'))
table_all <- merge(table_all,ftir,by=c('lot_id3'))
table_all <- merge(table_all,xrd,by=c('lot_id3'))
table_all <- merge(table_all,isotopes, by=c('lot_id3'))
table_all <- merge(table_all,hydro,by=c('lot_id3'))
table_all <- merge(table_all,caco3, by=c('lot_id3'))
table_all <- merge(table_all,icpms, by=c('lot_id3'))
table_all <- merge(table_all,eDNA_16S, by=c('lot_id3'))
table_all <- merge(table_all,ipa, by=c('lot_id3'))
table_all <- merge(table_all,exci, by=c('lot_id3'))


tables_compare=c('table_all')
titles=c("MORPH, NIR,FTIR, XRD, ISO-CNO, ISO-H, ISO-caco3,ICPMS,16S,IPAs,exci")
inclusion=c('all_tablets')


### PARAMETERS FOR REPRESENTATION ------

# date colors

order_dates=unique(dates$Mfg_date_year[order(dates$Mfg_date_year)])
color_dates=c('black','grey','grey','grey','grey','seashell2','seashell2','seashell2','seashell2')
names(color_dates) <-order_dates


# country colors
order_countries=c("Niger","Chad","Cameroon","Nigeria","Ghana","Burkina Faso","Liberia")
country_abbreviations=c('NE','TD','CM','NG','GH','BF','LR')
color_countries<-brewer.pal(7,"Set3")
names(color_countries) <-order_countries

color_countries[2] <- "goldenrod"
color_countries[6] <- 'magenta'

color_countries=color_countries[order(names(color_countries))]


### HEATMAP OF SHARED GROUPS FOR EACH PAIR OF SAMPLES, ACROSS ALL LABORATORY ANALYSES (Fig 4A, lower triangle) ------

# create crosstable for all shared group memberships

i=1
print(tables_compare[i])
get(tables_compare[i]) -> table_chosen

table_all2=table_chosen  %>% 
tibble::column_to_rownames(var = "lot_id3")
  
  table_all2= table_all2 %>%
    mutate(across(everything(), ~ LETTERS[.x])) %>%
    t()%>% as.data.frame()
  
  
  # Initialize an empty matrix to store the results
  col_names <- colnames(table_all2)
  cross_table <- matrix(0, ncol = ncol(table_all2), nrow = ncol(table_all2))
  rownames(cross_table) <- col_names
  colnames(cross_table) <- col_names
  
  # Calculate shared occurrences for each pair of columns
  combinations <- combn(col_names, 2, simplify = FALSE)
  
  for (pair in combinations) {
    # Extract the two columns
    col1 <- table_all2[[pair[1]]]
    col2 <- table_all2[[pair[2]]]
    
    # Count the number of rows with the same values in both columns
    shared_count <- sum(col1 == col2)
    
    # Assign the count symmetrically in the matrix
    cross_table[pair[1], pair[2]] <- shared_count
    cross_table[pair[2], pair[1]] <- shared_count
  }
  

cross_table_dist <- 1 - (cross_table / max(cross_table))
cross_table_dist <- as.dist(cross_table_dist)  # converts to "dist" object

col_hclust_dist=hclust(cross_table_dist, method='single')
col_hclust=as.dendrogram(col_hclust_dist)

order_groups=unique(dates$lot_id3[order(dates$Mfg_date_year,dates$collection_country)])
order_groups=order_groups[order_groups %in% labels(col_hclust)]

col_hclust <- rotate(col_hclust, order = order_groups)
plot(col_hclust)

order_labels2=c("W18","W23","W19","W11","W16","W22","W10","W21","W17-W","W17-Y","W01","W20","W13","W14","W15","W12","W03","W02","W04",
                "W05","W06","W07","W08","W09")

col_hclust2 <- rotate(col_hclust, order = order_labels2)
plot(col_hclust2)



# represent as heatmap

custom_row_names <- paste(rownames(cross_table),dates$Mfg_date_year[match(rownames(cross_table),dates$lot_id3)], sep='  ')
custom_col_names <- paste(metadata$col_country_ab[match(colnames(cross_table),metadata$lot_id3)],colnames(cross_table), sep='  ')

order_dates_rep=dates$Mfg_date_year[match(rownames(cross_table),dates$lot_id3)]
order_countries_rep=metadata$collection_country[match(colnames(cross_table),metadata$lot_id3)]

dates_bar=color_dates[match(order_dates_rep, names(color_dates))]
countries_bar=color_countries[match(order_countries_rep,names(color_countries))]
dates_contries_bars=as.matrix(cbind(dates_bar,countries_bar))


png(paste(inclusion[i],'heatmap_common_groupings.png',sep='_'),width = 600,height = 600)
par(mar = c(4, 4, 4, 4))
heatmap.2(cross_table,
        scale = "none",   
        Colv = col_hclust2,
        Rowv = col_hclust2,
        col = colorRampPalette(c("white", 'orange',"red")),
        key = FALSE, 
        trace="none",
        breaks = seq(0, max(cross_table),0.5),
        density.info = "none",
        sepcolor = "grey",  # Color of the separation lines
        rowsep = 0:nrow(cross_table),       # Lines between each row
        colsep = 0:ncol(cross_table),       # Lines between each column
        sepwidth = c(0.01, 0.01),
        lwid = c(1, 4),      # Width for layout adjustments
        lhei = c(1, 4),
        cellnote = ifelse(cross_table == 0, "", round(cross_table, 1)),
        notecol = "black",
        RowSideColors=countries_bar,
        ColSideColors=dates_bar,
        margins = c(7, 7))
dev.off()


### HEATMAP OF SHARED GROUPS FOR EACH PAIR OF SAMPLES, ACROSS PACKAGING ANALYSES (Fig 4A, upper triangle) ------

packaging=read.xlsx('summary_groupings_packaging3.xlsx')%>% 
  tibble::column_to_rownames('lot_id3')%>% 
  mutate(across(everything(), as.factor))%>%
  select(-cavity)

packaging2 = packaging %>%
  mutate(across(everything(), ~ LETTERS[.x])) %>%
  t()%>% as.data.frame()

# Initialize an empty matrix to store the results
col_names <- colnames(packaging2)
cross_table_pack <- matrix(0, ncol = ncol(packaging2), nrow = ncol(packaging2))
rownames(cross_table_pack) <- col_names
colnames(cross_table_pack) <- col_names

# Calculate shared occurrences for each pair of columns
combinations <- combn(col_names, 2, simplify = FALSE)

for (pair in combinations) {
  # Extract the two columns
  col1 <- packaging2[[pair[1]]]
  col2 <- packaging2[[pair[2]]]
  
  # Count the number of rows with the same values in both columns
  shared_count <- sum(col1 == col2)
  
  # Assign the count symmetrically in the matrix
  cross_table_pack[pair[1], pair[2]] <- shared_count
  cross_table_pack[pair[2], pair[1]] <- shared_count
}

cross_table_dist_pack <- 1 - (cross_table_pack/ max(cross_table_pack))
cross_table_dist_pack <- as.dist(cross_table_dist_pack)  # converts to "dist" object

col_hclust_dist_pack=hclust(cross_table_dist_pack, method='single')
col_hclust_pack=as.dendrogram(col_hclust_dist_pack)

plot(col_hclust_pack)

col_hclust_pack %>% sort() %>% plot()

order_names=c("W01","W02","W03","W04","W05","W06","W07-W","W07-Y","W08","W09","W12","W13","W14","W15","W16","W20","W10","W19","W21","W17","W11","W18","W22","W23")

col_hclust_pack2=rotate(col_hclust_pack,order=order_names)

plot(col_hclust_pack2)


#### represent as heatmap ####

custom_row_names <- paste(rownames(cross_table_pack),dates$Mfg_date_year[match(rownames(cross_table_pack),dates$lot_id3)], sep='  ')
custom_col_names <- paste(metadata$col_country_ab[match(colnames(cross_table_pack),metadata$lot_id3)],colnames(cross_table_pack), sep='  ')

order_dates_rep=dates$Mfg_date_year[match(rownames(cross_table_pack),dates$lot_id3)]
order_countries_rep=metadata$collection_country[match(colnames(cross_table_pack),metadata$lot_id3)]

dates_bar=color_dates[match(order_dates_rep, names(color_dates))]
countries_bar=color_countries[match(order_countries_rep,names(color_countries))]
dates_contries_bars=as.matrix(cbind(dates_bar,countries_bar))

png(paste(inclusion[i],'heatmap_common_groupings_packaging.png',sep='_'),width = 600,height = 600)
par(mar = c(4, 4, 4, 4))
heatmap.2(cross_table_pack,
          scale = "none",   
          Colv = col_hclust2,
          Rowv = col_hclust2,
          col = colorRampPalette(c("white", 'lightblue',"blue")),
          key = FALSE, 
          trace="none",
          breaks = seq(0, max(cross_table_pack),0.5),
          density.info = "none",
          sepcolor = "grey",  # Color of the separation lines
          rowsep = 0:nrow(cross_table_pack),       # Lines between each row
          colsep = 0:ncol(cross_table_pack),       # Lines between each column
          sepwidth = c(0.01, 0.01),
          lwid = c(1, 4),      # Width for layout adjustments
          lhei = c(1, 4),
          cellnote = ifelse(cross_table_pack == 0, "", round(cross_table_pack, 1)),
          notecol = "black",
          RowSideColors=countries_bar,
          ColSideColors=dates_bar,
          margins = c(7, 7))
dev.off()



### REPRESENTATION OF AFRICA MAP (Fig4B) -----

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



### sum of values to represent lines in map (lines added manually in ppt)

cross_table2=melt(cross_table) %>% 
  filter(Var1 != Var2) %>% 
  mutate(country1=metadata$collection_country[match(as.character(Var1),metadata$lot_id3)],
         country2=metadata$collection_country[match(as.character(Var2),metadata$lot_id3)],
         Var1 = as.character(Var1), Var2 =as.character(Var2),
         date_var1=dates$Mfg_date_year[match(Var1,dates$lot_id3)]) %>% 
  dplyr::select(country1,country2, everything()) %>% 
  arrange(country1,date_var1,Var1)


cross_table2$Var1[cross_table2$Var1 %in% c('W04','W05','W06','W07','W08','W09')] <- 'W04-09'
cross_table2$Var2[cross_table2$Var2 %in% c('W04','W05','W06','W07','W08','W09')] <- 'W04-09'

cross_table2$Var1[cross_table2$Var1 %in% c('W13','W14','W15')] <- 'W13-15'
cross_table2$Var2[cross_table2$Var2 %in% c('W13','W14','W15')] <- 'W13-15'


#accommodate the data##
cross_table4=cross_table2 %>%
  distinct()

cross_table4$pairs_countries <- apply(cross_table4, 1, function(x) {
  paste(sort(c(x['country1'], x['country2'])), collapse = " - ")
})

cross_table4 = cross_table4 %>%
  dplyr::select(pairs_countries,value)%>%
  distinct()%>%
  group_by(pairs_countries)%>%
  dplyr::summarise(sum = sum(value, na.rm = TRUE))%>%
  mutate(country1=word(pairs_countries,1,sep=' - '),
         country2=word(pairs_countries,2,sep=' - '))%>%
  dplyr::select(country1,country2,sum)%>%
  rename(from=country1,to=country2)%>%
  filter(sum>0)%>%
  filter(from!=to)










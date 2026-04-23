# UNIFIED SCRIPT FOR REPRESENTATION OF HCA FOR ALL SIDE DATA ###
### general working directories###

path_work="M:/postdoc_2021/who_extra_analyses/"

setwd(path_work)

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

packages <- c("ggplot2","openxlsx","reshape2","dplyr",
              "ggpubr","dendextend","factoextra","tidyr",
              'pROC','dplyr','tibble','HotellingEllipse','ggpubr','ggrepel',
              'RColorBrewer','lubridate','paletteer','MASS','ggord','scales',
              'gridExtra','patchwork','vegan','colorspace')

ipak(packages)

load("all_hcas_hist.RData")

# load data ------
metadata=read.xlsx("who_metadata.xlsx", sheet= "Sheet1")

dates=read.xlsx('who_dates.xlsx')
dates$manuf_mont_year=paste(dates$Mfg_date_month,dates$Mfg_date_year,sep='-')
dates$manuf_mont_year=as.Date(paste0("01-", dates$manuf_mont_year), format = "%d-%b-%Y")
dates$lot_id3=metadata$lot_id3[match(dates$lot_id,metadata$lot_id)]

##### weight #####

data_weight=read.xlsx("who_morphometry/tablets_weight.xlsx", sheet= "Sheet1")

data_weight2=data_weight %>%
  group_by(lot_id3) %>%
  mutate(lot_id3_rep=paste0(lot_id3, "_", row_number())) %>%
  dplyr::select(lot_id3, lot_id3_rep, everything())%>%
  filter(lot_id3_rep!='W16_3') %>%  ##remove w16_3 so we do not have an outlier in the HCAs
  column_to_rownames(var = "lot_id3_rep") %>%
  dplyr::select(-lot_id3, -measurement)

data_scaled_weight<-scale(data_weight2) 
mean(data_scaled_weight[,1]) #mean 0
sd(data_scaled_weight[,1])


#### morph #############
data_morph=read.xlsx("who_morphometry/allmorphology_forR.xlsx", sheet= "Sheet1") %>%
  pivot_longer(cols = c(tablet_1, tablet_2, tablet_3, tablet_4,tablet_5,tablet_6), 
               names_to = "tablet_rep", values_to = "diameter") %>%
  mutate(lot_id2 = sub("-[^-]*$", "", lot_id),
         rep=gsub('tablet_','',tablet_rep),
         lot_id3_rep=paste(lot_id3,rep,sep='_'))%>%
  filter(lot_id2 != 'W24')

data_morph2 <- data_morph %>%
  mutate(score_line_type = case_when(
    score_line_type == 'none' ~ 0,
    score_line_type == 'break' ~ 1,
    score_line_type == 'score' ~ 2,
    score_line_type == 'crossed' ~ 3),
    tablet_shape= case_when(
      tablet_shape == 'convex' ~ 0,
      tablet_shape == 'Flat on top and bottom' ~ 1))%>%
  column_to_rownames(var = "lot_id3_rep") %>%
  dplyr::select(diameter,tablet_shape,score_line_type) #score_line_present #score_line_type #tablet_shape

data_scaled_morph=scale(data_morph2,center = TRUE, scale = TRUE) #so that all variables in a same scale
mean(data_scaled_morph[,1]) #mean 0
sd(data_scaled_morph[,1]) #sd 1


###alternatively including weight

data_morph3=data_morph2
data_morph3$weight=data_weight2$weight[match(rownames(data_morph3),rownames(data_weight2))]

data_scaled_morph3=scale(data_morph3,center = TRUE, scale = TRUE) #so that all variables in a same scale
mean(data_scaled_morph3[,1]) #mean 0
sd(data_scaled_morph3[,1]) #sd 1


### FTIR ######

data_ftir_ad=read.xlsx("who_ftir/WHO EDP samples ATR-FT-IR SCREEN Compiled results.xlsx", sheet= "Sheet1")

data_ftir_ad1=read.xlsx("who_ftir/8Tabs_FTIR_4reps_4SK_edited.xlsx", sheet= "Sheet1")
data_ftir_ad1=data_ftir_ad1[data_ftir_ad1$SpleID!='AMOX-VYTH-T',]

data_ftir_ad2=read.xlsx("who_ftir/Additional 2 tabs_27-06-2025_edited.xlsx", sheet= "Sheet1")

data_ftir=rbind(data_ftir_ad,data_ftir_ad1,data_ftir_ad2)

data_ftir<- data_ftir %>%
  mutate(lot_id3=metadata$lot_id3[match(SpleID,metadata$lot_id)])%>%
  group_by(lot_id3) %>%
  mutate(lot_id3_rep=paste0(lot_id3, "_", row_number())) %>%
  dplyr::select(lot_id3, lot_id3_rep, everything())

data_ftir2=as.data.frame(data_ftir[,!colnames(data_ftir) %in% c("lot_id3","PrimID","SpleID")])
rownames(data_ftir2)=data_ftir2$lot_id3_rep
data_ftir2=data_ftir2[,!colnames(data_ftir2) %in% c("lot_id3_rep")]


data_scaled_ftir<-scale(t(data_ftir2),center=TRUE,scale=TRUE) 
mean(data_scaled_ftir[,1]) #mean 0
sd(data_scaled_ftir[,1]) #sd 1

data_scaled_ftir=t(data_scaled_ftir) #transpose again so that samples are row and spectra peaks are columns for representation


## isotopes #####

data_iso=read.xlsx("who_isotopes/WHO_LOGO_IRMS_edited.xlsx", sheet= "Sheet1")

data_iso <- data_iso %>%
  mutate(lot_id3=metadata$lot_id3[match(sample.id,metadata$lot_id)])%>%
  group_by(lot_id3,element) %>%
  mutate(lot_id3_rep=paste0(lot_id3, "_", row_number())) %>%
  dplyr::select(-measurement) %>%
  pivot_wider(names_from = element, values_from = value) %>%
  filter(lot_id3!='W24')

data_iso2 <- data_iso %>%  
  group_by(lot_id3_rep) %>%
  column_to_rownames(var = "lot_id3_rep") %>%
  dplyr::select(-lot_id3, -sample.id, -S34) # stick to CNO for now


## means for table 

data_iso3=data_iso %>% mutate(lot_id3 = case_when(
  lot_id3 %in% c('W04','W05','W06','W07','W08','W09') ~ "W4-9",
  lot_id3 %in% c("W13", "W14",'W15') ~ "W13-15",
  lot_id3 %in% c("W17-W", "W17-Y") ~ "W17",
  TRUE ~ lot_id3))

data_iso_mean <- data_iso3 %>% 
  group_by(lot_id3) %>%
  dplyr::summarise(across(where(is.numeric), ~ round(mean(.x, na.rm = TRUE),1)))%>%
  melt()

data_iso_sd <- data_iso3 %>% 
  group_by(lot_id3) %>%
  dplyr::summarise(across(where(is.numeric), ~ round(sd(.x, na.rm = TRUE),1)))%>%
  melt()

data_iso_mean_sd=merge(data_iso_mean,data_iso_sd, by=c('lot_id3','variable'))
data_iso_mean_sd$mean_sd=paste(data_iso_mean_sd$variable,'=',data_iso_mean_sd$value.x," ± ",data_iso_mean_sd$value.y)

write.xlsx(data_iso_mean_sd,'who_isotopes/data_iso_mean_sd.xlsx')

#check 

data_iso_w17=data_iso_mean_sd[data_iso_mean_sd$lot_id3 %in% c('W17-W','W17-Y'),]

ggplot(data_iso_w17, aes(x=variable,y=value.x,fill=lot_id3))+
  geom_col(position = "dodge")+
  geom_errorbar(aes(ymin = value.x - value.y, ymax = value.x + value.y),width = .2,
                position = position_dodge(.9))


## hydrogen #### 

data_hydro=read.xlsx("who_isotopes/hydrogen_raw_allreps_edited.xlsx", sheet= "IRMS_Fera")

data_hydro <- data_hydro %>%
  pivot_longer(cols = c(Measurement.1, Measurement.2, Measurement.3), 
               names_to = "tablet_rep", values_to = "hydro")  %>%
  mutate(lot_id3=metadata$lot_id3[match(lot_id,metadata$lot_id)])%>%
  dplyr::select(lot_id3,hydro)%>%
  group_by(lot_id3) %>%
  mutate(lot_id3_rep=paste0(lot_id3, "_", row_number()))

#outliers_hydro=c('W19_3','W02_1','W08_1','W21_1') #max number
outliers_hydro=c('W19_3') #min number

data_hydro2 <- data_hydro %>%
  group_by(lot_id3_rep) %>%
  filter(!lot_id3_rep %in% outliers_hydro) %>% #for the HCAs
  column_to_rownames(var = "lot_id3_rep") %>%
  dplyr::select(-lot_id3)%>%
  na.omit()%>%
  round(2)

data_scaled_hydro<-scale(data_hydro2,center=TRUE,scale=TRUE) 
mean(data_scaled_hydro) #mean 0
sd(data_scaled_hydro)

#Alternatives 
#data_iso2 <- data_iso2 %>%  
#mutate(S34 = replace_na(S34, -20)) # take range of naturally occurring values according to https://tidyr.tidyverse.org/reference/replace_na.html

#to only include C, O and N
#setwd('./only_CNO')

#data3 <- data2 %>%
#select(-S34)

data_scaled_iso=scale(data_iso2)
mean(data_scaled_iso[,1])
sd(data_scaled_iso[,1])

#data_scaled_iso[is.na(data_scaled_iso)] <- 1.9

#### CACO3 ##########

data_caco3=read.xlsx("who_caco3/FORESFA_WHOlogo_CaCO3_and quanti_20250606_edited.xlsx", sheet= "IRMS CaCO3")

data_caco3=data_caco3 %>%
  dplyr::select(lot_id3,XRD,signalmassratio_tab1,signalmassratio_tab2,signalmassratio_tab3,
                c_tab1,c_tab2,c_tab3,o_tab1,o_tab2,o_tab3)%>%
  pivot_longer(cols = c(-lot_id3,-XRD),
  names_to = c("tablet_test",'tablet_rep'),
  names_pattern = "(.*)_(.*)")%>%
  pivot_wider(names_from = tablet_test, values_from = value)%>%
  mutate(rep=gsub('tab','',tablet_rep),
         lot_id3_rep=paste(lot_id3,rep,sep='_'),
         XRD=if_else(XRD=='positive',1,0))#,
         #signalmassratio=if_else(is.na(signalmassratio),0,signalmassratio))

data_caco3_2 <- data_caco3 %>%  
  group_by(lot_id3_rep) %>%
  column_to_rownames(var = "lot_id3_rep") %>%
  dplyr::select(-lot_id3,-tablet_rep,-rep) #same HCA plot is obtained if leaving XRD out


data_scaled_caco3=scale(data_caco3_2)
mean(data_scaled_iso[,1])
sd(data_scaled_iso[,1])

#check what if only using C, O and mass signal
data_caco3_3=na.omit(data_caco3_2) %>%
  dplyr::select(-XRD,-signalmassratio) #W12 associated to W02 instad of W13 if signalmassratio is removed

data_scaled_caco3_3=scale(data_caco3_3)
mean(data_scaled_iso[,1])
sd(data_scaled_iso[,1])

#data_scaled_caco3=na.omit(data_scaled_caco3)

#for table calcium carbonate

data_caco4=data_caco3 %>% 
  dplyr::select(lot_id3,c,o) %>% 
  mutate(lot_id3 = case_when(
  lot_id3 %in% c('W04','W05','W06','W07','W08','W09') ~ "W4-9",
  lot_id3 %in% c("W13", "W14",'W15') ~ "W13-15",
  TRUE ~ lot_id3))

data_caco_mean <- data_caco4 %>% 
  group_by(lot_id3) %>%
  dplyr::summarise(across(where(is.numeric), ~ round(mean(.x, na.rm = TRUE),1)))%>%
  melt()

data_caco_sd <- data_caco4 %>% 
  group_by(lot_id3) %>%
  dplyr::summarise(across(where(is.numeric), ~ round(sd(.x, na.rm = TRUE),1)))%>%
  melt()

data_caco_mean_sd=merge(data_caco_mean,data_caco_sd, by=c('lot_id3','variable'))
data_caco_mean_sd$mean_sd=paste(data_caco_mean_sd$variable,'=',data_caco_mean_sd$value.x," ± ",data_caco_mean_sd$value.y)

write.xlsx(data_caco_mean_sd,'who_caco3/data_caco_mean_sd.xlsx')


### NIR #####

data_nir=read.xlsx("who_NIR/NIR_Data_Clusters.xlsx", sheet= "Sheet1")  %>%
  dplyr::rename(lot_id2=PACS.ID, lot_id3=Code)

data_nir2=data_nir %>%
  mutate(across(-c(lot_id2, lot_id3, Sample.Name, Group), as.numeric))%>%
  group_by(lot_id3,Sample.Name) %>%
  dplyr::summarise(across(everything(), ~ if(is.numeric(.)) mean(.) else NA), .groups = "drop") %>%
  group_by(lot_id3) %>%
  mutate(lot_id3_rep=paste0(lot_id3, "_", row_number())) %>%
  ungroup() %>%
  dplyr::select(lot_id3, lot_id3_rep, everything()) %>%
  dplyr::select(-Measure_number,-Sample.Name,-Group,-lot_id2,-lot_id3) %>%
  column_to_rownames(var = "lot_id3_rep") 

data_scaled_nir=data_nir2


##ICPMS #####

data_icpms=read.xlsx("who_ICPMS/ICPMS_edited.xlsx", sheet= "Sheet1")%>%
  group_by(lot_id3) %>%
  mutate(lot_id3_rep=paste0(lot_id3, "_", row_number()))

data_icpms2<- data_icpms %>%
  dplyr::select(lot_id3, lot_id3_rep, everything()) %>%
  dplyr::select(-Description,-lot_id,-lot_id3) %>%
  column_to_rownames(var = "lot_id3_rep")  %>%
  dplyr::select(-lot_id3)
  
data_icpms2=data_icpms2[,!colnames(data_icpms2) %in% c('lot_id3')] 

icpms_limits=read.xlsx("who_ICPMS/report_limits.xlsx", sheet= "Sheet1")
icpms_limits_half=icpms_limits/2

## using half values
data_icpms3 <- data_icpms2 %>%
  mutate(across(where(is.numeric), function(col) {
    name <- cur_column()
    thr <- icpms_limits[name]
    rep_val <- icpms_limits_half[name]
    col[col < thr] <- rep_val
    as.numeric(col)
  }))


## subsituting limits for 0 -> gives the same as taking the half LOC
#data_icpms2 <- data_icpms %>%
  #mutate(across(where(is.numeric), function(col) {
    #name <- cur_column()
    #thr <- icpms_limits[name]
    #col[col < thr] <- 0
    #col
  #}))



data_scaled_icpms=scale(data_icpms3) 
mean(data_scaled_icpms[,1])
sd(data_scaled_icpms[,1])

data_scaled_icpms <- data_scaled_icpms[, colSums(!is.na(data_scaled_icpms)) > 0]

#remove W03_4 for better grouping in the HCA analyses
data_scaled_icpms=data_scaled_icpms[!rownames(data_scaled_icpms) %in% c('W03_4'),]


## XRD #####

data_xrd=read.xlsx("who_xrd/summary_with_replicates.xlsx", sheet= "Sheet1") %>%
  mutate(lot_id3=metadata$lot_id3[match(lot_id2,metadata$lot_id2)], 
    lot_id3_rep=paste(lot_id3,replicate, sep='_')) %>%
  dplyr::select(lot_id3, lot_id3_rep, everything())

#outliers_xrd=c('W17-Y_1','W17-W_1','W11_1')
outliers_xrd=c('W18_2')

data_xrd2<- data_xrd %>%
  filter(!lot_id3_rep %in% outliers_xrd) %>%
  column_to_rownames(var = "lot_id3_rep") %>%
  dplyr::select(-lot_id3, -replicate, -lot_id2) 

data_xrd2=data_xrd2/rowSums(data_xrd2)

data_scaled_xrd=scale(data_xrd2) 
mean(data_scaled_xrd[,1])
sd(data_scaled_xrd[,1])


#### Parameters for representation #######

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
  else if (x=="11") 'deeppink'
  else if (x=="12") 'lawngreen'
  else if (x=="13") 'cyan'
  else if (x=="14") 'tan1'
  else if (x=="15") 'aquamarine'
  else if (x=="16") 'lemonchiffon2'
  else if (x=="17") 'royalblue1'
  }

#alternative to yellow is turquoise3

colors_hca_groups=c("blue4", "burlywood3", "azure4",'chartreuse3',"blueviolet",'brown','yellow2','red','lightgreen','black',
                    'deeppink','lawngreen','cyan','tan1','aquamarine','lemonchiffon2','royalblue1')


### date colors ####

order_dates=unique(dates$Mfg_date_year[order(dates$Mfg_date_year)])
#color_dates=hue_pal()(9)
color_dates=c('black','grey','grey','grey','grey','seashell2','seashell2','seashell2','lightpink')
names(color_dates) <-order_dates


## HCAs ####

data_input=c('data_scaled_morph','data_scaled_ftir','data_scaled_iso','data_scaled_hydro',
             'data_scaled_caco3','data_scaled_caco3_3','data_scaled_nir','data_scaled_icpms',
             'data_scaled_weight','data_scaled_morph3','data_scaled_xrd')
data_input2=c('data_morph','data_ftir','data_iso','data_hydro',
              'data_caco3','data_caco3','data_nir','data_icpms',
              'data_weight','data_morph3','data_xrd')
directories_output=c('who_morphometry/hca_allmorphology','who_ftir/average','who_isotopes','who_isotopes/hydrogen',
                     'who_caco3','who_caco3/only_co','who_NIR','who_icpms',
                     'who_morphometry/hca_weight','who_morphometry/hca_allmorpho_and_weight','who_xrd')

cut_val_all=c(0.22,10.5,0.94,1.35, #1 for hydro if removing all outliers, 1.35 for minimum outliers
              0.8,0.8,0.38,9.3,
              0,0.365,58) #0.8 for caco3, 42 if leaving outlier for XRD
y_shift_bars=c(-0.75,-5.5,-0.6,-0.3,
               -0.6,-0.6,-0.6,-1.6,
               -0.75,-0.75,-20)

## probably better to not include weight in the morphological parameters

for (i in 1:length(data_input)) {
  print(data_input[i])
  get(data_input[i]) -> data_scaled
  get(data_input2[i]) -> data
  
  #data_scaled=as.matrix(na.omit(data_scaled))

  setwd(paste(path_work,directories_output[i],sep=''))
 
  dist.matrix=dist(data_scaled, method='euclidean')
  hclust_avg <- hclust(dist.matrix, method = 'average')
  
  png('hca1.png',width = 1000, height = 600)
  plot(hclust_avg,cex=0.7,main = "HCA, euclidean dist - average hclust")
  dev.off()  
  
  dend <- as.dendrogram(hclust_avg)
  labels(dend)=gsub('_[0-9]*','',labels(dend))
  labels_cex(dend) <- 0.7
  #labels_colors(dend)=color_countries[match(metadata$collection_country[match(labels(dend),metadata$lot_id3)],names(color_countries))]
  #dend= dend %>% set("leaves_pch", shapes_vector[match(metadata$lot_id3[match(labels(dend),metadata$lot_id3)],names(shapes_vector))])
  #dend=sort(dend)

  png('hca2.png',width = 1000, height = 600)
  plot(dend, main = "HCA, euclidean dist - average hclust", xaxt="n")
  dev.off()
  
  unique(round(hclust_avg$height,4))
  
  cut_val=cut_val_all[i]
  cut_avg <- cutree(hclust_avg, h = cut_val)
  
  # parameters for representation
  
  labels_cutavg=cut_avg
  names(labels_cutavg)=gsub('_[0-9]','',names(labels_cutavg))
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
  
  
  #save dendrogram plot with rectangles groups and bars ###
  
  par(mar=c(5.1, 4.1, 2, 0))
  plot(dend,  ylab='Distance') #xaxt="n" 
  abline(h = cut_val, col = "blue", lty = 2, lwd = 1.5) #plot the abline of the height you think should be correct
  rect.dendrogram(dend, k=max(cut_avg), border = 'grey', xpd=0,cluster = cut_avg) #k=max(cut_avg)
  text(x=labels_position2, y=cut_val, col=color_vector, labels=labels_cutavg, font=2, cex=2)
  colored_bars(colors = dates_contries_bars, dend = dend, sort_by_labels_order = FALSE,rowLabels=c('',''),y_shift=y_shift_bars[i])
  legend("topright", title=expression(bold("Country col.")),
         legend = names(color_countries), 
         fill= color_countries, bty = "n",  pt.cex = 0.5, cex = 1.2, x.intersp = 0.5, 
         text.col = "black", horiz = FALSE, inset = c(0.08, 0),
         title.adj=c(0.1, 0.1))
  legend("topright",
         title=expression(bold("Period mfg.")),
         legend = c('2009','2011-2014','2017-2019','2021'), 
         fill= c('black','grey','seashell2','lightpink'), bty = "n",  pt.cex = 0.5, cex = 1.2, x.intersp = 0.5,
         text.col = "black", horiz = FALSE, inset = c(-0.01,0),
         title.adj=c(0.1, 0.1))
  p0 <- recordPlot()
  dev.off()
  
  png('hca3.png',width = 1000, height = 600)
  p0
  dev.off()
 
   
  
  # create pairwise distance matrix and plot
  pairwise_dist <- as.data.frame(as.table(as.matrix(dist.matrix)))%>%
    rename(distance=Freq) %>%
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
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.8),
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
                  palette = colors_hca_groups, ggtheme =theme_pubr())+ggtitle(NULL)  #ellipse.type='confidence'

  
  if(data_input[i]=='data_scaled_morph') {
  
    data$collection_country=metadata$collection_country[match(data$lot_id3, metadata$lot_id3)]
    data$hca_grouping=as.factor(cut_avg[match(data$lot_id3_rep,names(cut_avg))])
    data$diameter=round(data$diameter,2)
    
    order_lotid=unique(data[order(data$hca_grouping,data$lot_id3),]$lot_id3)
    order_collection=unique(data[order(data$hca_grouping,data$lot_id3),]$collection_country)
    
    p4=ggplot(data, aes(x = factor(lot_id3, levels=order_lotid), y = diameter, color=hca_grouping)) +
      geom_violin()+
      geom_point()+
      scale_color_manual(values=color_vector, name="Clusters")+
      facet_grid(~ factor(collection_country, levels=order_collection), scales='free') + 
      xlab('')+
      ylab('Diameter')+
      scale_y_continuous(breaks = seq(min(data$diameter), max(data$diameter), by = 0.25)) +  # Set x-axis breaks every 0.5
      theme_pubr()+
      theme(
        panel.grid.major = element_line(color = "grey90", size = 0.25),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", size = 0.25), # Minor gridlines
        panel.grid.minor.x = element_line(color = "grey90", size = 0.25), # Minor gridlines on x-axis
        panel.grid.minor.y = element_line(color = "grey90", size = 0.25)  # Minor gridlines on y-axis
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7),
            axis.text.y = element_text(size=10),
            legend.position = 'none')
    
    
    ggsave(p4,filename='scaterplot_allbottles.png', width = 18, height = 10.5)
    
    plots_below=ggarrange(p4,p3,p1,p2,ncol = 4, labels = c('B','',"C", "D"), widths = c(1, 0.6, 0.5, 0.5))
 
  } else {
   plots_below=ggarrange(p3,p1,p2,ncol = 3, labels = c('B',"C", "D"), widths = c(1, 0.6, 0.6))
  } 
 
  all_plots=ggarrange(p0,                                                
                      plots_below +
                        theme(plot.margin = margin(0.1, 2, 0, 0.4, "cm")), 
                      nrow = 2,
                      labels='A', heights = c(1,0.6)) + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))     
  
  if(data_input[i]=='data_scaled_nir'){
  ggsave(all_plots,filename='hca_allplots_average.png', width = 20, height = 12)  
  } else {
  ggsave(all_plots,filename='hca_allplots_average.png', width = 20, height = 10.5)  
  }
  
  
  #create file to save
  
  data$hca_grouping=cut_avg[match(data$lot_id3_rep,names(cut_avg))]
  write.xlsx(data,'summary_hca_groupings.xlsx')
  
  assign(paste('hca_clusters',data_input2[i],sep = '_'),hca_clusters)
  
  setwd(path_work)
  
}


#### #### saving cophenetic distances for later groupings comparisons ####

data_input=c('data_scaled_morph','data_scaled_ftir','data_scaled_iso','data_scaled_caco3','data_scaled_nir','data_scaled_icpms')
name_dist=c('morph_dist','ftir_dist','iso_dist','caco3_dist','nir_dist','icpms_dist')

pairwise_cophdist_all=list()

for (i in 1:length(data_input)) {
  get(data_input[i]) -> data_scaled

  data_scaled=na.omit(data_scaled)
  
  dist.matrix=dist(data_scaled, method='euclidean')
  hclust_avg <- hclust(dist.matrix, method = 'average')

  coph.matrix=cophenetic(hclust_avg)
  
  pairwise_cophdist_all[[i]] <- as.data.frame(as.table(as.matrix(coph.matrix)))%>%
    filter(Var1 != Var2) %>%
    mutate(Freq=rescale(Freq))
      
  pairwise_cophdist_all[[i]]$unique_combination=apply(pairwise_cophdist_all[[i]],1, function(x) {paste(sort(c(x['Var1'], x['Var2'])), collapse = "-")})
    
  pairwise_cophdist_all[[i]]<-pairwise_cophdist_all[[i]]%>%
  dplyr::select(unique_combination,Freq)%>%
  distinct()
  
  names(pairwise_cophdist_all[[i]])[names(pairwise_cophdist_all[[i]])=='Freq'] <- name_dist[i]
  

}

pairwise_cophdist_all2 <- Reduce(function(x, y) merge(x, y, by = "unique_combination"), pairwise_cophdist_all)  

write.xlsx(pairwise_cophdist_all2,'cophdist_all_variables.xlsx')


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



#loop

for (i in 1:length(data_input)) {
  get(data_input[i]) -> data_scaled
  get(data_input2[i]) -> data
  
  setwd(paste(path_work,directories_output[i],sep=''))
  
  dist.matrix=dist(data_scaled, method='euclidean')
  
  pairwise_dist <- as.data.frame(as.table(as.matrix(dist.matrix)))%>%
    dplyr::rename(distance=Freq) %>%
    filter(Var1 != Var2) %>%
    mutate(Var1=gsub('_[0-9]','',Var1),
           Var2=gsub('_[0-9]','',Var2),
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
  
  pcoa_df$lot_id3=gsub('_[0-9]*','',rownames(pcoa_df))
  pcoa_df$source=metadata$source[match(pcoa_df$lot_id3,metadata$lot_id3)]
  pcoa_df$collection_country=metadata$collection_country[match(pcoa_df$lot_id3,metadata$lot_id3)]
  pcoa_df$mfg_year=as.factor(dates$Mfg_date_year[match(pcoa_df$lot_id3,metadata$lot_id3)])
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
  
  ggsave(  
    free(rep_time,type = "label") + free(pcoa_time,type = "label") + free(rep_place,type = "label") + free(pcoa_place,type = "label") +
      plot_layout(ncol = 2) +
      plot_layout(widths = c(3, 1,3,1))+
      plot_layout(guides = "collect")  & guides(color = guide_legend(ncol = 3)), 
    filename='time_pcoa_plots.png', width = 20, height = 10.5)
  
  
    setwd(path_work)
  
}


### look at the loads of the PCA for the ICPMS data -> LDA

hca_clusters_data_icpms

colnames(data_scaled_icpms)=gsub("^[^.]*\\.",'',colnames(data_scaled_icpms))

res.pca<- prcomp(data_scaled_icpms,  scale = FALSE) #PCoA using distances rather than raw data, labelling the groups


pca_ind=fviz_pca_ind(res.pca,
             habillage=hca_clusters_data_icpms,
             palette = colors_hca_groups,
             repel = TRUE     # Avoid text overlapping
)

pca_var=fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             geom='text'
)

ggarrange(pca_ind,pca_var)

fviz_pca_biplot(res.pca,
                geom.var=c('text'),
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE)


pca_var_data=pca_var$data

pca_var_data=pca_var_data[order(-pca_var_data$contrib),]

pca_var_data



## other ----

# xrd iden ####

xrd_ident=read.xlsx('xrd_identified.xlsx') %>%
  tibble::column_to_rownames('lot_id3')


dist.matrix=dist(xrd_ident, method='binary')
hclust_avg <- hclust(dist.matrix, method = 'single')

plot(as.dendrogram(hclust_avg))

xrd_ident$hca_grouping=cut_avg[match(data$lot_id3_rep,names(cut_avg))]
write.xlsx(data,'summary_hca_groupings.xlsx')


#others

#xrd2=api2=read.xlsx('../who_xrd/others/who_xrd_results_edited.xlsx') %>%
#mutate(lot_id3=metadata$lot_id3[match(lot_id,metadata$lot_id)])

#write.xlsx(xrd2,'../who_xrd/others/who_xrd_results_edited2.xlsx')


#api2=read.xlsx('../who_orbitrap/WHO_orbitrap.xlsx')%>%
#mutate(lot_id2=gsub(".*-", "", lot_id2),
#lot_id3=metadata$lot_id3[match(lot_id2,metadata$lot_id2)])

#write.xlsx(api2,'../who_orbitrap/WHO_orbitrap_edited.xlsx')

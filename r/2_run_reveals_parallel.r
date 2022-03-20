library(raster)
library(sp)
#library(DISQOVER)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
#library(neotoma)
library(maptools)
library(fields)
library(maps)
library(rgdal)
library(truncnorm)

#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)

# # load package
# load('DISQOVER/R/sysdata.rda')

source('r/utils/main.R')

version="8.0"
LC6k=FALSE

args = commandArgs(trailingOnly=TRUE)
id_to_run = as.integer(args[1])

##################################################################################################################################################
## read in pollen data for NA
##################################################################################################################################################

pollen_trans = readRDS(paste0('data/pollen_bin_v', version, '.RDS'))
pollen_trans = pollen_trans[, -ncol(pollen_trans)]
taxa = colnames(pollen_trans)[9:ncol(pollen_trans)]

# taxa_list = data.frame(latin=taxa)
# write.csv(taxa_list, 'data/taxa_latin2common.csv', row.names=FALSE)

colnames(pollen_trans)[which(colnames(pollen_trans) == 'age_calbp')] = 'age'

# make time bins
if (LC6k){
# # start with: 0.1-0.35k BP, 0.35-0.7k BP, 5.5-6.2k BP
# breaks = c(0.1, 0.35, 0.7, 2.7, 3.2, 5.7, 6.2)
# slice_bins = seq(1, 6)
# slice_labels = c(50, 200, 500, 1500, 3000, 6000)

breaks = c(-74, 0.1, 0.35, 0.7, seq(1.2, 11.7, by=0.5))
slice_bins = seq(1, 25)
slice_labels = c(50, 200, seq(500, 11500, by=500))

# cuts = c(5.7, 6.2)
} else {
  # 2000 year intervals
  breaks = seq(0, 22000, by=2000)#c(0, 0.35, 0.7, 2.7, 3.2, 5.7, 6.2)
  slice_bins = seq(1, 11)
  slice_labels = seq(1000, 21000, by=2000)#c(50, 200, 500, 1500, 3000, 6000)
}

pollen_trans$slice_bin = cut(pollen_trans$age, breaks*1000, labels=FALSE)
colnames(pollen_trans) = tolower(colnames(pollen_trans))
pollen_trans = pollen_trans[which(!is.na(pollen_trans$slice_bin)), ]

# sum samples within a time bin for each site 
pollen_bin = ddply(pollen_trans, .(dataset_id, lat, long, sitename, lake_size, slice_bin),  
                   function(x) colSums(x[tolower(taxa)]))

# make grid for NA (or ENA)
source('r/make_grid.R')
grid <- make_grid(pollen_bin, coord_fun = ~ long + lat, projection = '+init=epsg:4326', resolution = 1)
  
cell_id <- raster::extract(grid, pollen_bin[,c('long', 'lat')])

pollen_bin <- data.frame(cell_id, pollen_bin)

# # remove other classes for now
# pollen_bin <- pollen_bin[,which(!(colnames(pollen_bin) %in% c('other.hardwood', 'other.conifer')))]
# taxa <- taxa[which(!(taxa %in% c('Other hardwood', 'Other conifer')))]

##################################################################################################################################################
## fill in missing lake sizes
##################################################################################################################################################

# library(fitdistrplus)

lake_sizes = as.numeric(unique(pollen_bin$lake_size))
lake_sizes = lake_sizes[which(!is.na(lake_sizes))]
lake_sizes = lake_sizes[which(lake_sizes<20000)]
lake_sizes = data.frame(area=lake_sizes)

# # hist(lake_sizes)
# # 
# # lake_sizes = unique(pollen_bin[which(pollen_bin$lake_size<10000),'lake_size'])
# # hist(lake_size, breaks=30)
# # 
# fit = fitdist(lake_sizes$area/100, 'exp')
# 
# lake_dist = rexp(1000, rate=fit$estimate[1])
# # lake_dist = rnorm(1000, mean=fit$estimate[1], sd=fit$estimate[2])
# # lake_dist = rweibull(1000, shape=fit$estimate[1], scale=fit$estimate[2])
# lake_dist = data.frame(area=lake_dist)
# 
# ggplot() +
#   geom_histogram(data=lake_dist, aes(x=area*100, y=..density..)) +
#   geom_histogram(data=lake_sizes, aes(x=area, y=..density..), fill='lightblue', alpha=0.5) 
# 
# 

# use 50 hectares are medium lake size
# assigned to all lakes with missing lake size
medium_radius = sqrt(50*0.01 / pi)*1000

idx_na = which(is.na(pollen_bin$lake_size))

# pollen_bin$lake_size[idx_na] = rexp(length(idx_na), rate=fit$estimate[1])*100
pollen_bin$lake_size[idx_na] = medium_radius


# 
# # lake sizes
# lake_sizes = read.csv('data/areas_with_datasetID.csv')
# lake_sizes = lake_sizes[!duplicated(lake_sizes$DatasetID),] 
# 
# lake_size = lake_sizes$AREAHA[match(pollen_bin$dataset, lake_sizes$DatasetID)]
# lake_size[which(is.na(lake_size))] = lake_sizes$Area[match(pollen_bin$dataset[which(is.na(lake_size))], lake_sizes$DatasetID)]
# 
# 
# random_lake_size <- function(lake_size){
#   library(truncnorm)
#   lake_radius = rtruncnorm(1,
#                            a=5,
#                            mean=mean(lake_size, na.rm=TRUE),
#                            sd=sd(lake_size, na.rm=TRUE))
#   
#   return(lake_radius)
# }

# in HA; convert to radius in m
# pi * r * r
# lake_size = sqrt(lake_size*0.01 / pi)*1000
# pollen_bin <- data.frame(lake_size, pollen_bin)

# pollen_bin$lake_size = sqrt(pollen_bin$lake_size*0.01 / pi)*1000


##################################################################################################################################################
## sort some taxa stuff
##################################################################################################################################################
taxa_trans = read.csv('data/taxa_latin2common_filled.csv')

colnames(pollen_bin)[8:ncol(pollen_bin)] = taxa_trans$common[match(colnames(pollen_bin[,8:ncol(pollen_bin)]),taxa_trans$latin)]
colnames(pollen_bin)[is.na(colnames(pollen_bin))] <- 'unknown'

# pollen_bin = pollen_bin %>%
#   pivot_longer(cols=colnames(pollen_bin)[8:ncol(pollen_bin)]) %>%
#   group_by(cell_id, dataset_id, lat, long, sitename, lake_size, slice_bin, name) %>% 
#   summarize(value=sum(value), .groups='keep') %>%
#   pivot_wider(names_from = name, values_from = value)


pollen_bin = pollen_bin %>%
  pivot_longer(cols=colnames(pollen_bin)[8:ncol(pollen_bin)]) %>%
  group_by(cell_id, dataset_id, lat, long, sitename, lake_size, slice_bin, name) %>% 
  summarize(value=sum(value), .groups='keep') %>%
  pivot_wider(names_from = name, values_from = value)

# # sum repeated columns
# dat.dup <- foo
# dups <- unique(names(dat.dup)[duplicated(names(dat.dup))])
# for (i in dups) {
#   dat.dup[[i]] <- rowSums(dat.dup[names(dat.dup) == i])
# }
# pollen_bin <- dat.dup[!duplicated(names(dat.dup))]
# 
# 
matched = taxa_trans$common[match(colnames(pollen_bin[,8:ncol(pollen_bin)]),taxa_trans$common)]
taxa_common = sort(matched[which(!is.na(matched))])

pollen_bin = data.frame(pollen_bin[,1:7], pollen_bin[,taxa_common])

# pollen_bin[,8:ncol(pollen_bin)] = pollen_bin[,taxa_common]
# 
# #
# colnames(pollen_bin)[8:ncol(pollen_bin)] = taxa_trans$common[match(colnames(pollen_bin[,8:ncol(pollen_bin)]),taxa_trans$latin)]
# colnames(pollen_bin)[is.na(colnames(pollen_bin))] <- 'unknown'
# 
# # sum repeated columns
# dat.dup <- pollen_bin
# dups <- unique(names(dat.dup)[duplicated(names(dat.dup))])
# for (i in dups) {
#   dat.dup[[i]] <- rowSums(dat.dup[names(dat.dup) == i])
# }
# pollen_bin <- dat.dup[!duplicated(names(dat.dup))]
# 
# 
# matched = taxa_trans$common[match(colnames(pollen_bin[,8:ncol(pollen_bin)]),taxa_trans$common)]
# taxa_common = sort(matched[which(!is.na(matched))])
# pollen_bin[,8:ncol(pollen_bin)] = pollen_bin[,taxa_common]
##################################################################################################################################################
## read in and prep ppes and svs
##################################################################################################################################################

N_taxa = length(taxa_common)
pars = data.frame(species = taxa_common, 
                  PPEs = numeric(N_taxa), 
                  PPE.errors = numeric(N_taxa), 
                  fallspeed = numeric(N_taxa)) 

ppes_NA = readRDS('data/PPEs_grass.RDS')

idx_match = match(pars$species, tolower(ppes_NA$taxon))
pars$PPEs = ppes_NA[idx_match, 'ppe']
pars$PPE.errors = ppes_NA[idx_match, 'error']

svs_NA = read.csv('data/svs_LC6K.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
svs_NA$taxon = tolower(svs_NA$taxon)

# svs_agg = aggregate(sv ~ taxon, svs, median, na.rm=TRUE)
# 
# svs_agg = svs_agg[match(taxa_common, svs_agg$taxon),]
# svs_agg[which(is.na(svs_agg$sv)), 'sv'] = 0.01

 
idx_match = match(pars$species, svs_NA$taxon)
pars$fallspeed = svs_NA[idx_match, 'sv']

ppes_NH = read.csv('data/PPE/RPP_FS_wieczorek.csv', stringsAsFactors=FALSE)

for (k in 1:N_taxa){
  taxon = taxa_common[k]
  
  idx_row = which(pars$species == taxon)
  
  if (is.na(pars$PPEs[idx_row])){
    idx_match = match(taxon, tolower(ppes_NH$taxon))
    
    if ((!is.na(ppes_NH[idx_match, "RPP.v2..America."])) & (!is.na(ppes_NH[idx_match, "SE..America."]))){
      pars[idx_row, 'PPEs'] = ppes_NH[idx_match, "RPP.v2..America."]
      pars[idx_row, 'PPE.errors'] = ppes_NH[idx_match,  "SE..America." ]
    } else {
      pars[idx_row, 'PPEs'] = ppes_NH[idx_match, "ppe"]
      pars[idx_row, 'PPE.errors'] = ppes_NH[idx_match,  "error" ]
    }
    
  }
  
  if (is.na(pars$fallspeed[idx_row])){
    idx_match = match(taxon, tolower(ppes_NH$taxon))
    
    if (taxon == 'pea') {
      idx_match = match('legume', tolower(ppes_NH$taxon))
    }
    
    if (!is.na(ppes_NH[idx_match,  "FS..m.s...America."])){
      pars[idx_row, 'fallspeed'] = ppes_NH[idx_match,  "FS..m.s...America."]
    } else {
      pars[idx_row, 'fallspeed'] = ppes_NH[idx_match, "fallspeed"]
      
    }
    
  }
  
}


pars$species = tolower(pars$species)

all(match(taxa_common, pars$species))

# 
# ppes$ppe = as.numeric(ppes$ppe)
# ppes$error = as.numeric(ppes$error)
# ppes$fallspeed = as.numeric(ppes$fallspeed)
# 
# ppes = ppes[,c('taxon', 'ppe', 'error', 'fallspeed')]
# ppes = rbind(ppes,
#              c('Ambrosia', 3.52, 0.81, 0.013))
# # missing fallspeed
# # ppes = rbind(ppes,
# #              c('Dogwood', 1.70, 0))
# ppes = rbind(ppes,
#              c('Hemlock', 2.31, 4.26, 0.08))

# ppes_old = readRDS('data/PPEs_agg.RDS')
# 
# 
# 
# hemlock_error <- ppes_old$ppe[ppes_old$taxon=='Hemlock']/ppes_old$error[ppes_old$taxon=='Hemlock']
# grass_hemlock_ratio <-  ppes_old$ppe[ppes_old$taxon=='Grass']/ppes_old$ppe[ppes_old$taxon=='Hemlock']
# hemlock_ppe <- ppes$ppe[ppes$taxon=='Grass']/grass_hemlock_ratio

# ppes$taxon = tolower(ppes$taxon)
# ppes = ppes[which(ppes$taxon %in% taxa_common),]
# 
# # ppes[which(ppes$error == 0),'error'] = 0.01
# 
# ppes = ppes[match(taxa_common, ppes$taxon),]

# 
# svs_agg = aggregate(sv ~ taxon, svs, median, na.rm=TRUE)
# 
# svs_agg = svs_agg[match(taxa_common, svs_agg$taxon),]
# svs_agg[which(is.na(svs_agg$sv)), 'sv'] = 0.01

# reveals_inputs = data.frame(species=taxa_common,
#                             fallspeed=ppes$fallspeed,
#                             PPEs=ppes$ppe,
#                             PPE.errors=ppes$error)
# rownames(reveals_inputs) = NULL

reveals_inputs = pars
rownames(reveals_inputs) = NULL


# # construct data frame of reveals inputs
# reveals_inputs = data.frame(species=taxa, 
#                             fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa)], 
#                             PPEs=ppe_stepps_stat[,2], 
#                             PPE.errors=ppe_stepps_stat[,3])
# rownames(reveals_inputs) = NULL
write.csv(reveals_inputs, "data/reveals_input/params.csv", row.names=FALSE)
# ##################################################################################################################################################
# ## read in and prep ppes and svs
# ##################################################################################################################################################
# 
# ppes = readRDS('data/PPEs_agg.RDS')
# 
# # taxa[which(!(taxa %in% ppes$taxon))]
# ppes$taxon = tolower(ppes$taxon)
# ppes = ppes[which(ppes$taxon %in% taxa_common),]
# 
# ppes[which(ppes$error == 0),'error'] = 0.01
# 
# ppes = ppes[match(taxa_common, ppes$taxon),]
# 
# # ###############################################################################################################################
# # # read in PPEs
# # ppes = read.csv('data/ppes_v2.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
# # # ppes$taxon[which(ppes$taxon == "OTHER.CONIFER")]  = "FIR"
# # # ppes$taxon[which(ppes$taxon == "OTHER.HARDWOOD")] = "OTHER HARDWOOD"
# # 
# # # think carefully about how to rescale these
# # maple_rescale = mean(ppes[which((ppes$taxon=='MAPLE')&(ppes$continent == 'North America')&(ppes$tag != 'calcote')),'ppe'])
# # ppes[which(ppes$tag == 'calcote'),'ppe'] = ppes[which(ppes$tag == 'calcote'),'ppe']*0.61
# # 
# # ppes_agg = aggregate(ppe ~ taxon + continent, ppes, function(x) quantile(x, c(0.05, 0.5, 0.95)))
# # ppes_sd = aggregate(ppe ~ taxon + continent, ppes, function(x) sd(x))
# # ppes_sd[which(is.na(ppes_sd$ppe)), 'ppe'] = 0.01
# # 
# # ppes_agg_NA = subset(ppes_agg, continent == 'North America')
# # 
# # ggplot(data=ppes_agg_NA) + geom_point(aes(x=ppe[,2], y=reorder(taxon, ppe[,2])))
# # 
# # taxa = unique(ppes_agg_NA$taxon)
# # 
# # # # read in STEPPS PPEs
# # # ppe_stepps = readRDS('data/ppe_post.RDS')
# # # # ppe_stepps = readRDS('data/ppe_stepps.RDS')
# # # ppe_stepps_ref = t(apply(ppe_stepps, 1, function(x) x/x[10]))
# # # ppe_stepps_stat  = data.frame(taxa, t(apply(ppe_stepps_ref, 2, function(x) c(mean(x), sd(x)))))
# 
# ###############################################################################################################################
# ## read in SVs
# svs = read.csv('data/svs_LC6K.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
# svs$taxon = tolower(svs$taxon)
# 
# svs_agg = aggregate(sv ~ taxon, svs, median, na.rm=TRUE)
# 
# 
# svs_agg = svs_agg[match(taxa_common, svs_agg$taxon),]
# svs_agg[which(is.na(svs_agg$sv)), 'sv'] = 0.01
# 
# # svs_agg = svs_agg[which(svs_agg$taxon %in% taxa_common), ]
# 
# ## construct param input data frame
# # reveals_inputs = data.frame(species=taxa_common,
# #                             fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa_common)],
# #                             PPEs=ppes$ppe_scaled[match(ppes$taxon, taxa_common)],
# #                             PPE.errors=ppes$error[match(ppes$taxon, taxa_common)])
# 
# reveals_inputs = data.frame(species=taxa_common,
#                             fallspeed=svs_agg$sv,
#                             PPEs=ppes$ppe,
#                             PPE.errors=ppes$error)
# rownames(reveals_inputs) = NULL
# 
# 
# # # construct data frame of reveals inputs
# # reveals_inputs = data.frame(species=taxa, 
# #                             fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa)], 
# #                             PPEs=ppe_stepps_stat[,2], 
# #                             PPE.errors=ppe_stepps_stat[,3])
# # rownames(reveals_inputs) = NULL
# write.csv(reveals_inputs, "data/reveals_input/params.csv", row.names=FALSE)

##################################################################################################################################################
## run reveals
##################################################################################################################################################
# # start with only one slice
# slice = 6
# # pol_dat = pollen_bin[which(pollen_bin$slice_bin == 6),]
# ids     = unique(pol_dat$dataset)

# pollen_bin = pollen_bin[which(!is.na(pollen_bin$slice_bin)),]

ids = unique(pollen_bin$dataset)
pol_dat = pollen_bin

print(length(ids))


veg_pred = data.frame(id=numeric(0), 
                      long=numeric(0), 
                      lat=numeric(0), 
                      cell_id = numeric(0),
                      basin_radius=numeric(0),
                      ages=numeric(0), 
                      taxon=character(0), 
                      meansim=numeric(0), 
                      mediansim=numeric(0), 
                      q10sim=numeric(0), 
                      q90sim=numeric(0), 
                      sdsim=numeric(0))



  
  id = ids[id_to_run]
  
  counts_site = data.frame(ages=slice_labels[pol_dat[which(pol_dat$dataset_id == id),'slice_bin']], 
                           pol_dat[which(pol_dat$dataset_id == id),which(colnames(pol_dat) %in% tolower(taxa_common))])
  colnames(counts_site) = c('ages', taxa_common)
  # counts_site = counts_site[which(rowSums(counts_site[,2:ncol(counts_site)])!=0),]
  
  basin_radius = pol_dat[which(pol_dat$dataset_id == id), c('lake_size')][1]
  # if (is.na(basin_radius)){
  #   basin_radius = random_lake_size(lake_size)
  # }
  
  coords_site = data.frame(pol_dat[which(pol_dat$dataset_id == id), c('dataset_id', 'long', 'lat', 'cell_id')][1,], basin_radius)
  rownames(coords_site) = NULL
  
  # prob don't need this step - fix later
  write.csv(counts_site, 'data/reveals_input/reveals_input.csv', row.names=FALSE)

  # 
  # # cycle through and estimate background veg
  # # # with english csv files
  # a <- REVEALSinR(pollenFile = counts_site,
  #                 pf         = reveals_inputs,
  #                 filetype   = "object",
  #                 dwm        = "gpm neutral", #"LSM unstable",
  #                 tBasin     = "lake",
  #                 dBasin     = 2*round(lake_rad_site), # diameter!
  #                 regionCutoff = 100000,
  #                 repeats      = 1000)


  
  # cycle through and estimate background veg
  # with english csv files
  a <- REVEALSinR(pollenFile = "data/reveals_input/reveals_input.csv",
                  pf         = "data/reveals_input/params.csv",
                  filetype   = "csv",
                  dwm        = "gpm neutral",
                  tBasin     = "lake",
                  dBasin     = 2*round(basin_radius), # diameter!
                  regionCutoff = 100000,
                  repeats      = 1000)
  
  veg = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
  # veg$variable = gsub('OTHER.', 'OTHER_', veg$variable)
  
  veg$taxon = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
  veg$type  = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))
  
  veg_cast = dcast(veg, ages + taxon ~ type)
  
  # ggplot(data=veg_cast) + 
  #   geom_line(aes(x=ages, y=mediansim, colour=taxon)) + 
  #   geom_point(aes(x=ages, y=mediansim, colour=taxon)) + 
  #   geom_errorbar(aes(x=ages, ymin=q10sim, ymax=q90sim, colour=taxon), width=1.5) +
  #   theme_bw()
  
  veg_pred = rbind(veg_pred, data.frame(coords_site, veg_cast))
  
if (LC6k){
  fname = sprintf('data/cache/veg_pred_LC6k_%06d.RDS', id_to_run)
} else {
  fname = sprintf('data/cache/veg_pred_LGM_%06d.RDS', id_to_run)
}
saveRDS(veg_pred, fname)

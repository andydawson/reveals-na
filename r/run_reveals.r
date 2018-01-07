library(raster)
library(sp)
library(DISQOVER)
library(reshape2)
library(ggplot2)
library(dplyr)
library(neotoma)
library(maptools)
library(fields)
library(maps)
library(rgdal)

devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)

# # load package
load('DISQOVER/R/sysdata.rda')

# source('r/utils/main.R')

ena <- TRUE

latlimits  <- c(40, 80) 
longlimits <- c(-165, -50) 
bounding_box = c(-165, 40, -50, 80)

na_shp <- readOGR('data/map_data/world/ne_50m_admin_1_states_provinces_lakes.shp')
proj4string(na_shp)
na_fort=na_shp

##################################################################################################################################################
## pull pollen data for NA
##################################################################################################################################################

# if(!'na_downloads.rds' %in% list.files('data/cache/')) {
#   canada <- get_dataset(gpid='Canada', datasettype = 'pollen') %>% get_download
#   usa    <- get_dataset(gpid='United States', datasettype = 'pollen') %>% get_download
#   na_pollen <- neotoma::bind(canada, usa)
# 
#   saveRDS(na_pollen, file = 'data/cache/na_downloads.rds')
# } else {
#   na_pollen <- readRDS('data/cache/na_downloads.rds')
# }

if(!'ena_downloads.rds' %in% list.files('data/cache/')) {
  ena_pollen <- get_dataset(loc=bounding_box, datasettype = 'pollen') %>% get_download
  saveRDS(ena_pollen, file = 'data/cache/ena_downloads.rds')
} else {
  ena_pollen <- readRDS('data/cache/ena_downloads.rds')
}

# ena_taxa <- lapply(taxa(ena_pollen, collapse = FALSE), as.data.frame) %>% bind_rows

# generate_tables(ena_pollen, output = 'things.csv')

if (ena) {
  datasets = ena_pollen
} else {
  datasets = na_pollen
}

# calib_dialect <- pull_timebins(downloads, calib_range = c(150, 350))
pollen_dialect <- compile_downloads(ena_pollen)

# # compile the pollen taxa
# generate_tables(na_pollen, output = 'data/pol_trans.csv')

# for now
pol_trans_edited <- read.csv('data/pol_trans_edited_LC6K.csv', stringsAsFactors=FALSE)
taxa = sort(unique(pol_trans_edited[which(!is.na(pol_trans_edited$match)),'match']))
taxa = taxa[which(!(taxa %in% c('Other conifer', 'Other hardwood')))]

# translate taxa
pollen_trans <- translate_taxa(pollen_dialect, 
                              pol_trans_edited,
                              id_cols = colnames(pollen_dialect)[1:10])
colnames(pollen_trans) = tolower(colnames(pollen_trans))
# colnames(pollen_trans)[which(colnames(pollen_trans)=='other hardwood')] = 'other.hardwood'
# colnames(pollen_trans)[which(colnames(pollen_trans)=='other conifer')] = 'other.conifer'

# calibrate radiocarbon years BP ages to calendar years BP
library(Bchron)

pollen_trans$age_calBP = pollen_trans$age
idx_rc = which((pollen_trans$date.type == 'Radiocarbon years BP') & (pollen_trans$age >= 71)& (pollen_trans$age <= 46401))
print(paste0('Calibrating all ', length(idx_rc), ' radiocarbon ages'))

for (i in 1:length(idx_rc)){
  print(i)
  ages = BchronCalibrate(ages   = pollen_trans[idx_rc[i], 'age'], 
                        ageSds = rep(200, length(idx_rc[i])), 
                        calCurves = rep('intcal13', length(idx_rc[i])))
  age_samples = SampleAges(ages)
  pollen_trans[idx_rc[i],'age_calBP'] = apply(age_samples, 2, quantile, prob=c(0.5))
  
}

# make time bins
# start with: 0.1-0.35k BP, 0.35-0.7k BP, 5.5-6.2k BP
breaks = c(-74, 0.35, 0.7, 2.7, 3.2, 5.5, 6.5)
slice_bins = seq(1, 6)
slice_labels = c(50, 200, 500, 1500, 3000, 6000)
# cuts = c(5.7, 6.2)

pollen_trans$slice_bin = cut(pollen_trans$age_calBP, breaks*1000, labels=FALSE)
colnames(pollen_trans) = tolower(colnames(pollen_trans))
pollen_trans = pollen_trans[which(!is.na(pollen_trans$slice_bin)), ]

# sum samples within a time bin for each site 
library(plyr)
pollen_bin = ddply(pollen_trans, c('dataset', 'lat', 'long', 'site.name', 'slice_bin'),  
                   function(x) colSums(x[tolower(taxa)]))

# make grid for NA (or ENA)
source('r/make_grid.R')
grid <- make_grid(pollen_bin, coord_fun = ~ long + lat, projection = '+init=epsg:4326', resolution = 1)
  
cell_id <- extract(grid, pollen_bin[,c('long', 'lat')])

pollen_bin <- data.frame(cell_id, pollen_bin)

# # remove other classes for now
# pollen_bin <- pollen_bin[,which(!(colnames(pollen_bin) %in% c('other.hardwood', 'other.conifer')))]
# taxa <- taxa[which(!(taxa %in% c('Other hardwood', 'Other conifer')))]

##################################################################################################################################################
## read in and prep lake sizes
##################################################################################################################################################
# lake_sizes = read.csv('data/data_lakes_20170807.csv')
# lake_sizes = read.csv('data/datasets_aug132017.csv')
# lake_sizes = read.csv('data/andrea_areas.csv')
lake_sizes_stid = read.csv('data/area_lakes_1.1.csv')

lake_sizes = data.frame(dsid=numeric(0), 
                        site.name=character(0), 
                        lat=numeric(0), 
                        long=numeric(0), 
                        area=numeric(0), 
                        type=character(0))    
for (i in 1:nrow(lake_sizes_stid)){  
  pol.datasets = get_dataset(lake_sizes_stid$stid[i], datasettype="pollen")
  pol.temp = data.frame(dsid=names(pol.datasets), lake_sizes_stid[i,][rep(1, length(pol.datasets)),2:ncol(lake_sizes_stid)])
  lake_sizes=rbind(lake_sizes, pol.temp)
}
  
pollen_bin <- data.frame(lake_size=rep(NA), pollen_bin)

pollen_bin$lake_size = lake_sizes$area[match(pollen_bin$dataset, lake_sizes$dsid)]
pollen_bin$lake_size_bool = !is.na(pollen_bin$lake_size)

# missing area
length(unique(pollen_bin$dataset[which(is.na(pollen_bin$lake_size))]))
# with area
length(unique(pollen_bin$dataset[which(!is.na(pollen_bin$lake_size))]))
# zero area
length(unique(pollen_bin$dataset[which(pollen_bin$lake_size==0)]))

p <- ggplot(data=pollen_bin)
p <- p + geom_point(aes(x=long, y=lat, colour=lake_size_bool))
p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + coord_cartesian(xlim = longlimits, ylim = latlimits)
print(p)
ggsave('figures/lake_size_bool_map.pdf')
# 
# # chron control
# # chron <- read.csv('data/chron_control_status.csv')
# # chron = data.frame(ncontrol=rowSums(lake_sizes[,6:ncol(chron)], na.rm=TRUE), chron)
# 
# # which ones match, have no lake size, and have 3 or more controls
# pollen_bin$ncontrol = lake_sizes$controls[match(pollen_bin$dataset, lake_sizes$dsid)]
# 
# match(pollen_bin$dataset[which((pollen_bin$lake_size_bool==FALSE)&(pollen_bin$ncontrol > 2))], lake_sizes$dsid)
# 
# length(unique(lake_sizes$dsid[match(pollen_bin$dataset[which((pollen_bin$lake_size_bool==FALSE)&(pollen_bin$ncontrol > 2))], 
#                                     lake_sizes$dsid)]))
# 
# saveRDS(unique(pollen_bin$dataset[which((pollen_bin$lake_size_bool==FALSE)&(pollen_bin$ncontrol > 2))]), 'data/ids_missing_lake_size.RDS')

# # lake sizes from neo as well
# lake_sizes_all = read.csv('data/areas_with_datasetID.csv')
# lake_sizes_all = lake_sizes_all[!duplicated(lake_sizes_all$DatasetID),] 
# 
# idx_na = which(is.na(pollen_bin$lake_size))
# # pollen_bin$lake_size[idx_na] = lake_sizes$AREAHA[match(pollen_bin$dataset[idx_na], lake_sizes_all$DatasetID)]
# lake_size = lake_sizes_all$AREAHA[match(pollen_bin$dataset[idx_na], lake_sizes_all$DatasetID)]
# lake_size[which(is.na(lake_size))] = lake_sizes$Area[match(pollen_bin$dataset[which(is.na(lake_size))], lake_sizes$DatasetID)]


random_lake_size <- function(lake_size){
  library(truncnorm)
  lake_radius = rtruncnorm(1,
                           a=5,
                           mean=mean(lake_size, na.rm=TRUE),
                           sd=sd(lake_size, na.rm=TRUE))
  
  return(lake_radius)
}

# in HA; convert to radius in m
# pi * r * r
pollen_bin$lake_size = sqrt(pollen_bin$lake_size*0.01 / pi)*1000

saveRDS(pollen_bin, 'data/pollen_bin.RDS')

##################################################################################################################################################
## read in and prep ppes and svs
##################################################################################################################################################
ppes = readRDS('data/PPEs_agg.RDS')
ppes = ppes[which(ppes$taxon %in% taxa),]

ppes[which(ppes$error == 0),'error'] = 0.01

# ###############################################################################################################################
# # read in PPEs
# ppes = read.csv('data/ppes_v2.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
# # ppes$taxon[which(ppes$taxon == "OTHER.CONIFER")]  = "FIR"
# # ppes$taxon[which(ppes$taxon == "OTHER.HARDWOOD")] = "OTHER HARDWOOD"
# 
# # think carefully about how to rescale these
# maple_rescale = mean(ppes[which((ppes$taxon=='MAPLE')&(ppes$continent == 'North America')&(ppes$tag != 'calcote')),'ppe'])
# ppes[which(ppes$tag == 'calcote'),'ppe'] = ppes[which(ppes$tag == 'calcote'),'ppe']*0.61
# 
# ppes_agg = aggregate(ppe ~ taxon + continent, ppes, function(x) quantile(x, c(0.05, 0.5, 0.95)))
# ppes_sd = aggregate(ppe ~ taxon + continent, ppes, function(x) sd(x))
# ppes_sd[which(is.na(ppes_sd$ppe)), 'ppe'] = 0.01
# 
# ppes_agg_NA = subset(ppes_agg, continent == 'North America')
# 
# ggplot(data=ppes_agg_NA) + geom_point(aes(x=ppe[,2], y=reorder(taxon, ppe[,2])))
# 
# taxa = unique(ppes_agg_NA$taxon)
# 
# # # read in STEPPS PPEs
# # ppe_stepps = readRDS('data/ppe_post.RDS')
# # # ppe_stepps = readRDS('data/ppe_stepps.RDS')
# # ppe_stepps_ref = t(apply(ppe_stepps, 1, function(x) x/x[10]))
# # ppe_stepps_stat  = data.frame(taxa, t(apply(ppe_stepps_ref, 2, function(x) c(mean(x), sd(x)))))

###############################################################################################################################
## read in SVs
svs = read.csv('data/svs_LC6K.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
svs = svs[which(svs$use=='y'),]
svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)

svs_agg = svs_agg[which(svs_agg$taxon %in% taxa), ]

## construct param input data frame
reveals_inputs = data.frame(species=taxa,
                            fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa)],
                            PPEs=ppes$ppe,
                            PPE.errors=ppes$error)
rownames(reveals_inputs) = NULL


# # construct data frame of reveals inputs
# reveals_inputs = data.frame(species=taxa, 
#                             fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa)], 
#                             PPEs=ppe_stepps_stat[,2], 
#                             PPE.errors=ppe_stepps_stat[,3])
# rownames(reveals_inputs) = NULL
write.csv(reveals_inputs, "data/reveals_input/params.csv", row.names=FALSE)

##################################################################################################################################################
## run reveals
##################################################################################################################################################
# # start with only one slice
# slice = 6
# # pol_dat = pollen_bin[which(pollen_bin$slice_bin == 6),]
# ids     = unique(pol_dat$dataset)

# pollen_bin = pollen_bin[which(!is.na(pollen_bin$slice_bin)),]

# # load package
library(DISQOVER)
# load('DISQOVER/R/sysdata.rda')
source('DISQOVER/R/main.R')

ids = unique(pollen_bin$dataset)
pol_dat = pollen_bin

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


for (i in 1:length(ids)){
  
  print(i)
  id = ids[i]
  
  counts_site = data.frame(ages=slice_labels[pol_dat[which(pol_dat$dataset == id),'slice_bin']], 
                           pol_dat[which(pol_dat$dataset == id),which(colnames(pol_dat) %in% tolower(taxa))])
  colnames(counts_site) = c('ages', taxa)
  # counts_site = counts_site[which(rowSums(counts_site[,2:ncol(counts_site)])!=0),]
  
  basin_radius = pol_dat[which(pol_dat$dataset == id), c('lake_size')][1]
  if (is.na(basin_radius)){
    basin_radius = random_lake_size(pol_dat$lake_size)
  }
  if (basin_radius==0){
    basin_radius = random_lake_size(pol_dat$lake_size)
  }
  
  if (basin_radius<5){
    print(paste0("Skipping record ", id, "; basin radius <5 m."))
    next
  }
  
  coords_site = data.frame(pol_dat[which(pol_dat$dataset == id), c('dataset', 'long', 'lat', 'cell_id')][1,], basin_radius)
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
  if ((2*round(basin_radius))<40000){
    regionCutoff = 50000
  } else if ((40000<=(2*round(basin_radius))) & ((2*round(basin_radius))<90000)){
    regionCutoff = 100000
  } else if ((90000<=(2*round(basin_radius))) & ((2*round(basin_radius))<150000)){
    regionCutoff = 150000  
  } else {
    print("Lake size is uncharacteristically large.")
  }
    
  a <- REVEALSinR(pollenFile = "data/reveals_input/reveals_input.csv",
                  pf         = "data/reveals_input/params.csv",
                  filetype   = "csv",
                  dwm        = "gpm neutral",
                  tBasin     = "lake",
                  dBasin     = 2*round(basin_radius), # diameter!
                  regionCutoff = regionCutoff,
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
  
}

saveRDS(veg_pred, 'data/cache/veg_pred.RDS')

veg_pred <- readRDS('data/cache/veg_pred.RDS')


##################################################################################################################################################
## process output
##################################################################################################################################################

pdf('figures/test_map.pdf')
plot(pollen_bin$long, pollen_bin$lat)
map('world', add=TRUE)
dev.off()

# how to get the coords for grid cells
coords   = xyFromCell(grid, veg_pred$cell_id)
veg_pred = cbind(coords, veg_pred)


# # na_shp@data$id <- rownames(na_shp@data)
# # na_shp_ll <- spTransform(na_shp, CRS("+proj=longlat +datum=WGS84"))
# na_fort <- spTransform(na_shp, CRS("+init=epsg:4326"))
# # na_fort <- fortify(na_shp_ll, region='id')

# na_shp <- readOGR('data/map_data/na/NA_States_Provinces_Albers.shp')
# # na_shp@data$id <- rownames(na_shp@data)
# # na_shp_ll <- spTransform(na_shp, CRS("+proj=longlat +datum=WGS84"))
# na_fort <- spTransform(na_shp, CRS("+init=epsg:4326"))
# # na_fort <- fortify(na_shp_ll, region='id')

veg_grid = aggregate(mediansim ~ taxon + ages + cell_id + x+ y, veg_pred, sum)
veg_cast = dcast(veg_grid, cell_id + x + y + ages ~ taxon, value.var='mediansim')
veg_cast[,5:ncol(veg_cast)] = t(apply(veg_cast[,5:ncol(veg_cast)], 1, function(x) x/sum(x)))

veg_grid = melt(veg_cast, id.vars=c('cell_id', 'x', 'y', 'ages'))

saveRDS(veg_grid, 'data/cache/veg_grid.RDS')

veg_grid <- readRDS('data/cache/veg_grid.RDS')



# c(-100, 40, -60, 60)

# # make time bins
# # start with: 0.1-0.35k BP, 0.35-0.7k BP, 5.5-6.2k BP
# breaks = c(0.1, 0.35, 0.7, 2.7, 3.2, 5.7, 6.2)
# slice_bins = seq(1, 6)
# slice_labels = c(50, 200, 500, 1500, 3000, 6000)

slice_labels_plot = c(200, 500, 6000)
nslices_plot = length(slice_labels_plot)
veg_grid = veg_grid[which(veg_grid$ages %in% slice_labels_plot),]

# cuts for plotting
plot_cuts = c(0, 1e-5, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
n_plot_cuts = length(plot_cuts) - 1
veg_grid$value_bin = cut(veg_grid$value, plot_cuts, labels=FALSE)
labels = c("0", "0 - 2", "2 - 4", "4 - 6", "6 - 8", "8 - 10", "10 - 20", "20 - 40", "40 - 60", "60 - 80", "80 - 100")

veg_grid$value_bin[which(is.na(veg_grid$value_bin))] = 1

ntaxa = length(taxa)

pdf('figures/reveals_taxon_slices_by_taxon_binned.pdf')
for (i in 1:ntaxa){
  
  taxon   = taxa[i]
  veg_sub = veg_grid[which(veg_grid$variable==taxon),]
  
  p <- ggplot(data=veg_sub)
  p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
  p <- p + geom_tile(aes(x=x, y=y, fill=factor(value_bin)), data=veg_sub)
  p <- p + scale_fill_manual(values=tim.colors(n_plot_cuts), labels=labels, name="Percent")
  # p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + facet_grid(ages~.)
  p <- p + coord_equal(xlim = longlimits, ylim = latlimits, ratio=1)
  p <- p + ggtitle(paste0(taxon))
  print(p)
}
dev.off()


pdf('figures/reveals_taxon_slices_by_slice_binned.pdf')
for (i in 1:nslices_plot){
  
  age = slice_labels_plot[i]
  veg_sub = veg_grid[which(veg_grid$ages==age),]
  
  p <- ggplot(data=veg_sub)
  p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
  p <- p + geom_tile(aes(x=x, y=y, fill=factor(value_bin)), data=veg_sub)
  p <- p + scale_fill_manual(values=tim.colors(n_plot_cuts), labels=labels, name="Percent")
  # p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + facet_wrap(~variable, ncol=5)
  # p <- p + facet_grid(ages~.)
  # p <- p + coord_cartesian(xlim = longlimits, ylim = latlimits)
  p <- p + coord_equal(xlim = longlimits, ylim = latlimits, ratio=1)
  p <- p + ggtitle(paste0(age, ' YBP slice'))
  print(p)
}
dev.off()

# now aggregate to LCT
taxon2pft = read.csv('data/taxon2LCT_translation.csv')
veg_grid$LCT = taxon2pft$LCT[match(veg_grid$variable, taxon2pft$taxon)]

veg_lct = aggregate(value~cell_id+x+y+ages+LCT, veg_grid, sum, na.rm=TRUE)
veg_lct$value_bin = cut(veg_lct$value, plot_cuts, labels=FALSE)
veg_lct$value_bin[which(is.na(veg_lct$value_bin))] = 1

pdf('figures/reveals_pft_binned.pdf', width=12, height=6)
p <- ggplot(data=veg_lct)
p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
p <- p + geom_tile(aes(x=x, y=y, fill=factor(value_bin)), data=veg_lct)
p <- p + scale_fill_manual(values=tim.colors(n_plot_cuts), labels=labels, name="Percent")
# p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + facet_grid(ages~LCT)
p <- p + coord_equal(xlim = longlimits, ylim = latlimits, ratio=1)
print(p)
# ggsave('figures/reveals_pft_binned.pdf')
dev.off()

saveRDS(veg_lct, "data/cache/veg_lct.RDS")

######################################################################################################################################
## plot maps of pollen proportions
######################################################################################################################################

pol_melt = melt(pol_dat, id=c('cell_id', 'lat', 'long', 'slice_bin'), measure.vars=tolower(taxa))

pol_agg = aggregate(value~cell_id+lat+long+slice_bin+variable, pol_melt, sum, na.rm=TRUE)
pol_agg = dcast(pol_agg, cell_id + lat + long + slice_bin ~ variable, value.var='value')
pol_agg[,tolower(taxa)] = pol_agg[,tolower(taxa)]/rowSums(pol_agg[,tolower(taxa)])

pol_melt2 = melt(pol_agg, id=c('cell_id', 'lat', 'long', 'slice_bin'), measure.vars=tolower(taxa))
pol_melt2$value_bin = cut(pol_melt2$value, plot_cuts, labels=FALSE)
pol_melt2$value_bin[which(is.na(pol_melt2$value_bin))] = 1
pol_melt2$ages = slice_labels[pol_melt2$slice_bin]


nslices_plot = length(slice_labels_plot)
pol_melt2 = pol_melt2[which(pol_melt2$ages %in% slice_labels_plot),]

# how to get the coords for grid cells
coords   = xyFromCell(grid, pol_melt2$cell_id)
pol_melt2 = cbind(coords, pol_melt2)

# cuts for plotting
plot_cuts = c(0, 1e-5, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
n_plot_cuts = length(plot_cuts) - 1
pol_melt2$value_bin = cut(pol_melt2$value, plot_cuts, labels=FALSE)
pol_melt2$value_bin[which(is.na(pol_melt2$value_bin))] = 1


pdf('figures/pollen_taxon_slices_by_taxon_binned.pdf')
for (i in 1:ntaxa){
  
  taxon   = taxa[i]
  veg_sub = pol_melt2[which(pol_melt2$variable==tolower(taxon)),]
  # veg_sub = veg_grid[which(veg_grid$variable==taxon),]
  
  p <- ggplot(data=veg_sub)
  p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
  # p <- p + geom_tile(aes(x=long, y=lat, fill=factor(value_bin)), data=veg_sub)
  p <- p + geom_tile(aes(x=x, y=y, fill=factor(value_bin)), data=veg_sub)
  
  p <- p + scale_fill_manual(values=tim.colors(n_plot_cuts), labels=labels, name="Percent")
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + facet_grid(ages~.)
  p <- p + coord_equal(xlim = longlimits, ylim = latlimits, ratio=1)
  p <- p + ggtitle(paste0(taxon))
  print(p)
}
dev.off()


# pdf('figures/reveals_taxon_slices.pdf')
# for (i in 1:nplotslice){
#   
#   age = plot_slices[i]
#   veg_sub = veg_grid[which(veg_grid$ages==age),]
#   
#   p <- ggplot(data=veg_sub)
#   p <- p + geom_tile(aes(x=x, y=y, fill=value), data=veg_sub)
#   p <- p + scale_fill_gradientn(colours=tim.colors(10))
#   p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
#   # p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
#   p <- p + theme_bw()
#   p <- p + theme(axis.text = element_blank(),
#                  axis.title = element_blank(),
#                  axis.ticks = element_blank())
#   p <- p + facet_wrap(~variable, ncol=3)
#   p <- p + coord_cartesian(xlim = longlimits, ylim = latlimits)
#   p <- p + ggtitle(paste0(age, ' YBP slice'))
#   print(p)
#   # ggsave('figures/reveals_taxon_50.pdf')
# }
# dev.off()
# 
# # now by taxon
# ntaxa = length(taxa)
# 
# pdf('figures/reveals_taxon_slices_by_taxon.pdf')
# for (i in 1:ntaxa){
#   
#   taxon   = taxa[i]
#   veg_sub = veg_grid[which(veg_grid$variable==taxon),]
#   
#   p <- ggplot(data=veg_sub)
#   p <- p + geom_tile(aes(x=x, y=y, fill=value), data=veg_sub)
#   p <- p + scale_fill_gradientn(colours=tim.colors(10))
#   p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
#   # p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
#   p <- p + theme_bw()
#   p <- p + theme(axis.text = element_blank(),
#                  axis.title = element_blank(),
#                  axis.ticks = element_blank())
#   p <- p + facet_wrap(~ages, ncol=3)
#   p <- p + coord_cartesian(xlim = longlimits, ylim = latlimits)
#   p <- p + ggtitle(paste0(taxon, ' YBP slice'))
#   print(p)
#   # ggsave('figures/reveals_taxon_50.pdf')
# }
# dev.off()

# ###################################################################################################################################
# 
# 
# 
# library(scatterpie)
# library(fields)
# library(reshape2)
# library(tidyr)
# library(ggforce)
# 
# cols = unique(veg_pred$taxon)
# 
# # geom_scatterpie <- function(mapping=NULL, data, cols, ...) {
# #   if (is.null(mapping))
# #     mapping <- aes_(x=~x, y=~y)
# #   mapping <- modifyList(mapping, aes_(r0=0, fill=~type,
# #                                       amount=~value))
# #   
# #   if (!'r' %in% names(mapping)) {
# #     xvar <- as.character(mapping)["x"]
# #     size <- diff(range(data[, xvar]))/50
# #     mapping <- modifyList(mapping, aes_(r=size))
# #   }
# #   
# #   names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")
# #   
# #   df <- gather_(data, "type", "value", cols)
# #   ## df$type <- factor(df$type, levels=cols)
# #   geom_arc_bar(mapping, data=df, stat='pie', inherit.aes=FALSE, ...)
# # }
# 
# 
# veg_cast = dcast(veg_pred, id + x + y + ages ~ taxon, value.var='meansim')
# 
# ages_list = ages*100
# pdf(file='figures/reveals_UMW_v2.pdf')
# for (age in ages_list){
#   veg_sub = subset(veg_cast, ages == age)
#   veg_sub$x = veg_sub$x*1e6
#   veg_sub$y = veg_sub$y*1e6
#   
#   p <- ggplot()
#   p <- p + geom_scatterpie(aes(x=x, y=y),
#                            data=veg_sub, cols=cols, alpha=.8)
#   p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
#   p <- p + theme_bw()
#   p <- p + theme(axis.text = element_blank(),
#                  axis.title = element_blank(),
#                  axis.ticks = element_blank())
#   p <- p + ggtitle(paste0(age, " YBP"))
#   print(p)
#   # ggsave(file=paste0('figures/reveals_UMW_', age, 'ybp.pdf'))
# }
# dev.off()
# ####################################################################################################################################
# cols = unique(veg_pred$taxon)
# 
# veg_cast = dcast(veg_pred, id + x + y + ages ~ taxon, value.var='meansim')
# source('r/utils/pie_funs.R')
# 
# centers   = centers_veg*1e6
# colnames(centers) = c('x', 'y')
# 
# xlo = min(centers[,1])
# xhi = max(centers[,1])
# ylo = min(centers[,2])
# yhi = max(centers[,2])
# y_shift=40000
# x_shift=100000
# 
# 
# ages_list = ages*100
# pdf(file='figures/reveals_UMW_ppies.pdf', width=12, height=10)
# for (age in ages_list){
#   veg_sub = subset(veg_cast, ages == age)
#   veg_sub$x = veg_sub$x*1e6
#   veg_sub$y = veg_sub$y*1e6
#   
#   
#   veg_props  = t(apply(veg_sub[,5:ncol(veg_sub)], 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))
#   coords = veg_sub[,c('x', 'y')]
#   
#   par(mfrow=c(1,1))
#   pieMap(proportions = veg_props, 
#          centers  = coords,
#          restrict = FALSE,
#          inputRestricted = FALSE,
#          xlim   = c(xlo+x_shift, xhi-x_shift),
#          ylim   = c(ylo+y_shift, yhi-y_shift),
#          radius = 18000,
#          scale  = 1,
#          xlab   = 'x',
#          ylab   = 'y',
#          add_legend = TRUE, 
#          main_title='')
# }
# dev.off()
# 
# ages_list = ages*100
# pdf(file='figures/reveals_UMW_v2.pdf')
# for (age in ages_list){
#   veg_sub = subset(veg_cast, ages == age)
#   veg_sub$x = veg_sub$x*1e6
#   veg_sub$y = veg_sub$y*1e6
#   
#   p <- ggplot()
#   p <- p + geom_scatterpie(aes(x=x, y=y),
#                            data=veg_sub, cols=cols, alpha=.8)
#   p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
#   p <- p + theme_bw()
#   p <- p + theme(axis.text = element_blank(),
#                  axis.title = element_blank(),
#                  axis.ticks = element_blank())
#   p <- p + ggtitle(paste0(age, " YBP"))
#   print(p)
#   # ggsave(file=paste0('figures/reveals_UMW_', age, 'ybp.pdf'))
# }
# dev.off()
# 
# ################################################################################################################################################
# ## make the grid
# ################################################################################################################################################
# 
# make_grid(meta_pol)
# 
# # BROKEN AFTER HERE
# 
# # by grid cell
# veg_cast = aggregate(meansim ~ idx_veg + ages + taxon, veg_pred, function(x) mean(x))
# 
# # veg_cast = dcast(veg_pred, idx_veg + ages ~ taxon, value.var='meansim', fun.aggregate = function(x) sum(x))
# # veg_cast[, 3:ncol(veg_cast)] = t(apply(veg_cast[, 3:ncol(veg_cast)], 1, function(x) x/sum(x)))
# # 
# veg_grid = data.frame(x=centers_veg$x[veg_cast$idx_veg], y=centers_veg$y[veg_cast$idx_veg], veg_cast)
# # 
# # veg_gridm = melt(veg_grid)
# 
# 
# ages_list = ages*100
# age=ages_list[1]
# veg_sub = subset(veg_grid, ages == age)
# veg_sub$x = veg_sub$x*1e6
# veg_sub$y = veg_sub$y*1e6
# 
# p <- ggplot(data=veg_sub)
# p <- p + geom_tile(aes(x=x, y=y, fill=meansim), data=veg_sub) 
# p <- p + scale_fill_gradientn(colours=tim.colors())
# p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
# p <- p + theme_bw()
# p <- p + theme(axis.text = element_blank(),
#                axis.title = element_blank(),
#                axis.ticks = element_blank())
# p <- p + facet_wrap(~taxon, ncol=3)
# print(p)
# 
# 
# veg_sub = subset(veg_grid, (taxon %in% c('HEMLOCK', 'PINE', 'OAK', 'ELM')) & (ages %in% c(200, 500, 1000, 1500, 2000)))
# veg_sub = subset(veg_grid, (ages %in% c(200, 500, 1000, 1500, 2000)))
# 
# veg_sub$x = veg_sub$x*1e6
# veg_sub$y = veg_sub$y*1e6
# p <- ggplot(data=veg_sub)
# p <- p + geom_tile(aes(x=x, y=y, fill=meansim), data=veg_sub) 
# p <- p + geom_point(aes(x=x, y=y, colour=meansim), data=veg_sub, size=1.2) 
# p <- p + scale_fill_gradientn(colours=tim.colors())
# p <- p + scale_colour_gradientn(colours=tim.colors())
# p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
# p <- p + theme_bw()
# p <- p + theme(axis.text = element_blank(),
#                axis.title = element_blank(),
#                axis.ticks = element_blank())
# p <- p + facet_grid(ages~taxon)
# print(p)
# ggsave('figures/reveals_UMW_grid.pdf')
# 
# # # plot pie maps
# # # postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
# # pdf(paste0('figures/pie_plot_pls_ALL_v', v, '.pdf'), width=12, height=10)
# # par(mfrow=c(1,1))
# # pieMap(proportions = pls_props, 
# #        centers  = knots_in,
# #        restrict = FALSE,
# #        inputRestricted = FALSE,
# #        xlim   = c(xlo+shift, xhi-shift),
# #        ylim   = c(ylo+shift, yhi-shift),
# #        radius = 14000,
# #        scale  = 1,
# #        xlab   = 'x',
# #        ylab   = 'y',
# #        add_legend = FALSE, 
# #        main_title='')
# # dev.off()
# 
# 
# ##################################################################################################################################################
# 
# 
# 
# 
# dat = read.csv('reveals/REVEALS.Lake.BigJohnPond.csv', row.names=NULL, header=FALSE)
# dat = dat[-c(2,3),]
# dat_flip = t(dat[, -1])
# colnames(dat_flip) = as.vector(dat[,1])
# # write.table(dat_flip, 'revealsR/reveals_BigJohnPond.csv', col.names=FALSE, sep=',')
# write.csv(dat_flip, 'revealsR/reveals_BigJohnPond.csv', row.names=FALSE)
# 
# 
# # with english csv files
# a <- REVEALSinR(pollenFile = "revealsR/reveals_BigJohnPond.csv", 
#                 pf="revealsR/reveals_inputs.csv", 
#                 dwm="LSM unstable", 
#                 tBasin="lake", 
#                 dBasin=600, 
#                 regionCutoff=100000, 
#                 repeats=1000)
# 
# foo = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
# foo$variable = gsub('OTHER.', 'OTHER_', foo$variable)
# 
# foo$taxon = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
# foo$type = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))
# 
# ggplot(data=subset(foo, type %in% 'meansim')) + geom_line(aes(x=ages, y=value, colour=taxon)) + theme_bw()
# 
# 
# # with english csv files
# gpm <- REVEALSinR(pollenFile = "revealsR/reveals_BigJohnPond.csv", 
#                 pf="revealsR/reveals_inputs.csv", 
#                 dwm="GPM neutral", 
#                 tBasin="lake", 
#                 dBasin=600, 
#                 regionCutoff=100000, 
#                 repeats=1000)
# 
# foo = melt(gpm, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
# foo$variable = gsub('OTHER.', 'OTHER_', foo$variable)
# 
# foo$taxon = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
# foo$type = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))
# 
# ggplot(data=subset(foo, type %in% 'meansim')) + geom_line(aes(x=ages, y=value, colour=taxon)) + theme_bw()
# 
# 
# # after REVEALS on all sites do pie plot as in STEPPS calibration paper
# 
# 
# 
# 
# 
# 
# 
# # # with english csv files
# # a <- REVEALSinR(pollenFile = "REVEALSinR Lake Viitina_en.csv", 
# #                 pf="REVEALSinR PPEs_en.csv", 
# #                 dwm="LSM unstable", 
# #                 tBasin="lake", 
# #                 dBasin=600, 
# #                 regionCutoff=100000, 
# #                 repeats=1000)

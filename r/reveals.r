library(DISQOVER)
library(reshape2)
library(ggplot2)

# # load package
# load('DISQOVER/R/sysdata.rda')

source('DISQOVER/R/main.R')
source('generate_tables.R')

##################################################################################################################################################
## pulle pollen data for NA
##################################################################################################################################################
library(dplyr)
library(neotoma)
library(stepps)

if(!'na_downloads.rds' %in% list.files('data/cache/')) {
  canada <- get_dataset(gpid='Canada', datasettype = 'pollen') %>% get_download
  usa <-  get_dataset(gpid='United States', datasettype = 'pollen') %>% get_download

  na_pollen <- bind(canada, usa)

  saveRDS(na_pollen, file = 'data/cache/na_downloads.rds')
} else {
  na_pollen <- readRDS('data/cache/ena_downloads.rds')
}

if(!'ena_downloads.rds' %in% list.files('data/cache/')) {
  
  ena <- get_dataset(loc=c(-100, 40, -60, 60), datasettype = 'pollen') %>% get_download
  
  saveRDS(ena, file = '../data/cache/ena_downloads.rds')
  
} else {
  ena_pollen <- readRDS('../data/cache/ena_downloads.rds')
}

ena_taxa <- lapply(taxa(ena_pollen, collapse = FALSE), as.data.frame) %>% bind_rows

generate_tables(ena_pollen, output = 'things.csv')

















# read in pollen site meta data
# meta   = read.csv('data/site_meta.csv') # do we need this? don't think so
# counts = read.csv('data/sediment_ages_v1.0_varves.csv')

load('data/input.rdata')

# read in PPEs
ppes = read.csv('data/ppes_v2.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
# ppes$taxon[which(ppes$taxon == "OTHER.CONIFER")]  = "FIR"
# ppes$taxon[which(ppes$taxon == "OTHER.HARDWOOD")] = "OTHER HARDWOOD"

# think carefully about how to rescale these
maple_rescale = mean(ppes[which((ppes$taxon=='MAPLE')&(ppes$continent == 'North America')&(ppes$tag != 'calcote')),'ppe'])
ppes[which(ppes$tag == 'calcote'),'ppe'] = ppes[which(ppes$tag == 'calcote'),'ppe']*0.61

ggplot(data=subset(ppes, continent == 'North America')) + geom_point(aes(x=ppe, y=taxon))

ppes_agg = aggregate(ppe ~ taxon + continent, ppes, function(x) quantile(x, c(0.05, 0.5, 0.95)))
ppes_sd = aggregate(ppe ~ taxon + continent, ppes, function(x) sd(x))
ppes_sd[which(is.na(ppes_sd$ppe)), 'ppe'] = 0.01

ppes_agg_NA = subset(ppes_agg, continent == 'North America')

ggplot(data=ppes_agg_NA) + geom_point(aes(x=ppe[,2], y=reorder(taxon, ppe[,2])))

# read in STEPPS PPEs
ppe_stepps = readRDS('data/ppe_post.RDS')
# ppe_stepps = readRDS('data/ppe_stepps.RDS')
ppe_stepps_ref = t(apply(ppe_stepps, 1, function(x) x/x[10]))
ppe_stepps_stat  = data.frame(taxa, t(apply(ppe_stepps_ref, 2, function(x) c(mean(x), sd(x)))))

# read in SVs
svs = read.csv('data/svs_meta3.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
svs_agg = aggregate(sv ~ taxon, svs, median)

# construct data frame of reveals inputs
reveals_inputs = data.frame(species=taxa,
                            fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa)],
                            PPEs=ppes_agg_NA$ppe[,2],
                            PPE.errors=ppes_agg_NA$ppe[,2][match(ppes_agg_NA$taxon, taxa)])
rownames(reveals_inputs) = NULL

# # construct data frame of reveals inputs
# reveals_inputs = data.frame(species=taxa, 
#                             fallspeed=svs_agg$sv[match(svs_agg$taxon, taxa)], 
#                             PPEs=ppe_stepps_stat[,2], 
#                             PPE.errors=ppe_stepps_stat[,3])
# rownames(reveals_inputs) = NULL
write.csv(reveals_inputs, "revealsR/reveals_inputs.csv", row.names=FALSE)

# read in lake size, and only keep large lakes
# lake_size = read.csv('data/midwest_lake_areas.csv')
# lake_size2 = read.csv('data/neotoma_and_GNIS.csv') 
lake_size = read.csv('data/areas_with_datasetID.csv')
lake_size = lake_size[!duplicated(lake_size$DatasetID),]

lake_dists = lake_size$AREAHA[which(!is.na(lake_size$AREAHA))]

veg_pred = data.frame(id=numeric(0), 
                      x=numeric(0), 
                      y=numeric(0), 
                      idx_veg = numeric(0),
                      radius=numeric(0),
                      ages=numeric(0), 
                      taxon=character(0), 
                      meansim=numeric(0), 
                      mediansim=numeric(0), 
                      q10sim=numeric(0), 
                      q90sim=numeric(0), 
                      sdsim=numeric(0))

# 
# id = 29
ids = unique(meta_pol$id)

for (i in 1:length(ids)){
  
  print(i)
  id = ids[i]
  
  counts_site = data.frame(meta_pol$age[which(meta_pol$id == id)]*100, y[which(meta_pol$id == id), ])
  colnames(counts_site) = c('ages', taxa)
  counts_site = counts_site[which(rowSums(counts_site[,2:ncol(counts_site)])!=0),]
  
  coords_site = data.frame(meta_pol[which(meta_pol$id == id), c('id', 'x', 'y')][1,])
  rownames(coords_site) = NULL
  
  # prob don't need this step - fix later
  write.csv(counts_site, 'revealsR/reveals_test.csv', row.names=FALSE)
  
  if (!any(lake_size$DatasetID == id)){
    print("Lake not in lake area file. Skipping reconstruction.")
    next
  } else if (is.na(lake_size$AREAHA[lake_size$DatasetID == id]) & is.na(lake_size$Area[lake_size$DatasetID == id])) {
    print("Lake area missing. Skipping reconstruction.")
    next
  }
  
  # in HA; convert to radius in m
  # pi * r * r
  if (!is.na(lake_size$AREAHA[lake_size$DatasetID == id])){
    lake_rad_site = sqrt(lake_size$AREAHA[which(lake_size$DatasetID == id)]*0.01 / pi)*1000
  } else {
    lake_rad_site = sqrt(lake_size$Area[which(lake_size$DatasetID == id)]*0.01 / pi)*1000
  }
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
  a <- REVEALSinR(pollenFile = "revealsR/reveals_test.csv",
                  pf         = "revealsR/reveals_inputs.csv",
                  dwm        = "lsm unstable",
                  tBasin     = "lake",
                  dBasin     = 2*round(lake_rad_site), # diameter!
                  regionCutoff = 100000,
                  repeats      = 1000)
  
  veg = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
  veg$variable = gsub('OTHER.', 'OTHER_', veg$variable)
  
  veg$taxon = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
  veg$type = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))
  
  veg_cast = dcast(veg, ages + taxon ~ type)
  
  ggplot(data=veg_cast) + 
    geom_line(aes(x=ages, y=mediansim, colour=taxon)) + 
    geom_point(aes(x=ages, y=mediansim, colour=taxon)) + 
    geom_errorbar(aes(x=ages, ymin=q10sim, ymax=q90sim, colour=taxon), width=1.5) +
    theme_bw()
  
  veg_pred = rbind(veg_pred, data.frame(coords_site, idx_veg=idx_cores[i], radius=lake_rad_site, veg_cast))
  
}


veg_agg = aggregate(meansim~ages+taxon, veg_pred, mean)
ggplot(data=veg_agg) + geom_point(aes(x=ages, y=meansim, colour=taxon))

library(maptools)

us.shp <- readShapeLines('data/map_data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))
us.shp@data$id <- rownames(us.shp@data)
us.fort <- fortify(us.shp, region='id') 


add_map_albers <- function(plot_obj, map_data=us.fort, limits){
  p <- plot_obj + geom_path(data=map_data, aes(x=long, y=lat, group=group),  colour='grey55') + 
    #     scale_x_continuous(limits = c(min(umw.coord$x, na.rm=TRUE), max(umw.coord$x, na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(umw.coord$y, na.rm=TRUE), max(umw.coord$y, na.rm=TRUE)))#, colour = "black", size = 1, fill = "white", aes(x=long, y=lat, group = group))
    # #   
    #     scale_x_continuous(limits = c(min(dat[,1], na.rm=TRUE), max(dat[,1], na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(dat[,2], na.rm=TRUE), max(dat[,2], na.rm=TRUE)))
    scale_x_continuous(limits = limits$xlims*1000000) +
    scale_y_continuous(limits = limits$ylims*1000000) #+ coord_map("albers")
  return(p)
  
}


get_limits <- function(centers){
  xlo = min(centers[,1])
  xhi = max(centers[,1])
  
  ylo = min(centers[,2])
  yhi = max(centers[,2])
  
  return(list(xlims=c(xlo,xhi),ylims=c(ylo, yhi)))
}  

limits = get_limits(centers_pls)

library(scatterpie)
library(fields)
library(reshape2)
library(tidyr)
library(ggforce)

cols = unique(veg_pred$taxon)

# geom_scatterpie <- function(mapping=NULL, data, cols, ...) {
#   if (is.null(mapping))
#     mapping <- aes_(x=~x, y=~y)
#   mapping <- modifyList(mapping, aes_(r0=0, fill=~type,
#                                       amount=~value))
#   
#   if (!'r' %in% names(mapping)) {
#     xvar <- as.character(mapping)["x"]
#     size <- diff(range(data[, xvar]))/50
#     mapping <- modifyList(mapping, aes_(r=size))
#   }
#   
#   names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")
#   
#   df <- gather_(data, "type", "value", cols)
#   ## df$type <- factor(df$type, levels=cols)
#   geom_arc_bar(mapping, data=df, stat='pie', inherit.aes=FALSE, ...)
# }


veg_cast = dcast(veg_pred, id + x + y + ages ~ taxon, value.var='meansim')

ages_list = ages*100
pdf(file='figures/reveals_UMW_v2.pdf')
for (age in ages_list){
  veg_sub = subset(veg_cast, ages == age)
  veg_sub$x = veg_sub$x*1e6
  veg_sub$y = veg_sub$y*1e6
  
  p <- ggplot()
  p <- p + geom_scatterpie(aes(x=x, y=y),
                           data=veg_sub, cols=cols, alpha=.8)
  p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + ggtitle(paste0(age, " YBP"))
  print(p)
  # ggsave(file=paste0('figures/reveals_UMW_', age, 'ybp.pdf'))
}
dev.off()
####################################################################################################################################
cols = unique(veg_pred$taxon)

veg_cast = dcast(veg_pred, id + x + y + ages ~ taxon, value.var='meansim')
source('r/utils/pie_funs.R')

centers   = centers_veg*1e6
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])
y_shift=40000
x_shift=100000


ages_list = ages*100
pdf(file='figures/reveals_UMW_ppies.pdf', width=12, height=10)
for (age in ages_list){
  veg_sub = subset(veg_cast, ages == age)
  veg_sub$x = veg_sub$x*1e6
  veg_sub$y = veg_sub$y*1e6
  
  
  veg_props  = t(apply(veg_sub[,5:ncol(veg_sub)], 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))
  coords = veg_sub[,c('x', 'y')]
  
  par(mfrow=c(1,1))
  pieMap(proportions = veg_props, 
         centers  = coords,
         restrict = FALSE,
         inputRestricted = FALSE,
         xlim   = c(xlo+x_shift, xhi-x_shift),
         ylim   = c(ylo+y_shift, yhi-y_shift),
         radius = 18000,
         scale  = 1,
         xlab   = 'x',
         ylab   = 'y',
         add_legend = TRUE, 
         main_title='')
}
dev.off()

ages_list = ages*100
pdf(file='figures/reveals_UMW_v2.pdf')
for (age in ages_list){
  veg_sub = subset(veg_cast, ages == age)
  veg_sub$x = veg_sub$x*1e6
  veg_sub$y = veg_sub$y*1e6
  
  p <- ggplot()
  p <- p + geom_scatterpie(aes(x=x, y=y),
                           data=veg_sub, cols=cols, alpha=.8)
  p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + ggtitle(paste0(age, " YBP"))
  print(p)
  # ggsave(file=paste0('figures/reveals_UMW_', age, 'ybp.pdf'))
}
dev.off()

################################################################################################################################################
## make the grid
################################################################################################################################################

make_grid(meta_pol)

# BROKEN AFTER HERE

# by grid cell
veg_cast = aggregate(meansim ~ idx_veg + ages + taxon, veg_pred, function(x) mean(x))

# veg_cast = dcast(veg_pred, idx_veg + ages ~ taxon, value.var='meansim', fun.aggregate = function(x) sum(x))
# veg_cast[, 3:ncol(veg_cast)] = t(apply(veg_cast[, 3:ncol(veg_cast)], 1, function(x) x/sum(x)))
# 
veg_grid = data.frame(x=centers_veg$x[veg_cast$idx_veg], y=centers_veg$y[veg_cast$idx_veg], veg_cast)
# 
# veg_gridm = melt(veg_grid)


ages_list = ages*100
age=ages_list[1]
veg_sub = subset(veg_grid, ages == age)
veg_sub$x = veg_sub$x*1e6
veg_sub$y = veg_sub$y*1e6

p <- ggplot(data=veg_sub)
p <- p + geom_tile(aes(x=x, y=y, fill=meansim), data=veg_sub) 
p <- p + scale_fill_gradientn(colours=tim.colors())
p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + facet_wrap(~taxon, ncol=3)
print(p)


veg_sub = subset(veg_grid, (taxon %in% c('HEMLOCK', 'PINE', 'OAK', 'ELM')) & (ages %in% c(200, 500, 1000, 1500, 2000)))
veg_sub = subset(veg_grid, (ages %in% c(200, 500, 1000, 1500, 2000)))

veg_sub$x = veg_sub$x*1e6
veg_sub$y = veg_sub$y*1e6
p <- ggplot(data=veg_sub)
p <- p + geom_tile(aes(x=x, y=y, fill=meansim), data=veg_sub) 
p <- p + geom_point(aes(x=x, y=y, colour=meansim), data=veg_sub, size=1.2) 
p <- p + scale_fill_gradientn(colours=tim.colors())
p <- p + scale_colour_gradientn(colours=tim.colors())
p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + facet_grid(ages~taxon)
print(p)
ggsave('figures/reveals_UMW_grid.pdf')

# # plot pie maps
# # postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
# pdf(paste0('figures/pie_plot_pls_ALL_v', v, '.pdf'), width=12, height=10)
# par(mfrow=c(1,1))
# pieMap(proportions = pls_props, 
#        centers  = knots_in,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 14000,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y',
#        add_legend = FALSE, 
#        main_title='')
# dev.off()


##################################################################################################################################################




dat = read.csv('reveals/REVEALS.Lake.BigJohnPond.csv', row.names=NULL, header=FALSE)
dat = dat[-c(2,3),]
dat_flip = t(dat[, -1])
colnames(dat_flip) = as.vector(dat[,1])
# write.table(dat_flip, 'revealsR/reveals_BigJohnPond.csv', col.names=FALSE, sep=',')
write.csv(dat_flip, 'revealsR/reveals_BigJohnPond.csv', row.names=FALSE)


# with english csv files
a <- REVEALSinR(pollenFile = "revealsR/reveals_BigJohnPond.csv", 
                pf="revealsR/reveals_inputs.csv", 
                dwm="LSM unstable", 
                tBasin="lake", 
                dBasin=600, 
                regionCutoff=100000, 
                repeats=1000)

foo = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
foo$variable = gsub('OTHER.', 'OTHER_', foo$variable)

foo$taxon = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
foo$type = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))

ggplot(data=subset(foo, type %in% 'meansim')) + geom_line(aes(x=ages, y=value, colour=taxon)) + theme_bw()


# with english csv files
gpm <- REVEALSinR(pollenFile = "revealsR/reveals_BigJohnPond.csv", 
                pf="revealsR/reveals_inputs.csv", 
                dwm="GPM neutral", 
                tBasin="lake", 
                dBasin=600, 
                regionCutoff=100000, 
                repeats=1000)

foo = melt(gpm, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
foo$variable = gsub('OTHER.', 'OTHER_', foo$variable)

foo$taxon = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
foo$type = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))

ggplot(data=subset(foo, type %in% 'meansim')) + geom_line(aes(x=ages, y=value, colour=taxon)) + theme_bw()


# after REVEALS on all sites do pie plot as in STEPPS calibration paper







# # with english csv files
# a <- REVEALSinR(pollenFile = "REVEALSinR Lake Viitina_en.csv", 
#                 pf="REVEALSinR PPEs_en.csv", 
#                 dwm="LSM unstable", 
#                 tBasin="lake", 
#                 dBasin=600, 
#                 regionCutoff=100000, 
#                 repeats=1000)

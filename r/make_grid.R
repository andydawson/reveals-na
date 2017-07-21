make_grid <- function(pollen, coord_fun = ~ x + y, projection = '+init=epsg:3175', resolution = 1) {
  sp::coordinates(pollen) <- coord_fun
  sp::proj4string(pollen) <- sp::CRS(projection)
  latlong <- sp::spTransform(pollen, '+init=epsg:4326')

  borders <- sp::bbox(latlong)

  out_rast <- raster::raster(xmn = floor(borders[1,1]),
                             xmx = ceiling(borders[1,2]),
                             ymn = floor(borders[2,1]),
                             ymx = ceiling(borders[2,2]),
                             resolution = 1,
                             crs = CRS('+init=epsg:4326'))

  out_rast <- raster::setValues(out_rast, 1:ncell(out_rast))
  return(out_rast)
}

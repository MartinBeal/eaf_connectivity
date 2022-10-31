## Group neighbouring hexcells into 'sites' 
# ***** CURRENTLY NOT WORKING *****

x <- st_intersects(grid_re)

# get all pairwise intersections between hexcells
ovrlp_grps <- lapply(x, function(y){
  # y <- x[[1]]
  ovrlps_list <- list()
  for(i in seq_along(y)){
    ovrlps <- x[[y[i]]]
    ovrlps_list[[i]] <- ovrlps
  }
  y_ovrlps <- paste(unique(do.call(c,ovrlps_list)), collapse=" ")
  return(y_ovrlps)
})

xx <- unique(do.call(c, ovrlp_grps))

xxdf <- data.frame(id = paste0("site_", seq_along(xx)), xx)

## identify unique groupings of neighboring cells ----------------------------
grid_re$siteID <- NA
for(i in seq_len(nrow(grid_re))){
  for(w in seq_along(xx)){
    onegrp <- as.numeric(str_split(xx[w], pattern = " ")[[1]])
    if(i %in% onegrp){
      grid_re$siteID[i] <- paste0("site_", w)
    } else {next}
  }
  
}

## make site IDs a numeric order (no skips)
grid_re$siteID <- paste0("site_", 
                         as.numeric(as.factor((grid_re$siteID))))

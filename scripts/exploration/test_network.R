## create simple UNDIRECTED network from test BT godwit cr-data

pacman::p_load(igraph, stringr, tictoc, tidygraph)


## choose data subset to run

# netdat <- subset(alldat, combo %in% unique(alldat$combo)[1:100])
netdat <- alldat ## all individuals

## Fix some bad country data

netdat$country <- ifelse(netdat$country == "Morroco", 
                         "Morocco", netdat$country)
netdat$country <- ifelse(
  netdat$country %in% c("Scotland", "England", 
                        "Channel Island", "Lincolnshire", "Northern Ireland", 
                        "Wales"), "United Kingdom", netdat$country)
netdat$country <- ifelse(netdat$country == "Faroe Island", 
                         "Denmark", netdat$country)

## choose which location column to use as 'site'
# netdat$site_num <- as.numeric(as.factor(netdat$site))
# site_summ <- netdat %>% group_by(site_num, site) %>% summarise()

netdat$site_num <- as.numeric(as.factor(netdat$country))
site_summ <- netdat %>% group_by(site_num, country) %>% summarise()




## remove sightings w/ no site information (i.e. NAs in site info)
netdat %<>% subset(!is.na(site_num))

## Edge list ------------------------------------------------------------------

## produce all combinations of sites visited by each individual
netdat_list <- split(netdat, netdat$combo)

tic()
netdat_list <- lapply(
  seq_along(netdat_list), 
  function(x){
    # print(x)
    one <- netdat_list[[x]]
    xx <- as.data.frame(
      RcppAlgos::comboGrid(
        one$site_num,
        one$site_num, repetition = TRUE)
    )
    xx$combo <- one$combo[1]
    return(xx)
  })
toc()

full <- data.table::rbindlist(netdat_list)

## remove self connections
noself <- full[-which(full$Var1 == full$Var2), ]

## combine sites into single variable for summarizing
noself$sitecomb <- paste(noself$Var1, noself$Var2)

## retain only unique site combos, summ how many individuals for connex
edgelist <- noself %>% group_by(sitecomb) %>% 
  summarise(
    n_id = n()
  )

## re-split site combos into separate columns
edgelist$site1 <- as.integer(
  do.call(
    rbind, str_split(edgelist$sitecomb, pattern = " "))[,1])
edgelist$site2 <- as.integer(
  do.call(
    rbind, str_split(edgelist$sitecomb, pattern = " "))[,2])

edgelist %<>% dplyr::select(site1, site2, n_id)


## Vertex list ---------------------------------------------------------------
nodelist <- netdat %>% group_by(site_num) %>% 
  summarise(
    n_id = n_distinct(combo),
    n_obs  = n()
    )

## Remove sites visited by <= X bird(s)
edgelist %<>% subset(n_id > 2)
nodelist %<>% subset(n_id > 2)

## igraph network object
routes_igraph <- graph_from_data_frame(
  d = edgelist, 
  vertices = nodelist,
  directed = F)

plot(routes_igraph, vertex.size = 5, vertex.label=NA)

## convert to tidygraph tbl_graph object

tgraph <- as_tbl_graph(routes_igraph)


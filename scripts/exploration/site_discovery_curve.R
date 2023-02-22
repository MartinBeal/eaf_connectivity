## Assessing representativeness of individual relocation data for 
# identifying sites

pacman::p_load(dplyr, foreach, stringr)

## choose data subset to run --------------------------------------------------

## Load combined data, with both IBA overlay info and hexgrid outsites -------
### SAVE 
alldat <- readRDS(
  paste0("data/analysis/combined/alldatatypes_ibas_outsites_all.rds"))

alldat <- mutate(alldat, 
  datatype = ifelse(datatype == "trax", "tracking", datatype),
  datatype = factor(datatype, levels=c("metal", "color", "tracking"))
)

iteration <- 10

dtypes <- c("metal", "color", "tracking")

listedResults <- list()

## loop through each datatype
for(i in seq_along(dtypes)){

  dtype <- dtypes[i]
  
  ## filter to one datatype
  sitedat <- subset(alldat, str_detect(alldat$datatype, dtype))
  
  sitedat %>% group_by(bird_id) %>% 
    summarise(n_site = n_distinct(SitRecID)) %>% 
    ungroup() %>% 
    summarise(
      mn_nsite = mean(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  ## reduce number of individuals to speed it up
  # sitedat <- filter(sitedat, bird_id %in% unique(alldat$bird_id)[1:20])
  
  
  ###
  ids  <- unique(sitedat$bird_id) # get IDs, after any filtered out
  nids <- length(ids)
  
  Nloop <- seq(1, (nids - 3), 1)
  if (nids >= 100) {
    Nloop <- c(seq(1, 25, 1), seq(26, 49, 3), seq(50, (nids - 1), 6))
  }
  if (nids >= 200) {
    Nloop <- seq(1, (nids - 3), 10)
  }
  DoubleLoop <- data.frame(
    sample_size = rep(Nloop, each = iteration),
    iteration = rep(seq_len(iteration), length(Nloop))
  )
  LoopNr <- seq_len(dim(DoubleLoop)[1])
  
  
  ### Parallel or sequential processing?----------------------------------------
  nCores <- 10 # of 12
  
  if (nCores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE) |
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages \"parallel\" and \"doParallel\" needed for this function
        to work. Please install.",
           call. = FALSE)
    }
    maxCores <- parallel::detectCores()
    # ensure that at least one core is un-used
    nCores <- ifelse(nCores == maxCores, nCores - 1, nCores)
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
  } else { foreach::registerDoSEQ() }
  
  
  ## loop
  tictoc::tic()
  
  Result <- data.frame()
  
  Result <- foreach::foreach(
    LoopN = LoopNr, .combine = rbind, .packages = c("sp", "dplyr", "raster")
  ) %dopar% {
    
    N <- DoubleLoop$sample_size[LoopN]
    i <- DoubleLoop$iteration[LoopN]
    
    Output <- data.frame(
      datatype = sitedat$datatype[1],
      sample_size = N, 
      n_site = 0, 
      iteration = i
    )
    
    RanNum <- sample(ids, N, replace = FALSE)
    Selected <- sitedat[sitedat$bird_id %in% RanNum, ]
    
    Output$n_site <- n_distinct(Selected$SitRecID)
    
    return(Output)
  }
  tictoc::toc()
  
  if (nCores > 1) {on.exit(parallel::stopCluster(cl))} # stop the cluster
  
  ## summarize results
  ressumm <- Result %>% group_by(datatype, sample_size) %>% 
    summarise(
      mn_nsite = mean(n_site),
      sd       = sd(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
    
  listedResults[[i]] <- list(ressumm, Result)
  
}

allsumm    <- do.call(rbind, lapply(listedResults, function(x) x[[1]]))
allresults <- do.call(rbind, lapply(listedResults, function(x) x[[2]]))


##

saveRDS(listedResults, "data/analysis/site_discover_nidreloc.rds")

## plot
ggplot() + 
  geom_point(data=allresults, aes(x=sample_size, y=n_site),
             color = "grey60", size=.5, alpha=.5) +
  geom_point(data=allsumm, aes(x=sample_size, y=mn_nsite)) +
  geom_linerange(
    data=allsumm, 
    aes(x=sample_size, ymin=mn_nsite - sd, ymax=mn_nsite + sd)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    strip.text.x = element_text(size = 14)
  ) + xlab("N birds relocated") + ylab("N sites discovered") +
  facet_grid(~datatype, scales="free_x")

#
ggsave("figures/site_discover_curve_10i.png", width=13, height=5)
  

## Run exercise within latitudinal bands to show differing degrees of survey
# completeness for each datatype ----------------------------------------------

alldat$lat_group <- as.numeric(cut(alldat$latitude, 4)) # equal-length

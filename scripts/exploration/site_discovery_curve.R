## Assessing representativeness of individual relocation data for 
# identifying sites

pacman::p_load(dplyr, foreach, stringr, ggplot2)

## fxn to run iterative sub-sampling
source("C:/Users/Martim Bill/Documents/R/source_scripts/site_discover_curve.R")

## Load combined data, with both IBA overlay info and hexgrid outsites -------
### SAVE 
alldat <- readRDS(
  "data/analysis/combined/alldatatypes_ibas_stpovrs_outsites_all_cntry.rds"
)

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
  
  ## loop
  tictoc::tic()
  
  Result <- site_discover(
    id_col = sitedat$bird_id,
    site_col = sitedat$SitRecID,
    iterations = 10, nCores = 10
  )
  
  tictoc::toc()
  
  Result$datatype <- dtype
  
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

##
saveRDS(listedResults, "data/analysis/site_discover_nidreloc.rds")

## split
# listedResults <- readRDS("data/analysis/site_discover_nidreloc.rds")

allsumm    <- do.call(rbind, lapply(listedResults, function(x) x[[1]]))
allresults <- do.call(rbind, lapply(listedResults, function(x) x[[2]]))


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


# -----------------------------------------------------------------------------
## Run exercise within latitudinal bands to show differing degrees of survey
# completeness for each datatype ----------------------------------------------
# -----------------------------------------------------------------------------

## cut data into latitudinal bands
alldat$lat_group <- as.numeric(cut(alldat$latitude, 4)) # equal-length

iterations <- 50
nCores <- 10

dtypes <- c("metal", "color", "tracking")

listedResults <- list()

## loop through each datatype
for(i in seq_along(dtypes)){
  dtype <- dtypes[i]
  print(dtype)
  
  ## filter to one datatype
  sitedat <- subset(alldat, str_detect(alldat$datatype, dtype))
  
  sitedat %>% group_by(bird_id, lat_group) %>% 
    summarise(n_site = n_distinct(SitRecID)) %>% 
    ungroup() %>% 
    summarise(
      mn_nsite = mean(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  latgrouplist <- alist()
  
  ## loop through each latitude band
  for(k in sort(unique(sitedat$lat_group))){
    print(k)
    
    ## filter to one lat band
    latdat <- subset(sitedat, sitedat$lat_group == k)
    
    ## loop
    tictoc::tic()
    
    Result_lat <- site_discover(
      id_col = latdat$bird_id,
      site_col = latdat$SitRecID,
      iterations = iterations, nCores = nCores,
      logSeg = TRUE, # even steps on a log scale
      nSteps = 100
    )
    
    tictoc::toc()
    
    Result_lat$lat_group <- k
    
    ## summarize results PER lat group
    ressumm_lat <- Result_lat %>% group_by(lat_group, sample_size) %>% 
      summarise(
        mn_nsite = mean(n_site),
        sd       = sd(n_site),
        min_nsite = min(n_site),
        max_nsite = max(n_site)
      )
    
    latgrouplist[[k]] <- list(ressumm_lat, Result_lat)
  }
  
  ## split
  ressumm_lat <- do.call(rbind, lapply(latgrouplist, function(x) x[[1]]))
  results_lat <- do.call(rbind, lapply(latgrouplist, function(x) x[[2]]))
  
  ressumm_lat$datatype <- dtype
  results_lat$datatype <- dtype
  
  ## summarize results ACROSS lat groups
  ressumm <- results_lat %>% group_by(datatype, sample_size) %>% 
    summarise(
      mn_nsite = mean(n_site),
      sd       = sd(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  listedResults[[i]] <- list(ressumm, ressumm_lat, results_lat)
  
}

## SAVE
saveRDS(listedResults, "data/analysis/site_discover_nidreloc_latgroups.rds")

# listedResults <- readRDS("data/analysis/site_discover_nidreloc_latgroups.rds")


## split
allsumm        <- do.call(rbind, lapply(listedResults, function(x) x[[1]]))
allsumm_lat    <- do.call(rbind, lapply(listedResults, function(x) x[[2]]))
allresults_lat <- do.call(rbind, lapply(listedResults, function(x) x[[3]]))


## plot

allsumm_lat$dtype_lat    <- paste(allsumm_lat$datatype, allsumm_lat$lat_group)
allresults_lat$dtype_lat <- paste(allresults_lat$datatype, allresults_lat$lat_group)


ggplot() + 
  geom_point(data=allresults_lat, aes(x=sample_size, y=n_site),
             color = "grey60", size=.5, alpha=.5) +
  geom_point(data=allsumm_lat, aes(x=sample_size, y=mn_nsite)) +
  geom_linerange(
    data=allsumm_lat, 
    aes(x=sample_size, ymin=mn_nsite - sd, ymax=mn_nsite + sd)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    strip.text.x = element_text(size = 14)
  ) + xlab("N birds relocated") + ylab("N sites discovered") +
  facet_wrap(~dtype_lat, scales="free_x")

#
ggsave("figures/site_discover_curve_latgroup_10i.png", width=14, height=12)


## split plots by lat group

latgrps_summ <- split(allsumm_lat, allsumm_lat$lat_group)
latgrps_res  <- split(allresults_lat, allresults_lat$lat_group)

nlatgrp <- n_distinct(unique(allsumm_lat$lat_group))

plotlist <- list()

for(i in seq_len(nlatgrp)){
  onesumm <- latgrps_summ[[i]]
  oneres  <- latgrps_res[[i]]
  
  latplot <- ggplot() + 
    # geom_point(data=oneres, aes(x=sample_size, y=n_site),
    #            color = "grey60", size=.5, alpha=.5) +
    geom_point(data=onesumm, aes(x=sample_size, y=mn_nsite)) +
    geom_point(data=onesumm, aes(x=log(sample_size), y=log(mn_nsite))) +
    geom_linerange(
      data=onesumm, 
      aes(x=sample_size, ymin=mn_nsite - sd, ymax=mn_nsite + sd)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      axis.text = element_text(size=14),
      axis.title = element_blank(),
      strip.text.x = element_text(size = 14)
    ) +
    facet_wrap(~dtype_lat, scales="free_x")
  
  plotlist[[i]] <- latplot
}


## patch together

library(patchwork)

source("C:/Users/Martim Bill/Documents/R/misc.scripts/add_global_label_mingsun.R")


patched <- plotlist[[4]] /  plotlist[[3]] /  plotlist[[2]] /  plotlist[[1]] & ylab("POOP")

## use hack fxn to combine axis labels
patched %>%
  add_global_label(
    Ylab = "N sites discovered",
    Xlab = "N birds relocated",
    size = 6
    # Ygap = 0.04
  )

ggsave("figures/site_discover_curve_latgroup_10iX.png", width=11, height=13)


### transform data to put it all on linear scale and get slopes

allsumm_lat$lat_group <- as.factor(allsumm_lat$lat_group)
allsumm_lat$datatype <- factor(allsumm_lat$datatype, levels = c("metal", "color", "tracking"))
fctlbls <- as_labeller(c("1" = "8°N – 23°N", "2" = "23°N – 38°N", "3" = "38°N – 52°N", "4" = "52° – 67°N"))
  
## same turbo colors as site coincidence maps/plot
turboscl <- viridis::viridis(n = 6, option = "H")

ggplot() + 
  # geom_point(data=allresults_lat, aes(x=sample_size, y=n_site),
  #            color = "grey60", size=.5, alpha=.5) +
  geom_point(data=allsumm_lat, aes(x=log(sample_size), y=log(mn_nsite), color = datatype)) +
  # geom_linerange(
  #   data=allsumm_lat, 
  #   aes(x=sample_size, ymin=mn_nsite - sd, ymax=mn_nsite + sd)) +
  # scale_color_manual(
  #   values=c("metal"="red", "color"="yellow", "tracking"= "blue"),
  #   labels=c("metal"="Metal", "color"="Color", "tracking"= "Tracking")
  # ) +
  scale_color_manual(
    values=c("metal"=turboscl[6], "color"=turboscl[4], "tracking"= turboscl[2]),
    labels=c("metal"="Metal", "color"="Color", "tracking"= "Tracking")
  ) +
  # viridis::scale_color_viridis( 
  #   option = "turbo", discrete = T
  # ) +
  scale_y_continuous(
    breaks = log(c(1.01, 10, 50, 150, 400)),
    labels = c(1, 10, 50, 150, 400)
    ) +
  scale_x_continuous(
    breaks = log(c(1.01, 10, 100, 1000, 9000)),
    labels = c(1, 10, 100, 1000, 9000)
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    strip.text.x = element_text(size = 14)
  ) + xlab("N birds relocated") + ylab("N sites discovered") +
  # facet_wrap(~forcats::fct_rev(lat_group), ncol=1) +
  facet_wrap(~lat_group,
             nrow=1,
             labeller = fctlbls)

# ggsave("figures/site_discover_curve_latgroup_logscale.png", width=5, height=13)
ggsave("figures/site_discover_curve_latgroup_logscaleX.png", height=4, width=15)



# -----------------------------------------------------------------------------
## Run exercise at country level to show how complete picture is of site used 
# within a country ------------------------------------------------------------ 
# -----------------------------------------------------------------------------

iterations <- 10
nCores <- 10

dtypes <- c("metal", "color", "tracking")

listedResults <- list()

## loop through each datatype
for(i in seq_along(dtypes)){
  dtype <- dtypes[i]
  print(dtype)
  
  ## filter to one datatype
  sitedat <- subset(alldat, str_detect(alldat$datatype, dtype))
  
  sitedat %>% group_by(bird_id, country) %>% 
    summarise(n_site = n_distinct(SitRecID)) %>% 
    ungroup() %>% 
    summarise(
      mn_nsite = mean(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  cntrygrouplist <- alist()
  
  ## loop through each latitude band
  for(k in sort(unique(sitedat$country))){
    print(k)
    
    ## filter to one lat band
    cntrydat <- subset(sitedat, sitedat$country == k)
    
    nids <- n_distinct(cntrydat$bird_id)
    
    if(nids < 4) next
    
    ## loop
    tictoc::tic()
    
    Result_cntry <- site_discover(
      id_col   = cntrydat$bird_id,
      site_col = cntrydat$SitRecID,
      iterations = iterations, nCores = nCores,
      logSeg = TRUE, # even steps on a log scale
      nSteps = 100
    )
    
    tictoc::toc()
    
    Result_cntry$country <- cntrydat$country[1]
    
    ## summarize results PER lat group
    ressumm_cntry <- Result_cntry %>% group_by(country, sample_size) %>% 
      summarise(
        mn_nsite = mean(n_site),
        sd       = sd(n_site),
        min_nsite = min(n_site),
        max_nsite = max(n_site)
      )
    
    cntrygrouplist[[k]] <- list(ressumm_cntry, Result_cntry)
  }
  
  ## split
  ressumm_cntry <- do.call(rbind, lapply(cntrygrouplist, function(x) x[[1]]))
  results_cntry <- do.call(rbind, lapply(cntrygrouplist, function(x) x[[2]]))
  
  ressumm_cntry$datatype <- dtype
  results_cntry$datatype <- dtype
  
  ## summarize results ACROSS lat groups
  ressumm <- results_cntry %>% group_by(datatype, sample_size) %>% 
    summarise(
      mn_nsite = mean(n_site),
      sd       = sd(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  listedResults[[i]] <- list(ressumm, ressumm_cntry, results_cntry)
  
}



## split
allsumm          <- do.call(rbind, lapply(listedResults, function(x) x[[1]]))
allsumm_cntry    <- do.call(rbind, lapply(listedResults, function(x) x[[2]]))
allresults_cntry <- do.call(rbind, lapply(listedResults, function(x) x[[3]]))


## plot

allsumm_cntry$dtype_cntry    <- paste(allsumm_cntry$datatype, allsumm_cntry$country)
allresults_cntry$dtype_cntry <- paste(allresults_cntry$datatype, allresults_cntry$country)

## remove countries w/ fewer than 10 inds

cntrytoplot <- allsumm_cntry %>% group_by(country) %>% 
  summarise(
    max_ss = max(sample_size)
  ) %>% filter(max_ss >= 10)

allsumm_plot <- filter(allsumm_cntry, country %in% cntrytoplot$country)

## same turbo colors as site coincidence maps/plot
turboscl <- viridis::viridis(n = 6, option = "H")

ggplot() + 
  # geom_point(data=allresults_cntry, aes(x=sample_size, y=n_site),
  #            color = "grey60", size=.5, alpha=.5) +
  geom_point(data=allsumm_plot, aes(x=sample_size, y=mn_nsite, color=datatype)) +
  scale_color_manual(
    values=c("metal"=turboscl[6], "color"=turboscl[4], "tracking"= turboscl[2]),
    labels=c("metal"="Metal", "color"="Color", "tracking"= "Tracking")
  ) +  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    strip.text.x = element_text(size = 14)
  ) + 
  xlab("N birds relocated") + ylab("N sites discovered") +
  facet_wrap(~country, scales="free")

## save
ggsave("figures/site_discover/site_discover_curve_incntrysites.png", 
       height=8, width=12)

# -----------------------------------------------------------------------------
## Assess how complete picutre is of site use for any birds that pass thru
# a country AT ANY GIVEN TIME (e.g. PT's birds in PT and elsewhere) -----------  
# -----------------------------------------------------------------------------

iterations <- 10
nCores <- 10

dtypes <- c("metal", "color", "tracking")

listedResults <- list()

## loop through each datatype
for(i in seq_along(dtypes)){
  dtype <- dtypes[i]
  print(dtype)
  
  ## filter to one datatype
  sitedat <- subset(alldat, str_detect(alldat$datatype, dtype))
  
  sitedat %>% group_by(bird_id, country) %>% 
    summarise(n_site = n_distinct(SitRecID)) %>% 
    ungroup() %>% 
    summarise(
      mn_nsite = mean(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  cntrygrouplist <- alist()
  
  ## loop through each latitude band
  for(k in sort(unique(sitedat$country))){
    print(k)
    
    idsincntry <- sitedat %>% group_by(bird_id, country) %>% 
      summarise() %>% filter(country == k)
    
    ## filter to one lat band
    cntrydat <- subset(sitedat, bird_id %in% unique(idsincntry$bird_id))
    
    nids <- n_distinct(cntrydat$bird_id)
    
    if(nids < 4) next
    
    ## loop
    tictoc::tic()
    
    Result_cntry <- site_discover(
      id_col   = cntrydat$bird_id,
      site_col = cntrydat$SitRecID,
      iterations = iterations, nCores = nCores,
      logSeg = TRUE, # even steps on a log scale
      nSteps = 100
    )
    
    tictoc::toc()
    
    Result_cntry$country <- cntrydat$country[1]
    
    ## summarize results PER lat group
    ressumm_cntry <- Result_cntry %>% group_by(sample_size) %>% 
      summarise(
        country  = k,
        mn_nsite = mean(n_site),
        sd       = sd(n_site),
        min_nsite = min(n_site),
        max_nsite = max(n_site)
      )
    
    cntrygrouplist[[k]] <- list(ressumm_cntry, Result_cntry)
  }
  
  ## split
  ressumm_cntry <- do.call(rbind, lapply(cntrygrouplist, function(x) x[[1]]))
  results_cntry <- do.call(rbind, lapply(cntrygrouplist, function(x) x[[2]]))
  
  ressumm_cntry$datatype <- dtype
  results_cntry$datatype <- dtype
  
  ## summarize results ACROSS lat groups
  ressumm <- results_cntry %>% group_by(datatype, sample_size) %>% 
    summarise(
      mn_nsite = mean(n_site),
      sd       = sd(n_site),
      min_nsite = min(n_site),
      max_nsite = max(n_site)
    )
  
  listedResults[[i]] <- list(ressumm, ressumm_cntry, results_cntry)
  
}

## split
allsumm          <- do.call(rbind, lapply(listedResults, function(x) x[[1]]))
allsumm_cntry2    <- do.call(rbind, lapply(listedResults, function(x) x[[2]]))
allresults_cntry2 <- do.call(rbind, lapply(listedResults, function(x) x[[3]]))


## plot

allsumm_cntry2$dtype_cntry    <- paste(allsumm_cntry2$datatype, allsumm_cntry2$country)
allresults_cntry2$dtype_cntry <- paste(allresults_cntry2$datatype, allresults_cntry2$country)

## remove countries w/ fewer than 10 inds

cntrytoplot2 <- allsumm_cntry2 %>% group_by(country) %>% 
  summarise(
    max_ss = max(sample_size)
  ) %>% filter(max_ss >= 10)

allsumm_plot2 <- filter(allsumm_cntry2, country %in% cntrytoplot2$country)

ggplot() + 
  # geom_point(data=allresults_cntry, aes(x=sample_size, y=n_site),
  #            color = "grey60", size=.5, alpha=.5) +
  geom_point(data=allsumm_plot2, aes(x=sample_size, y=mn_nsite, color=datatype)) +
  scale_color_manual(
    values=c("metal"=turboscl[6], "color"=turboscl[4], "tracking"= turboscl[2]),
    labels=c("metal"="Metal", "color"="Color", "tracking"= "Tracking")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    strip.text.x = element_text(size = 14)
  ) + 
  xlab("N birds relocated") + ylab("N sites discovered") +
  facet_wrap(~country, scales="free")

## save
ggsave("figures/site_discover/site_discover_curve_cntrypop.png", height=8, width=12)

## Log scale ## 

ggplot() + 
  geom_point(
    data=allsumm_plot2, 
    aes(x=log(sample_size), y=log(mn_nsite), 
        color = datatype)) +
  scale_color_manual(
    values=c("metal"=turboscl[6], "color"=turboscl[4], "tracking"= turboscl[2]),
    labels=c("metal"="Metal", "color"="Color", "tracking"= "Tracking")
  ) +
  scale_y_continuous(
    breaks = log(c(1.01, 10, 50, 150, 400)),
    labels = c(1, 10, 50, 150, 400)
  ) +
  scale_x_continuous(
    breaks = log(c(1.01, 10, 100, 1000, 9000)),
    labels = c(1, 10, 100, 1000, 9000)
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    strip.text.x = element_text(size = 14)
  ) + xlab("N birds relocated") + ylab("N sites discovered") +
  # facet_wrap(~forcats::fct_rev(lat_group), ncol=1) +
  facet_wrap(~country)


ggsave("figures/site_discover/site_discover_curve_cntrypop_logscaleX.png", 
       height=8, width=12)

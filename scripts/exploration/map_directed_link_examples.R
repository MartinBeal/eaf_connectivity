## run first part of map_direct_networks first

site <- filter(site_summ, site_poly == "Tejo estuary")

siteid <- site$SitRecID

node_plot  <- filter(node_prj, loc_num == siteid)
edge_plot <- filter(edge_prj, from_sid == siteid)

### Map 1a: edges colored by data origin, nodes by disconnectedness ---------
map5 <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
  # geom_sf(data = edge_prj,                 ## edges
  #         aes(size = n_bird_trax), col = "black") + #, alpha = 0.65
  geom_sf(data = arrange(edge_plot, prop_id_out_from),                 ## edges
          aes(color = prop_id_out_from)) + #, alpha = 0.65
  geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
  geom_sf(data = arrange(node_plot, n_id_out), 
          aes(fill = n_id_out),
          pch=21, stroke=.25,
          col="black") +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  # scale_fill_brewer(type="qual") +
  # labs(color = "Data type",
  #      fill = "Connected") +
  # scale_size(range = c(0.01, 2)) +
  maptheme
map5

ggsave(
  plot = map5, 
  paste0("figures/direct_networks/out_links_", site$site_poly, "_", season, ".png"), 
  width=5, height = 5)

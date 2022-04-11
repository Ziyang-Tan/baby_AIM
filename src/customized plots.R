library(ggplot2)

labeled_line_graph <- function(data, aes_line, ref, padding = 0.2){
  x_name = quo_name(aes_line$x)
  pad = max(data[x_name]) * padding
  col = names(ref)
  ref_name = ref[[col]]
  aes_label = aes_line
  aes_label$x <- aes(x = max(data[x_name]) + pad/2)[[1]]
  ggplot(data = data %>% filter(!!sym(col)!=ref_name), 
         mapping = aes_line) +
    geom_line() + 
    geom_point() +
    scale_colour_manual(values=wes_palette("FantasticFox1", 9, type = 'continuous')) +
    geom_label_repel(mapping = aes_line, 
                     nudge_x = -pad, 
                     na.rm = TRUE, 
                     direction = 'y',
                     segment.linetype = 6) +
    xlim(min(data[x_name])-pad, max(data[x_name])+2*pad) +
    geom_point(data = data %>% filter(!!sym(col)==ref_name), 
               mapping = aes_label, 
               shape='M', 
               size=5) +
    geom_label_repel(data = tmp %>% filter(!!sym(col)==ref_name), 
                     mapping = aes_label, 
                     nudge_x = max(data[x_name])+pad, 
                     na.rm = TRUE, 
                     direction = 'y',
                     segment.linetype = 6) +
    theme(legend.position = "none")
}




plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce') # for `geom_sina`
  plot1 <- ggplot(data = data_long)+
    coord_flip() +
    labs(title  = "SHapley Additive exPlanations", x = "",
         y = 'SHAP value (Impact on Model Output)',
         color = "Feature Value")+
    # sina plot:
    geom_sina(aes(x = variable, y = value, color = stdfvalue),size=1.7) +
    # print the mean absolute value:
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.4f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2,
              fontface = "bold") + # bold
    #theme_classic(base_rect_size = 1.5) +
	theme_bw(base_rect_size = 1.65) +
    guides(fill = guide_legend(title.position = 'right'))+
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) +
    # scale_color_gradient(low="#FFCC33", high="#6600CC",
    #                      breaks=c(0,1), labels=c("Low","High")) +
    scale_color_gradient2(low="#5AB4C2",mid = 'white', high="#E36258",
                          labels=c("Low",'',"High"),breaks=c(0,0.5,1),
                          midpoint = 0.5 #breaks=c(0,1), labels=c("Low",'0',"High")
                         ) +

    theme(axis.title.x = element_text(size = 10,colour = 'black',face = 'bold',vjust= .01),
          axis.title.y = element_text(size = 12,colour = 'grey15',face = 'bold',vjust= 2),
          plot.title = element_text(size = 12,colour = 'darkred',face = 'bold',hjust= .6),
          axis.text = element_text(size= 9 ,face = 'bold'),
          axis.ticks = element_blank(),

          axis.line.y = element_blank(),
          #axis.line.x = element_line( color = 'black',size = .8 ),

          panel.grid.major  = element_line(color = "white", size = .8, linetype = "solid"),
          panel.grid.minor = element_line(color = "white", size = .5, linetype = "solid"),
          
		  #panel.background = element_rect(fill='white'),
      panel.background = element_rect(  fill=  scales::alpha('#F0F8FF',.7) ),

          legend.key.height = unit(27,'mm'),
          legend.key.width  = unit(1.3,'mm'),
          legend.position="right") +
    geom_hline(yintercept = 0,color = scales::alpha('grey80',1),size=.65,linetype = "dashed") + # the vertical line
  #geom_hline(yintercept = 0,color = scales::alpha('white',1),size= 1.0,linetype = "solid") + # the vertical line
  
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable))
    )

  return(plot1)
}

?guide_legend


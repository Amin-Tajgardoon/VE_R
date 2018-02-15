## function to draw calibration curve

library(ggplot2)

save.calib.curve = function(yt, pt, nBin, title){
  
  calib = calib.xy(yt,pt, nBin)
  d = data.frame(obs = calib$observed,exp = calib$expected)

  max.limit = max(c(d$obs,d$exp), na.rm = TRUE)
  g = ggplot(data = d, aes( obs, exp)) + xlim(c(0,max.limit)) + ylim(c(0,max.limit))
  g = g + geom_line(colour = "blue", lwd = 1)
  g = g+ geom_point(size = 3) #+ geom_smooth(span = .8,se = FALSE)
  g = g + geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "longdash", lwd = 1)
  g = g + labs(x = "Fraction of observed positive",
             y = "Mean of predicted positive",
             title = title)
  g = g+theme(plot.title=element_text(face="bold",hjust=.5))
  g = g + coord_equal()
  plot(g)

  ggsave(filename = paste0(title,".png"), plot = g, device = "png", dpi = 300)
  return(g)

}

 #d = read.csv("test_yt_pt.csv")
 #yt = d$yt
 #pt = d$pt
 #nBin = 10
 #title = "Reliability digaram" 
 #g = save.calib.curve(yt = yt, pt = pt, nBin = 10, title )                     
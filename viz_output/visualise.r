require(ggplot2)
require(reshape2)
library(directlabels)
library(tools)
library(scatterplot3d)


options(echo=FALSE)

args <- commandArgs(trailingOnly = TRUE)

results <- read.csv(args[1])

png(paste("./output/",file_path_sans_ext(basename(args[1])),".png", sep=''), width=800, height=800)

# with(results, {
#     scatterplot2d(d1,   # x axis
#                   d2,     # y axis
#                   main=args[2])
# })


ggplot(results, aes(x=d1, y=d2)) + geom_point()

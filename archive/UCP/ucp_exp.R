## Heatmap of average expression.
firstCol <- function(df) {
    df <- df[ ,1:2]
    # colnames(df) <- c("Samples", funcName(df))
    return(df)
}

funcName <- function(func) {
    deparse(substitute(func))
}

df <- read.csv("UCP_average.csv")
# df$Samples <- with(df, reorder(Samples, UCP1))
library(ggplot2)
df.m <- melt(df)
df.m <- ddply(df.m, .(variable), transform, rescale=scale(value))

p <- ggplot(df.m, aes(variable, Samples)) + geom_tile(aes(fill = rescale),
     color = "white") + scale_fill_gradient(low = "white",
     high = "steelblue")

base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "",
    y = "") + scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none",
    axis.text.x = element_text(size = base_size * 0.8, angle = 330,
    hjust = 0, colour = "grey50"))
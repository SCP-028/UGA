library(dplyr)
setwd("~/data/serinesynthesis")
load("../Cancers_Express.RData")
annot <- data.frame(
  `shared name`=c("3-phospho-D-glycerate (PHGDH) 3-phosphonooxypyruvate",
                  "3-phosphonooxypyruvate (PSAT1) phosphoserine",
                  "phosphoserine (PSPH) serine"),
  source=c("3-phospho-D-glycerate", "3-phosphonooxypyruvate", "phosphoserine"),
  interaction=c("PHGDH", "PSAT1", "PSPH"),
  target=c("3-phosphonooxypyruvate", "phosphoserine", "serine"),
  check.names=F  # prevent changing `shared name` to `shared.name`
)
result <- vector("list", length(All.DE))
for(i in seq_along(All.DE)) {
  df <- as.data.frame(All.DE[[i]])
  df <- df[grep("(PHGDH)|(PSAT1$)|(PSPH$)", rownames(df)), ]
  df$interaction <- sub("^.*\\|(.*)$", "\\1", rownames(df))
  df$colors <- ifelse(df$p.fdr <= 0.05, ifelse(df$fold.change >= 1, "red", "blue"), "black")
  result[[i]] <- df
}
names(result) <- names(All.DE)
for(i in seq_along(result)) {
  df <- result[[i]]
  df <- left_join(df, annot)
  write.csv(df, file = paste0(names(result)[i], ".csv"), quote = F, row.names = F)
}

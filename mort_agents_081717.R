############################
#  Mortality agent prevalence
############################
library("VennDiagram")
library("vioplot")
opar <- par(no.readonly = TRUE)


plot_vars <- read.csv("all_vars_EXPORT.csv")
greenwood_trees <- read.csv("Individual_TreeData_AllData_SF_edits.csv")
  greenwood_trees$Unique_tree <- paste0(greenwood_trees$Plot, greenwood_trees$Number)
trees <- read.csv("all_trees_with_delta_and_ENN_041916.csv")



## Calculate infestation by particular mortality agents


### some numbers for juniper
  
nrow(trees[trees$Spp == "JUOS" & trees$mt == "Y", ]) #juos with mt

nrow(trees[trees$Spp == "JUOS" & trees$bb == "Y", ]) #juos with bark beetles

nrow(trees[trees$Spp == "JUOS", ]) #total juos

####################
# Mean canopy decline for different mort agent classes
####################

#no mort agents

#create convenient column for NoMort
trees$NoMort <- "N"

trees[which(trees$Spp == "PIMO"
            & trees$Live == "Y"
            & trees$ips == "N"
            & trees$rtb == "N"
            & trees$pmb == "N"
            & trees$mt == "N"
            & trees$ptb == "N"
            & trees$ptm == "N"
            & trees$pns == "N"
            & trees$psf == "N"), 
            "NoMort"] <- "Y"
  
mean(trees$Delta_pdc[trees$Spp == "PIMO"
              & trees$Live == "Y"
              & trees$ips == "N"
              & trees$rtb == "N"
              & trees$pmb == "N"
              & trees$mt == "N"
              & trees$ptb == "N"
              & trees$ptm == "N"
              & trees$pns == "N"
              & trees$psf == "N"])

sd(trees$Delta_pdc[trees$Spp == "PIMO"
                     & trees$Live == "Y"
                     & trees$ips == "N"
                     & trees$rtb == "N"
                     & trees$pmb == "N"
                     & trees$mt == "N"
                     & trees$ptb == "N"
                     & trees$ptm == "N"
                     & trees$pns == "N"
                     & trees$psf == "N"]) / 
  sqrt(length(trees$Delta_pdc[trees$Spp == "PIMO"
                         & trees$Live == "Y"
                         & trees$ips == "N"
                         & trees$rtb == "N"
                         & trees$pmb == "N"
                         & trees$mt == "N"
                         & trees$ptb == "N"
                         & trees$ptm == "N"
                         & trees$pns == "N"
                         & trees$psf == "N"]))

#bark beetles
trees$Beetles <- "N"

trees[which(trees$Spp == "PIMO"
            & (trees$ips == "Y"
            | trees$rtb == "Y")), 
      "Beetles"] <- "Y"


mean(trees$Delta_pdc[trees$Spp == "PIMO" 
                    & trees$ips == "Y"
                     | trees$rtb == "Y"])

sd(trees$Delta_pdc[trees$Spp == "PIMO" 
                     & trees$ips == "Y"
                     | trees$rtb == "Y"]) / sqrt(length(trees$Delta_pdc[trees$Spp == "PIMO" 
                                                                        & trees$ips == "Y"
                                                                        | trees$rtb == "Y"]))

## folivores
trees$Folivores <- "N"

trees[which(trees$Spp == "PIMO"
            & (trees$ptm == "Y"
               | trees$pns == "Y"
               | trees$psf == "Y")), 
      "Folivores"] <- "Y"

mean(trees$Delta_pdc[trees$Spp == "PIMO" 
                     & trees$ptm == "Y"
                     | trees$pns == "Y"
                     | trees$psf == "Y"])

sd(trees$Delta_pdc[trees$Spp == "PIMO" 
                   & trees$ptm == "Y"
                   | trees$pns == "Y"
                   | trees$psf == "Y"]) / sqrt(length(trees$Delta_pdc[trees$Spp == "PIMO" 
                                                                      & trees$ptm == "Y"
                                                                      | trees$pns == "Y"
                                                                      | trees$psf == "Y"]))

NoMort_t_test <- t.test(Delta_pdc ~ NoMort, data = trees)

NoMort_diff_mixed <- lmer(Delta_pdc ~ NoMort + (1|Plot), data = trees)
summary(NoMort_diff_mixed)

#####
  # Four set Venn Diagram

A <- which((trees$Spp == "PIMO" & trees$Live  == "Y" & (trees$ips == "Y" | trees$rtb == "Y" | trees$pmb == "Y")))
B <- which((trees$Spp == "PIMO" & trees$Live  == "Y" & (trees$ptb == "Y" | trees$ptm == "Y")))
C <- which((trees$Spp == "PIMO" & trees$Live  == "Y" & (trees$pns == "Y" | trees$psf == "Y")))
D <- which((trees$Spp == "PIMO" & trees$Live  == "Y" & (trees$mt == "Y")))
# E <- sample(1:1000, 375, replace = FALSE);

#number of trees with no mortality agents
nrow(trees[trees$Spp =="PIMO" & trees$Live  == "Y"& trees$ips == "N" & trees$rtb == "N" & trees$pmb == "N" & trees$ptb == "N"
           & trees$ptm == "N" & trees$pns == "N" & trees$psf == "N" & trees$mt == "N", ])

nrow(trees[trees$Spp == "PIMO" & trees$Live  == "Y", ])

length(A)
length(B)
length(C)
length(D)
venn.plot <- venn.diagram(
  x = list(
    "Bark (n = 903)" = A,
    "MT (n = 320)" = D,
    "Twig (n = 1173)" = B,
    "Foli (n = 1178)"= C
  ),
  filename = "mort_agent_venn_111916.tiff",
  col = "transparent",
  fill = c("grey30", "grey50", "grey65", "grey85"),
  alpha = 0.4,
  label.col = c("grey30", "white", "grey30", "white",
                "white", "white", "white", "white", "grey30", "white",
                "white", "white", "white", "grey30", "white"),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.09,
  cat.fontfamily = "serif",
  rotation.degree = 270,
  margin = 0,
  height = 6,
  width = 6,
  resolution = 600,
  imagetype = "tiff",
  units = "in",
  compression = "lzw",
  antialias = "gray"
);


#######################################################
# Pie chart of dead trees ~ ips
par <- opar
par(mfrow = c(2,1))
par(mar = c(1,1,1,1))
pie(c(nrow(trees[trees$Live == "Y" & trees$Spp == "PIMO" & trees$ips == "Y", ]), 
      nrow(trees[trees$Live == "Y" & trees$Spp == "PIMO" & trees$ips == "N",])),
    c("Ips", "No Ips"), radius = 1
)
text(x = 1.5, y = 1, "Live Trees", cex = 1.3)
pie(c(nrow(trees[trees$Live == "N" & trees$Spp == "PIMO" & trees$ips == "Y", ]), 
      nrow(trees[trees$Live == "N" & trees$Spp == "PIMO" & trees$ips == "N",])),
    c("Ips", "No Ips"), radius = 1
    )
text(x = 1.5, y = 1, "Dead Trees", cex = 1.3)


########################################################
# Boxplots for PDC ~ mortality agents

## For this to work, we need to use the "trees <- tree_vars" section of code which uses a tree dataset with Delta_pdc
par <- opar

#applies mean to each type
mort_means <- unlist(lapply(list(trees$Delta_pdc[trees$Spp == "PIMO" & trees$ips == "Y"], 
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$rtb == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$pmb == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$mt == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$ptb == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$ptm == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$pns == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$psf == "Y"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$ips == "N"
                          & trees$Live == "N"
                          & trees$rtb == "N"
                          & trees$pmb == "N"
                          & trees$mt == "N"
                          & trees$ptb == "N"
                          & trees$ptm == "N"
                          & trees$pns == "N"
                          & trees$psf == "N"],
          trees$Delta_pdc[trees$Spp == "PIMO" & trees$ips == "N"
                          & trees$Live == "Y"
                          & trees$rtb == "N"
                          & trees$pmb == "N"
                          & trees$mt == "N"
                          & trees$ptb == "N"
                          & trees$ptm == "N"
                          & trees$pns == "N"
                          & trees$psf == "N"]), FUN = mean))

mean(trees$Delta_pdc[trees$Spp == "PIMO" & (trees$ips == "Y" | trees$rtb == "Y")])
mean(trees$Delta_pdc[trees$Spp == "PIMO" & (trees$pns == "Y" | trees$psf == "Y")])

tiff(filename="mort_agent_vioplot_111916.tiff", 
    type="cairo",
    units="in", 
    width = 6, 
    height = 4, 
    pointsize=10, 
    res=600,
    compression = "lzw",
    antialias = "gray")

par(oma = c(0,0,0,0),
    mar = c(5.1,4.1,1,1),
    family = "serif",
    bty = "n")


vioplot(trees$Delta_pdc[trees$Spp == "PIMO" & trees$ips == "Y"], 
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$rtb == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$pmb == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$mt == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$ptb == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$ptm == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$pns == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$psf == "Y"],
        trees$Delta_pdc[trees$Spp == "PIMO" & trees$ips == "N"
                        & trees$Live == "Y"
                        & trees$rtb == "N"
                        & trees$pmb == "N"
                        & trees$mt == "N"
                        & trees$ptb == "N"
                        & trees$ptm == "N"
                        & trees$pns == "N"
                        & trees$psf == "N"],
        names = c("IPS", "RTB", "PMB", "MT", "PTB", "PTM", "PNS", "PSF", ""),
        col = "grey80",
        ylim = c(-100, 105))
title(ylab = "Change in percent live crown", xlab = "Mortality agent")
text(x = c(1:9), y = 105, c("597", "792", "109", "453", "1167", "130", "733", "619", "162"))
axis(1, at = 9, labels = "None,\nLive", line = 0.7, tick = FALSE)
points(x = c(1:9), y = mort_means[c(1:8,10)], pch = 24, bg = "white")
abline(h = median(trees$Delta_pdc[trees$Spp == "PIMO"]))
abline(h =mean((trees$Delta_pdc[trees$Spp == "PIMO"])),
       lty = 2)

dev.off()
##############################
# Percentage infestation 2005 vs 2015

comp_2015 <- c(length(which(trees$ips == "Y" & trees$Live == "Y")) / length(which(trees$Spp == "PIMO"& trees$Live == "Y")),
               length(which(trees$rtb == "Y" & trees$Live == "Y")) / length(which(trees$Spp == "PIMO"& trees$Live == "Y")),
               length(which(trees$mt == "Y" & trees$Spp == "PIMO")) / length(which(trees$Spp == "PIMO"& trees$Live == "Y")))

comp2005 <- c(length(which(greenwood_trees$MortalityCause1 == "IPS" & greenwood_trees$LiveDead == "L")) / 
                length(which(greenwood_trees$Spp == "Pimo"& greenwood_trees$LiveDead == "L")),
              length(which(greenwood_trees$MortalityCause2 == "RWB" & greenwood_trees$LiveDead == "L")) / 
                length(which(greenwood_trees$Spp == "Pimo"& greenwood_trees$LiveDead == "L")),
              length(which(greenwood_trees$MortalityCause3 == "MT" & greenwood_trees$LiveDead == "L")) / 
                length(which(greenwood_trees$Spp == "Pimo"& greenwood_trees$LiveDead == "L")))

#############################
## Infestation 2005 vs mortality risk

trees$Died <- "N"
for(i in 1:nrow(trees)){
  if (trees$Live[i] == "N" & trees$LiveDead[i] == "L"){
    trees$Died[i] <- "Y"
  }
}
trees$Died <- as.factor(trees$Died)

trees$Ips05 <- "N"
for (i in 1:nrow(trees)){
  trees$Ips05[i] <- ifelse(trees$Unique_tree[i] %in% greenwood_trees$Unique_tree & 
                             greenwood_trees[greenwood_trees$Unique_tree == trees$Unique_tree[i], "MortalityCause1"] == "IPS",
                           "Y", "N")
}

trees$RTB05 <- "N"
for (i in 1:nrow(trees)){
  trees$RTB05[i] <- ifelse(trees$Unique_tree[i] %in% greenwood_trees$Unique_tree & 
                             greenwood_trees[greenwood_trees$Unique_tree == trees$Unique_tree[i], "MortalityCause2"] == "RWB",
                           "Y", "N")
}

trees$MT05 <- "N"
for (i in 1:nrow(trees)){
  trees$MT05[i] <- ifelse(trees$Unique_tree[i] %in% greenwood_trees$Unique_tree & 
                            (greenwood_trees[greenwood_trees$Unique_tree == trees$Unique_tree[i], "MortalityCause3"] == "MT"
                             | greenwood_trees[greenwood_trees$Unique_tree == trees$Unique_tree[i], "DMR"] > 0),
                          "Y", "N")
}

par(mfrow = c(2,1))
pie(c(nrow(trees[trees$Died == "N" & trees$Spp == "PIMO" & trees$Ips05 == "Y", ]), 
      nrow(trees[trees$Died == "N" & trees$Spp == "PIMO" & trees$Ips05 == "N",])),
    c("Ips", "No Ips")
)
pie(c(nrow(trees[trees$Died == "Y" & trees$Spp == "PIMO" & trees$Ips05 == "Y", ]), 
      nrow(trees[trees$Died == "Y" & trees$Spp == "PIMO" & trees$Ips05 == "N",])),
    c("Ips", "No Ips")
)
####
pie(c(nrow(trees[trees$Died == "N" & trees$Spp == "PIMO" & trees$RTB05 == "Y", ]), 
      nrow(trees[trees$Died == "N" & trees$Spp == "PIMO" & trees$RTB05 == "N",])),
    c("RTB", "No RTB")
)
pie(c(nrow(trees[trees$Died == "Y" & trees$Spp == "PIMO" & trees$RTB05 == "Y", ]), 
      nrow(trees[trees$Died == "Y" & trees$Spp == "PIMO" & trees$RTB05 == "N",])),
    c("RTB", "No RTB")
)
####
pie(c(nrow(trees[trees$Died == "N" & trees$Spp == "PIMO" & trees$MT05 == "Y", ]), 
      nrow(trees[trees$Died == "N" & trees$Spp == "PIMO" & trees$MT05 == "N",])),
    c("MT", "No MT")
)
pie(c(nrow(trees[trees$Died == "Y" & trees$Spp == "PIMO" & trees$MT05 == "Y", ]), 
      nrow(trees[trees$Died == "Y" & trees$Spp == "PIMO" & trees$MT05 == "N",])),
    c("MT", "No MT")
)

#######################
summary(lm(Delta_pdc ~ ips + rtb + pmb + mt + 
               ptb + ptm + pns + psf, data = trees[trees$Spp == "PIMO" & trees$Live == "Y", ]))


summary(lmer(Delta_pdc ~ ips + rtb + pmb + mt + 
     ptb + ptm + pns + psf + (1|Plot), data = trees[trees$Spp == "PIMO" & trees$Live == "Y", ]))

summary(glmer(Died ~ ips + rtb + pmb + (1|Plot), 
             data = trees[trees$Spp == "PIMO", ],
             family = "binomial"))

summary(glmer(Died ~ ips + (1|Plot), 
              data = trees[trees$Spp == "PIMO", ],
              family = "binomial"))


## compare different agents versus no agent
ips + rtb + pmb + mt + 
  ptb + ptm + pns + psf


t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$ips == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$rtb == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$pmb == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$mt == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$ptb == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$ptm == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$pns == "Y", "Delta_pdc"])

t.test(trees[trees$Spp == "PIMO" & trees$NoMort == "Y", "Delta_pdc"], 
       trees[trees$Spp == "PIMO" & trees$psf == "Y", "Delta_pdc"])

# #########################
# # 3 set scaled Euler plot
# 
# 
# A <- which((trees$ips == "Y" | trees$rtb == "Y" | trees$pmb == "Y"))
# B <- which((trees$ptb == "Y" | trees$ptm == "Y" | trees$pns == "Y" | trees$psf == "Y"))
# C <- which((trees$mt == "Y" | trees$DMR > 0))
# 
# nab <- intersect(A, B)
# nbc <- intersect(B, C)
# nac <- intersect(A, C)
# nabc <- intersect(nab, C)
# 
# # A more complicated diagram
# venn.plot <- draw.triple.venn(
#   area1 = length(A),
#   area2 = length(B),
#   area3 = length(C),
#   n12 = length(nab),
#   n23 = length(nbc),
#   n13 = length(nac),
#   n123 = length(nabc),
#   euler.d = T,
#   overrideTriple = T,
#   scaled = T,
#   category = c("Bark", "Foli", "MT"),
#   col = "transparent",
#   fill = c("cornflowerblue", "green", "darkorchid1"),
#   alpha = 0.50,
#   label.col = c("white", "white", "white",
#                  "white", "white",
#                 "white", "white"),
#   cex = 1.5,
#   fontfamily = "serif",
#   fontface = "bold",
#   cat.col = c("darkblue", "darkgreen", "darkorchid4"),
#   cat.cex = 1.5,
#   cat.pos = 0,
#   cat.dist = 0.07,
#   cat.fontfamily = "serif",
#   rotation.degree = 270,
#   margin = 0
# )
# grid.draw(venn.plot)
# 
# 
# 
# 
# v3 <- get.venn.partitions(x = list(A, B, C))
# 
# draw.triple.venn(
#   v3,
#   filename = "Euler_3set_scaled.tiff",
#   cex = 2.5,
#   cat.cex = 2.5,
#   cat.pos = 0,
#   scaled = T,
#   overrideTriple = T
# );
# 
# 
# 

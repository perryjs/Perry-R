require(vegan)

# kmeans
kclus <- kmeans(d, centers = 2, iter.max=1000, nstart=10000)

# distance matrix
d_dist <- dist(d)

# Multidimensional scaling
cmd <- cmdscale(d_dist)

# plot MDS, with colors by groups from kmeans
groups <- levels(factor(kclus$cluster))
fig <- ordiplot(cmd, type = "points")
cols <- c("black", "black")
for(i in seq_along(groups)){
  points(cmd[factor(kclus$cluster) == groups[i], ], col = cols[i], pch = 16)
}
cmd.spp  <- scores(cmd, display = "species", 
                    scaling = 2)
lab.xloc <- cmd.spp[,1] * 0.90
lab.yloc <- cmd.spp[,2] * 1.0
text(lab.xloc, lab.yloc, rownames(cmd.spp), col = "firebrick2", cex = 0.8, pos = 3, offset = 0.75)


# add spider and hull
ordispider(cmd, factor(kclus$cluster), label = FALSE)
ordihull(cmd, factor(kclus$cluster), lty = "solid")

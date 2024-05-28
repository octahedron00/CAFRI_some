# by Octo Moon, octahedron00@gmail.com
# wrote in 2024.02, with R 4.3.2
# HyunLab CAFRI-Viola Project

# to use library in this project only
.libPaths(paste(R.home(), "/library"))

if(!require(BiocManager)){
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}
if(!require(reshape2)){
  install.packages("reshape2", repos = "http://cran.us.r-project.org")
}
if(!require(ggtree)){
  BiocManager::install("ggtree")
}
if(!require(viridis)){
  BiocManager::install("viridis")
}




library(ggtree)
library(treeio)
library(reshape2)
library(viridis)

source("src\\R_script\\gheatmap_custom.R")

tree_file <- "var\\_tree.dnd"
df_file <- "var\\_df.csv"
md_file <- "var\\_md.csv"


COL_W = 0.7
TREE_W = 0.03
ROW_H = 0.166

MARGIN_ERROR_CORRECTION = 0.035

TEXT_W = 0.08
TEXT_SIZE = 4
# size = 4

# Text adjustment is not completed yet...
# G and . have different distances, and...
# ...Which is determined by the font using.

# All works same from here.



tree <- read.tree(tree_file)

df <- read.csv(df_file)
length_text = max(nchar(df$Gene)) + 1
rownames(df) <- df$Gene
df <- df[, -1]

df[is.na(df)] <- 0

md <- read.csv(md_file)
md_list = list()

for(row in 1:nrow(md)){

  if(is.null(md_list[[md[row, 'group']]])){
    md_list[[md[row, 'group']]] = append(c(),  md[row, 'name'])
  }
  else{
    md_list[[md[row, 'group']]] = append(md_list[[md[row, 'group']]],  md[row, 'name'])
  }

  md_list[[md[row, 'group']]] = append(md_list[[md[row, 'group']]],  md[row, 'name'])
}

tree <- groupOTU(tree, md_list)


p_tree <- ggtree(tree, branch.length = 'none') +
  geom_tiplab(aes(color=group),size=TEXT_SIZE, align=TRUE) +
  scale_color_manual(
    values=c( '0' = "black", 'origin' = "red", 'duplicated' = "#777777",
              'none' = "black", 'align' = "black", 'add' = "blue", 'align+add' = '#000077'))

df_tree <- p_tree$data
depth_tree <- max(df_tree$x, na.rm=TRUE)
print(df_tree)

plot <- new_gheatmap(p_tree, df, width = (ncol(df) * COL_W) / (nrow(df) * TREE_W),
                 offset = (length_text * TEXT_W + nrow(df) * TREE_W * MARGIN_ERROR_CORRECTION) * depth_tree / (nrow(df) * TREE_W) - 1)

# + scale_x_ggtree()

plot <- plot + scale_fill_viridis(option="viridis", discrete=FALSE)

# 300 dpi, inch unit is working...
# +1 for the ending
# +1

ggsave("result.png", plot, limitsize = FALSE,
       height = max(nrow(df), 16) * ROW_H + 0.5,
       width = ((nrow(df) * TREE_W) + (length_text * TEXT_W + nrow(df) * TREE_W * MARGIN_ERROR_CORRECTION) + (ncol(df) * COL_W)) * 1.1 + 4/3)

print((((nrow(df) * TREE_W) + (length_text * TEXT_W + nrow(df) * TREE_W * MARGIN_ERROR_CORRECTION) + (ncol(df) * COL_W)) * 1.1 + 4/3)*300)



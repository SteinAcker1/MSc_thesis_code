if(!requireNamespace("devtools")) install.packages("devtools")

packages_list <- list(
  "Seurat" = "4.1.0",
  "tidyverse" = "1.3.1",
  "ggpubr" = "0.4.0",
  "viridis" = "0.6.2",
  "clusterProfiler" = "4.2.2",
  "org.Hs.eg.db" = "3.14.0"
)

for(i in names(packages_list)) {
  packages_list[[i]] <- package_version(packages_list[[i]])
  devtools::install_version(package = i, version = packages_list[[i]], repos = "http://cran.us.r-project.org")
}
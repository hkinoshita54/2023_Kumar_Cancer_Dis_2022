# install Seurat V5

# first install 'remotes' and 'devtools'
install.packages("remotes")
install.packages("devtools")

# upon first run of the next line, may be prompted to install 'command line tools' (outside of Rstudio)
remotes::install_github("satijalab/seurat", "seurat5", quiet = F)

# set github personal access token to avoid errors related to "API rate limit"
usethis::create_github_token()
gitcreds::gitcreds_set()
usethis::edit_r_environ() # GITHUB_PAT=XXXXXXX...(git hub personal access token)

remotes::install_github("satijalab/seurat-data", "seurat5", quiet = F)
remotes::install_github("mojaveazure/seurat-object", "seurat5", quiet = F)
BiocManager::install("Rsamtools") # Rsamtools should be installed before Azimuth
remotes::install_github("stuart-lab/signac", "seurat5", quiet = F) # Signac should be installed before Azimuth
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") # BSgenome.Hsapiens.UCSC.hg38 should be installed before Azimuth
remotes::install_github("satijalab/azimuth", "seurat5", quiet = F)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = F)

# before installing BPCells, you need to install HDF5. refer to the developers github
# https://github.com/bnprks/BPCells
remotes::install_github("bnprks/BPCells", quiet = F)

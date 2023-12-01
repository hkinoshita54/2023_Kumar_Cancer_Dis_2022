# install Seurat V5
# first installed (and loaded) 'remotes' and 'devtools'
# upon first run of the next line, prompted to install 'command line tools' (outside of Rstudio) -> installed it 
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

# set github personal access token to avoid errors related to "API rate limit"
usethis::create_github_token()
gitcreds::gitcreds_set()
usethis::edit_r_environ() # GITHUB_PAT=XXXXXXX...(git hub personal access token)

remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("mojaveazure/seurat-object", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE) # Signac should be installed before Azimuth
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)

remotes::install_github("bnprks/BPCells", quiet = TRUE)

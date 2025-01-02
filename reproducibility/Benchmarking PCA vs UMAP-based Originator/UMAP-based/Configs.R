hca_lung <- list(
  name = "lung",
  uri = "https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0",
  ann_col = "ann_level_4",
  query = c(
    "B cells",
    "NK cells",
    "Classical monocytes",
    "Non-classical monocytes",
    "Alveolar macrophages",
    "Interstitial macrophages",
    "DC1",
    "DC2",
    "Plasma cells",
    "CD4 T cells",
    "CD8 T cells"
  ),
  ref = c(
    "B cell",
    "nk cell",
    "monocyte",
    "monocyte",
    "macrophage",
    "macrophage",
    "DC",
    "DC",
    "plasma cell",
    "T-cell",
    "T-cell"
  ),
  unified = c(
    "B-cell",
    "NK-cell",
    "Monocyte",
    "Monocyte",
    "Macrophage",
    "Macrophage",
    "DC",
    "DC",
    "Plasma-cell",
    "T-cell",
    "T-cell"
  )
)

global_liver <- list(
  name = "global_liver",
  uri = "https://explore.data.humancellatlas.org/projects/04e4292c-f62f-4098-ae9b-fd69ae002a90",
  ann_col = "Predicted_labels_CellTypist",
  query = c(
    "B cells",
    "Cycling B cells",
    "Follicular B cells",
    "Germinal center B cells",
    "Memory B cells",
    "Naive B cells",
    "Transitional B cells",
    "CD16+ NK cells",
    "CD16- NK cells",
    "NK cells",
    "Transitional NK",
    "Cycling NK cells",
    "Classical monocytes",
    "Cycling monocytes",
    "Non-classical monocytes",
    "Monocytes",
    "Macrophages",
    "DC",
    "Transitional DC",
    "DC1",
    "DC2",
    "DC3",
    "Cycling DCs",
    "Migratory DCs",
    "Plasma cells",
    "Cycling T cells",
    "Cycling gamma-delta T cells",
    "Cytotoxic T cells",
    "Follicular helper T cells",
    "Helper T cells",
    "Memory CD4+ cytotoxic T cells",
    "Regulatory T cells",
    "T cells",
    "Type 1 helper T cells",
    "Type 17 helper T cells",
    "Tcm/Naive cytotoxic T cells",
    "Tcm/Naive helper T cells",
    "Tem/Effector cytotoxic T cells",
    "Tem/Effector helper T cells",
    "gamma-delta T cells"
  ),
  ref = c(
    "B cell",
    "B cell",
    "B cell",
    "B cell",
    "B cell",
    "B cell",
    "B cell",
    "nk cell",
    "nk cell",
    "nk cell",
    "nk cell",
    "nk cell",
    "monocyte",
    "monocyte",
    "monocyte",
    "monocyte",
    "macrophage",
    "DC",
    "DC",
    "DC",
    "DC",
    "DC",
    "DC",
    "DC",
    "plasma cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell"
  ),
  unified = c(
    "B-cell",
    "B-cell",
    "B-cell",
    "B-cell",
    "B-cell",
    "B-cell",
    "B-cell",
    "NK-cell",
    "NK-cell",
    "NK-cell",
    "NK-cell",
    "NK-cell",
    "Monocyte",
    "Monocyte",
    "Monocyte",
    "Monocyte",
    "Macrophage",
    "DC",
    "DC",
    "DC",
    "DC",
    "DC",
    "DC",
    "DC",
    "Plasma-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell",
    "T-cell"
  )
)

global_lung <-
  modifyList(global_liver, list(name = "global_lung"))

global_spleen <-
  modifyList(global_liver, list(name = "global_spleen"))

humanliver <- list(
  name = "humanliver",
  uri = "https://explore.data.humancellatlas.org/projects/4d6f6c96-2a83-43d8-8fe1-0f53bffd4674",
  ann_col = "cell_type",
  query = c(
    "mature B cell",
    "natural killer cell",
    "inflammatory macrophage",
    "non-inflammatory macrophage",
    "plasma cell",
    "alpha-beta T cell",
    "gamma-delta T cell"
  ),
  ref = c(
    "B cell",
    "nk cell",
    "macrophage",
    "macrophage",
    "plasma cell",
    "T-cell",
    "T-cell"
  ),
  unified = c(
    "B-cell",
    "NK-cell",
    "Macrophage",
    "Macrophage",
    "Plasma-cell",
    "T-cell",
    "T-cell"
  )
)

liver_psc_normal <- list(
  name = "psc_normal",
  uri = "https://cellxgene.cziscience.com/collections/0c8a364b-97b5-4cc8-a593-23c38c6f0ac5",
  ann_col = "cell_type",
  query = c(
    "CD4-positive, alpha-beta T cell",
    "macrophage",
    "monocyte",
    "plasma cell"
  ),
  ref = c(
    "T-cell",
    "macrophage",
    "monocyte",
    "plasma cell"
  ),
  unified = c(
    "T-cell",
    "Macrophage",
    "Monocyte",
    "Plasma-cell"
  )
)

liver_psc_biliary <-
  modifyList(liver_psc_normal, list(name = "psc_biliary"))

liver_psc_sclerosing <-
  modifyList(liver_psc_normal, list(name = "psc_sclerosing"))

lung_3ca_xing_normal <- list(
  name = "3ca_lung_xing_normal",
  uri = "https://www.science.org/doi/10.1126/sciadv.abd9738",
  ann_col = "cell_type",
  query = c(
    "B_cell",
    "Dendritic",
    "Macrophage",
    "Monocyte",
    "NK_cell",
    "T_cell"
  ),
  ref = c(
    "B cell",
    "DC",
    "macrophage",
    "monocyte",
    "nk cell",
    "T-cell"
  ),
  unified = c(
    "B-cell",
    "DC",
    "Macrophage",
    "Monocyte",
    "NK-cell",
    "T-cell"
  )
)

lung_3ca_xing_luad <-
  modifyList(lung_3ca_xing_normal, list(name = "3ca_lung_xing_luad"))

lung_3ca_xing_nodule <-
  modifyList(lung_3ca_xing_normal, list(name = "3ca_lung_xing_nodule"))

lung_3ca_bischoff_normal <- list(
  name = "3ca_lung_bischoff_normal",
  uri = "https://www.nature.com/articles/s41388-021-02054-3",
  ann_col = "cell_type",
  query = c(
    "B_cell",
    "Dendritic",
    "Macrophage",
    "Monocyte",
    "NK_cell",
    "Plasma",
    "T_cell"
  ),
  ref = c(
    "B cell",
    "DC",
    "macrophage",
    "monocyte",
    "nk cell",
    "plasma cell",
    "T-cell"
  ),
  unified = c(
    "B-cell",
    "DC",
    "Macrophage",
    "Monocyte",
    "NK-cell",
    "Plasma-cell",
    "T-cell"
  )
)

lung_3ca_bischoff_tumor <-
  modifyList(lung_3ca_bischoff_normal, list(name = "3ca_lung_bischoff_tumor"))

kidney_3ca_krishna_normal <- list(
  name = "3ca_kidney_krishna_normal",
  uri = "https://www.sciencedirect.com/science/article/pii/S1535610821001653",
  ann_col = "cell_type",
  query = c(
    "B_cell",
    "Dendritic",
    "Macrophage",
    "Monocyte",
    "NK_cell",
    "T_cell"
  ),
  ref = c(
    "B cell",
    "DC",
    "macrophage",
    "monocyte",
    "nk cell",
    "T-cell"
  ),
  unified = c(
    "B-cell",
    "DC",
    "Macrophage",
    "Monocyte",
    "NK-cell",
    "T-cell"
  )
)

kidney_3ca_krishna_tumor <-
  modifyList(kidney_3ca_krishna_normal, list(name = "3ca_kidney_krishna_tumor"))

liver_guilliam <- list(
  name = "liver_guilliam",
  uri = "https://cellxgene.cziscience.com/collections/74e10dc4-cbb2-4605-a189-8a1cd8e44d8c",
  ann_col = "cell_type",
  query = c(
    # "monocyte",
    # "macrophage",
    # "conventional dendritic cell",
    # "plasmacytoid dendritic cell",
    # "B cell",
    # "natural killer cell",
    # "plasma cell",
    "T cell"
  ),
  ref = c(
    # "monocyte",
    # "macrophage",
    # "DC",
    # "DC",
    # "B cell",
    # "nk cell",
    # "plasma cell",
    "T-cell"
  ),
  unified = c(
    # "Monocyte",
    # "Macrophage",
    # "DC",
    # "DC",
    # "B-cell",
    # "NK-cell",
    # "Plasma-cell",
    "T-cell"
  )
)

lung_szabo <- list(
  name = "lung_szabo",
  ann_col = "cell_type",
  query = c("T-cell"),
  ref = c("T-cell"),
  unified = c("T-cell")
)

#!/bin/bash
#SBATCH --job-name=originator_TBannotation_simulatedData_mono_64
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=32g
#SBATCH --mail-user=thatchau@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/logs/%x-%j.log
#SBATCH --account=lgarmire99
#SBATCH --partition=standard,largemem

dir='/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/'
cd $dir

# Replace options accordingly
# --name Name of the job
# --path Path to the rds file (query data)
# --reference Path to the rds file (reference data)
# --output Output path
# --query-celltype-col Column in meta data for cell type annotation (query data)
# --ref-celltype-col Column in meta data for cell type annotation (reference data)
# --query-celltypes Cell type names (query data), separated by '|', must match same length as --ref-celltypes and --unified-celltypes, if not enough pad with 'None'
# --ref-celltypes Cell type names (reference data), separated by '|', must match same length as --query-celltypes and --unified-celltypes, if not enough pad with 'None'
# --unified-celltypes Cell type names (saved data), separated by '|', must match same length as --query-celltypes and --ref-celltypes
# --seed Random state seed
/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.3.1/bin/Rscript src/Main.R \
  --name 'TBannotation_simulatedData_mono_64' \
  --path '/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/data_merged_tissue_blood_residents_DEMO_noNKcellTissueBlood_withGroundtruth.rds' \
  --reference '/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/Tabular_sapiens_blood_editedGeneNames_noduplicatedGenes_withLevel1Annotation.rds' \
  --output '/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/simulatedData_seeds/mono/64' \
  --query-celltype-col 'final_annotation' \
  --ref-celltype-col 'level1_annotation' \
  --query-celltypes 'Monocyte' \
  --ref-celltypes 'monocyte' \
  --unified-celltypes 'Monocyte' \
  --seed '64'

#!/bin/bash
#SBATCH --job-name=originator_kidney_krishna_TBannotation_cleanedccRCC_wholeblood_ref_PCAbased_10PC_b_0
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
/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.3.1/bin/Rscript /nfs/dcmb-lgarmire/thatchau/originator_GB_revision/src/Main.R \
  --name 'kidney_krishna_TBannotation_cleanedccRCC_wholeblood_ref_PCAbased_10PC_b_0' \
  --path '/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/krishna_2021_complete_response_with_paired_PBMC_removeBloodContaminationByOriginator_wholebloodRef_with_groudtruth.rds' \
  --reference '/nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/blood.rds' \
  --output '/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/TBannotation_cleanedccRCC_wholeblood_ref_PCAbased_10PC/b/0' \
  --query-celltype-col 'cell_type' \
  --ref-celltype-col 'level1_annotation' \
  --query-celltypes 'B_cell' \
  --ref-celltypes 'B cell' \
  --unified-celltypes 'B-cell' \
  --seed '0'

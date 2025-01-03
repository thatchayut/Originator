runPopscle <- function(
  vcf_path,
  barcodes_path,
  bam_path,
  popscle_path,
  output_path = "./result"
) {
  message(paste0("PATH=", Sys.getenv("PATH")))
  cmd <- paste(
    "./exec/popscle_pnas.sh",
    shQuote(vcf_path),
    shQuote(barcodes_path),
    shQuote(bam_path),
    shQuote(popscle_path),
    shQuote(output_path),
    "> ./logs/system.log 2> ./logs/error.log"
  )
  message(paste0("bash: ", cmd))
  system(cmd)
}

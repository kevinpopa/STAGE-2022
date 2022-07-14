args = commandArgs(trailingOnly = TRUE)

# Check that at least the input file is specified.
if (length(args)==0) {
  stop("The script input should be Rscript script_prepareHomer.R (input file : required) (output file name : optional)", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste(gsub('.{4}$', '', args[1]), "HOMER.txt", sep = ".")
}
df <- read.delim(args[1])

df_1 <- df[, c("RegionChrome", "RegionStart", "RegionEnd")]
df_1$peak_name <- paste(df_1$RegionChrome, df_1$RegionStart, df_1$RegionEnd, sep = "_")
df_1$NV <- "NV"
df_1$Strand <- "*"

write.table(df_1, args[2], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
pairs <- args[1] # .pairparse.txt or .pairdedup.txt to summarize
pname <- substr(pairs, 0, nchar(pairs) - 4)

tmp <- read.table(paste0(pairs), header = F, sep = "\t", row.names = 1)
tmp1 <- tmp %>% filter(row.names(tmp) %in% c("total", "total_unmapped", "total_single_sided_mapped", "total_mapped", "total_dups", "total_nodups", "cis_1kb+", "cis_20kb+", "cis", "trans")) %>% 
t() %>% as.data.frame() %>%
mutate(Sequenced_Read_Pairs = sum(total_unmapped, total_single_sided_mapped, total_dups, total_nodups),
      Mapped_Read_Pairs = sum(total_dups, total_nodups),
      Nonclonal_Read_Pairs = total_nodups,
      Intra_chromosomal = cis,
      short_range_1k = cis - `cis_1kb+`,
      long_range_1k = `cis_1kb+`,
      Inter_chromosomal = trans,
      Unmapped = sum(total_unmapped, total_single_sided_mapped),
      Removed_Duplicates = total_dups) %>% 
select(c(11:19)) %>% 
t() %>% as.data.frame() %>%
mutate(ratio = ifelse(row_number() %in% c(1,2,8), V2/V2[[1]], V2/V2[[3]])) %>%
mutate(ratio = ifelse(row_number() %in% c(3,9), V2/V2[[2]], ratio)) %>%
mutate(ratio = ifelse(row_number() %in% c(5,6), V2/V2[[4]], ratio))

write.table(tmp1, file = paste0(pname,".summary.txt"), row.names = T, col.names = F, sep = "\t", quote = F)

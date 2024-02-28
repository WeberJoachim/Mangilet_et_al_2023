library("tidyverse")


process_motif_data <- function(data, motif_column, filter_motifs) {
  motif = unlist(strsplit(motif_column, ".", fixed = T))[1]
  data %>%
    separate_rows(motif_column, sep = '\\),' ) %>%
    separate(., motif_column, into = c(paste0("Number_Motif_", motif), paste0("Strand_", motif), paste0("Value_", motif_column)), sep = ",") %>%
    separate(., paste0("Number_Motif_", motif), into = c(paste0("Distance_", motif), paste0("Motif_", motif)), sep = "\\(") %>%
    filter(., get(paste0("Strand_", motif)) == Strand & !motif %in% filter_motifs)
}







	args <- commandArgs(trailingOnly = TRUE)

	if (length(args) != 1) {
  		stop("Usage: Rscript your_script.R <input_file>")
	}

	input_file <- args[1]

	name = substr(input_file, start = 1, stop = nchar(input_file)-4)
	print(name)
        data = read.delim(input_file, header = T)
       

        data = process_motif_data(data, "AAUAAA.Distance.From.Peak.sequence.strand.conservation.", c("AAAAAA", "TTTTTT") )
        data = process_motif_data(data, "UUGUUU.Distance.From.Peak.sequence.strand.conservation.", c("TTTTTT", "AAAAAA") )
        data = process_motif_data(data, "TGTA.Distance.From.Peak.sequence.strand.conservation.", c(""))
        data = process_motif_data(data, "YA.Distance.From.Peak.sequence.strand.conservation.", c(""))
	
	
	colnames(data)

        data$Distance_AAUAAA 	= as.numeric(data$Distance_AAUAAA)
        data$Distance_UUGUUU 	= as.numeric(data$Distance_UUGUUU)
        data$Distance_TGTA 	= as.numeric(data$Distance_TGTA)
        data$Distance_YA 	= as.numeric(data$Distance_YA)

	print("4")

        positive = data %>% filter(Strand == "+" & 131 > (Distance_AAUAAA - Distance_TGTA) & (Distance_AAUAAA - Distance_TGTA) > 59 & 28 > (Distance_UUGUUU - Distance_AAUAAA) & (Distance_UUGUUU - Distance_AAUAAA) > 5 & Distance_YA - Distance_AAUAAA > 0 )
        negative = data %>% filter(Strand == "-" & 131 > (Distance_TGTA - Distance_AAUAAA) & (Distance_TGTA - Distance_AAUAAA) > 59 & 20 > (Distance_UUGUUU - Distance_AAUAAA) & (Distance_UUGUUU - Distance_AAUAAA) > 5 & Distance_AAUAAA - Distance_YA > 0 )


        perfekte = rbind(positive, negative)
        positionen_perfekt = perfekte %>% select(Chr, Start, End, Strand) %>% distinct()
	print("5")
        write.table(positionen_perfekt, paste0(name, "_perfect_motifs.bed"), quote = F, col.names = F, row.names = F, sep = "\t")


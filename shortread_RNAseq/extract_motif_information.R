	library("tidyverse")
	args <- commandArgs(trailingOnly = TRUE)

	if (length(args) != 1) {
  		stop("Usage: Rscript your_script.R <input_file>")
	}

	input_file <- args[1]

	name = substr(input_file, start = 1, stop = nchar(input_file)-4)
	print(name)
        data = read.delim(input_file, header = T)
        data = separate_rows(data, "AAUAAA.Distance.From.Peak.sequence.strand.conservation.", sep = '\\),' )
        data = separate(data, "AAUAAA.Distance.From.Peak.sequence.strand.conservation.", into = c("A_Number_Motif", "A_Strand", "A_Value"), sep = ",")
        data = separate(data, "A_Number_Motif", into = c("A_Distance", "A_Motif"), sep = "\\(")
        data = data %>% filter(Strand == A_Strand & A_Motif != "TTTTTT" & A_Motif != "AAAAAA")

	print("1")
	print(head(data))

        data = separate_rows(data, "UUGUUU.Distance.From.Peak.sequence.strand.conservation.", sep = '\\),' )
        data = separate(data, "UUGUUU.Distance.From.Peak.sequence.strand.conservation.", into = c("U_Number_Motif", "U_Strand", "U_Value"), sep = ",")
        data = separate(data, "U_Number_Motif", into = c("U_Distance", "U_Motif"), sep = "\\(")
        data = data %>% filter(Strand == U_Strand & U_Motif != "TTTTTT" & U_Motif != "AAAAAA")

	print("2")

        data = separate_rows(data, "TGTA.Distance.From.Peak.sequence.strand.conservation.", sep = '\\),' )
        data = separate(data, "TGTA.Distance.From.Peak.sequence.strand.conservation.", into = c("USE_Number_Motif", "USE_Strand", "USE_Value"), sep = ",")
        data = separate(data, "USE_Number_Motif", into = c("USE_Distance", "USE_Motif"), sep = "\\(")
        data = data %>% filter(Strand == USE_Strand)

	print("3")

        data = separate_rows(data, "YA.Distance.From.Peak.sequence.strand.conservation.", sep = '\\),' )
        data = separate(data, "YA.Distance.From.Peak.sequence.strand.conservation.", into = c("YA_Number_Motif", "YA_Strand", "YA_Value"), sep = ",")
        data = separate(data, "YA_Number_Motif", into = c("YA_Distance", "YA_Motif"), sep = "\\(")
        data = data %>% filter(Strand == YA_Strand)

        data$A_Distance = as.numeric(data$A_Distance)
        data$U_Distance = as.numeric(data$U_Distance)
        data$USE_Distance = as.numeric(data$USE_Distance)
        data$YA_Distance = as.numeric(data$YA_Distance)

	print("4")

        positive = data %>% filter(Strand == "+" & 131 > (A_Distance - USE_Distance) & (A_Distance - USE_Distance) > 59 & 28 > (U_Distance - A_Distance) & (U_Distance - A_Distance) > 5 & YA_Distance - A_Distance > 0 )
        negative = data %>% filter(Strand == "-" & 131 > (USE_Distance - A_Distance) & (USE_Distance - A_Distance) > 59 & 20 > (U_Distance - A_Distance) & (U_Distance - A_Distance) > 5 & A_Distance - YA_Distance > 0 )


        perfekte = rbind(positive, negative)
        positionen_perfekt = perfekte %>% select(Chr, Start, End, Strand) %>% distinct()
	print("5")
        write.table(positionen_perfekt, paste0(name, "_perfect_motifs.bed"), quote = F, col.names = F, row.names = F, sep = "\t")


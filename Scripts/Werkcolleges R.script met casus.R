#Dag 1

#Stel je working directory in.
setwd("C:/Users/hieke/OneDrive - NHL Stenden/Leerjaar 2/Periode 4/Transcriptomics/Werkcollege")
getwd()
# Vervang de bestandsnaam hieronder met je eigen zip-bestand
unzip("ethanol_data.zip", exdir = "ethanol_data") #Hiermee worden de bestanden uitgepakt in een submap 'ethanol_data'
install.packages('BiocManager')
BiocManager::install('Rsubread')
library(Rsubread)
#het referentiegenoom indexeren
buildindex(
  basename = 'ref_ecoli',
  reference = 'GCF_000005845.2_ASM584v2_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
# Ethanol monsters
align.eth1 <- align(index = "ref_ecoli", readfile1 = "SRR8394576_ethanol_12h_1.fasta.gz", output_file = "eth1.BAM")
align.eth2 <- align(index = "ref_ecoli", readfile1 = "SRR8394577_ethanol_12h_2.fasta.gz", output_file = "eth2.BAM")
align.eth3 <- align(index = "ref_ecoli", readfile1 = "SRR8394578_ethanol_12h_3.fasta.gz", output_file = "eth3.BAM")
align.ctrl1 <- align(index = "ref_ecoli", readfile1 = "SRR8394612_control_12h_1.fasta.gz", output_file = "ctrl1.BAM")
align.ctrl2 <- align(index = "ref_ecoli", readfile1 = "SRR8394613_control_12h_2.fasta.gz", output_file = "ctrl2.BAM")
align.ctrl3 <- align(index = "ref_ecoli", readfile1 = "SRR8394614_control_12h_3.fasta.gz", output_file = "ctrl3.BAM")
# Laad Rsamtools voor sorteren en indexeren
library(Rsamtools)

# Bestandsnamen van de monsters
samples <- c('eth1', 'eth2', 'eth3', 'ctrl1', 'ctrl2', 'ctrl3')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

#Dag 2

library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)
############
#Niet uitvoeren met casus!!!

# Inlezen en filteren van GFF3-bestand
gff <- read_tsv("Escherichia_coli_str_k_12_substr_mg1655.gff3.gz", comment = "#", col_names = FALSE)
# Kolomnamen toevoegen
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# Alleen genregels selecteren
gff_gene <- gff %>% filter(type == "gene")
# 'type' aanpassen naar 'exon' zodat featureCounts het accepteert
gff_gene$type <- "exon"
# Extraheer de chromosoomnaam uit BAM-header
bam_chr <- names(scanBamHeader("eth1.BAM")[[1]]$targets)[1]
gff_gene$seqid <- bam_chr
write_delim(gff_gene, "ecoli_ready_for_featureCounts.gtf", delim = "\t", col_names = FALSE)

#Niet uitvoeren met casus!!!
##########

#Deze stappen wel weer uitvoeren
# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- c("eth1.BAM", "eth2.BAM", "eth3.BAM", "ctrl1.BAM", "ctrl2.BAM", "ctrl3.BAM")

count_matrix <- featureCounts(
  files = allsamples ,
  annot.ext = "ecoli_ready_for_featureCounts.gtf",
  isPairedEnd = FALSE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

#Resultaten bekijken
head(count_matrix$annotation)
head(count_matrix$counts)

#Countmatrix opslaan en inladen
# Bekijk eerst de structuur van het object
str(count_matrix)

# Haal alleen de matrix met tellingen eruit
counts <- count_matrix$counts

#Stel duidelijke kolomnamen in, zodat ze overeenkomen met je samples. Dit helpt bij herkenning tijdens visualisaties en analyses.
colnames(counts) <- c("eth1", "eth2", "eth3", "ctrl1", "ctrl2", "ctrl3")

#De eerste kolom bevat de gen-ID’s. Je kunt deze beter als rij-namen gebruiken (in plaats van als gewone kolom).


#Nu de matrix klaar is, sla je deze op. Zo kun je hem in Werkcollege 3 direct opnieuw inlezen.
write.csv(counts, "bewerkt_countmatrix.csv")
head(counts)

#Dag 3

#Bestand inladen
counts <- read.csv("bewerkt_countmatrix.csv", row.names = 1)

#In deze stap maken we een tabel die beschrijft welke monsters behandeld zijn met ethanol en welke als controle dienden. Deze informatie is essentieel voor de differentiële expressieanalyse.
treatment <- c("ethanol", "ethanol", "ethanol", "control", "control", "control")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c('eth1', 'eth2', 'eth3', 'ctrl1', 'ctrl2', 'ctrl3')

install.packages('BiocManager')
BiocManager::install('DESeq2')
BiocManager::install('KEGGREST')

library(DESeq2)
library(KEGGREST)

# Maak DESeqDataSet aan
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.

write.table(resultaten, file = 'ResultatenWC3.csv', row.names = TRUE, col.names = TRUE)

#Stap 1: Hoeveel genen zijn er écht veranderd? Hier tellen we hoeveel genen er significant op- of neer-gereguleerd zijn.
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

#Stap 2: Welke genen springen eruit? Nu sorteren we het resultaat om te kijken naar de opvallendste genen:
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

#Bekijk nu welke genen het belangrijkst zijn volgens de analyse.
head(laagste_p_waarde)

#Een volcano plot laat zien welke genen significant veranderen in expressie. Op de x-as staat de log2 fold change en op de y-as de -log10 van de aangepaste p-waarde (padj). Alleen genen met padj < 0,05 en log2fc van > 1 of < -1 worden gelabeled.
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

#We kunnen ook plots maken met een andere opmaak.
# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

#Sla het figuur op
dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()


#Installeer en laad het pathview pakket
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

#Visualiseer een KEGG-pathway
#Met pathview() maken we een plot van een specifiek KEGG-pathway. In dit voorbeeld gebruiken we “eco02026”, wat biofilmvorming representeert.
resultaten[1] <- NULL
resultaten[2:5] <- NULL

pathview(
  gene.data = resultaten,
  pathway.id = "eco02026",  # KEGG ID voor Biofilm formation – E. coli
  species = "eco",          # 'eco' = E. coli in KEGG
  gene.idtype = "KEGG",     # Geef aan dat het KEGG-ID's zijn
  limit = list(gene = 5)    # Kleurbereik voor log2FC van -5 tot +5
) 

#Een ander voorbeeld: Werken met KEGG pathways in R
#Met de functie keggLink() uit het KEGGREST-pakket kun je in R koppelingen leggen tussen genen en pathways in de KEGG-database.

#Gen → Pathways
keggLink("pathway", "eco:b0221")
#✅ Handig om te achterhalen in welke netwerken jouw gen een rol speelt.

#Pathway → Genen
keggLink("eco", "path:eco00010")
#✅ Dit gebruik je om te zien welke andere genen mogelijk ook veranderen in dezelfde pathway.

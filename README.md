# Transcriptomics
Het toepassen van transcriptomics om genexpressiepatronen op te sporen die wijzen op vroege stadia van reumatoïde artritis.
# Inhoud

# Inleiding
Reumatoïde artritis (RA) is een chronische auto-immuunziekte waarbij gewrichten ontstoken raken. Dit veroorzaakt schade aan kraakbeen en bot. Typische klachten zijn ochtendstijfheid, gewrichtspijn, vermoeidheid en gewichtsverlies.

De oorzaak van RA is een combinatie van genetische aanleg en omgevingsfactoren, zoals roken, infecties en het darmmicrobioom. Bij RA raken afweercellen overactief en ontstaat er een aanhoudende ontstekingsreactie. Autoantilichamen kunnen daarbij extra schade aanrichten. Er is geen remedie voor RA. Om de klachten te verlichten en gewrichtsschade te voorkomen wordt vaak gebruik gemaakt van ontstekingsremmers. [(Bron)](https://github.com/vdmH46/Transcriptomics/blob/main/Overview%20of%20Rheumatoid%20Arthritis%20and%20Scientific%20Understanding%20of%20the%20Disease.pdf)

Om RA beter te begrijpen, vroegtijdig te herkennen en ernstige schade aan gewrichten en botten te voorkomen, wordt gebruikgemaakt van transcriptomics. Deze techniek onderzoekt welke genen verhoogd of verlaagd tot expressie komen bij mensen met RA en welke metabole routes anders functioneren.
# Methode
In dit onderzoek werd RNA-sequencing data geanalyseerd van acht synoviumbiopten: vier van gezonde individuen en vier van personen met reumatoïde artritis (RA). RA-patiënten hadden een bevestigde diagnose van meer dan 12 maanden en testten positief op ACPA. De gezonde individuen waren ACPA-negatief. De data is afkomstig uit een eerder gepubliceerd onderzoek (Platzer et al., 2019).

De ruwe sequencing data werd verwerkt in RStudio (versie 4.4.1). Eerst werd met het *Rsubread*-package (versie 2.20.0) een referentie-index opgebouwd op basis van het humane referentiegenoom GRCh38 (Homo_sapiens.GRCh38.dna.toplevel_1.fa). De reads werden uitgelijnd met de align()-functie van *Rsubread*. De verkregen BAM-bestanden werden gesorteerdmet behulp van het *Rsamtools*-package (versie 2.22.0).

De reads werden geteld met featureCounts() op basis van een GTF-annotatiebestand dat informatie bevat over de ligging van genen (Homo_sapiens.GRCh38.114.gtf). De countmatrix werd geanalyseerd met *DESeq2* (versie 1.46.0) om differentieel tot expressie komende genen tussen RA en controle te identificeren (p-waarde < 0,05; log2 fold change ≥ 1). De resultaten werden gevisualiseerd in een *volcano plot* met het *EnhancedVolcano*-package (1.24.0).

Vervolgens werd een Gene Ontology (GO)-verrijkingsanalyse uitgevoerd met *goseq* (1.58.0). Significante termen werden gevisualiseerd met *ggplot2* (3.5.2). Tot slot werd KEGG-pathway analyse uitgevoerd met het *pathview*-package (1.46.0), waarbij onder andere pathway hsa04670 (“Leukocyte transendothelial migration”) werd gevisualiseerd.
# Resultaten
Om inzicht te krijgen in welke genen significant verschillend tot expressie kwamen tussen de RA-groep en gezonde controles, werden 29.407 genen geanalyseerd. De resultaten zijn gevisualiseerd in een *volcano plot* (Figuur 1). Hierin staat de log₂ fold change (> 1) op de x-as, en de –log₁₀(p-waarde < 0,05) op de y-as.

Om beter te begrijpen welke biologische processen geassocieerd zijn met de differentieel tot expressie komende genen bij RA, werd een Gene Ontology (GO)-verrijkingsanalyse uitgevoerd met het *goseq*-package. Uit deze analyse kwamen meerdere verrijkte termen naar voren die sterk gerelateerd zijn aan het immuunsysteem (Figuur 2). 

Ook werd een pathway-analyse uitgevoerd met behulp van het *pathview*-package. Hierbij werd de KEGG-pathway "Leukocyte transendothelial migration" (hsa04670) gevisualiseerd (Figuur 3). In deze pathway werd een duidelijke opregulatie waargenomen van meerdere cell adhesion molecules (CAMs), waaronder VCAM1, ICAM1 en PECAM1. 

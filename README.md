# Transcriptomics
Het toepassen van transcriptomics om genexpressiepatronen op te sporen die wijzen op vroege stadia van reumatoïde artritis.
# Inhoud
* `Bronnen/`- Gebruikte bronnen
* `Data/`- Verwerkte datasets
* `Data_R_Raw/`- De verkregen ruwe data
* `Dataset/`- Artikel met dataset
* `Resultaten/`- Figuren uit de R analyse
* `Scripts/`- Script met de uitgevoerde stappen
* `Flowchart/`- Stroomschema van methode
* `README.md/`- Het document om de tekst hier te genereren
* `Data_stewardship/`- Voor de competentie beheren
# Inleiding
Reumatoïde artritis (RA) is een chronische auto-immuunziekte waarbij gewrichten ontstoken raken. Dit veroorzaakt schade aan kraakbeen en bot. Typische klachten zijn ochtendstijfheid, gewrichtspijn, vermoeidheid en gewichtsverlies.

De oorzaak van RA is een combinatie van genetische aanleg en omgevingsfactoren, zoals roken, infecties en het darmmicrobioom. Bij RA raken afweercellen overactief en ontstaat er een aanhoudende ontstekingsreactie. Autoantilichamen kunnen daarbij extra schade aanrichten. Er is geen remedie voor RA. Om de klachten te verlichten en gewrichtsschade te voorkomen wordt vaak gebruik gemaakt van ontstekingsremmers. [(Jahid et al., 2023)](https://github.com/vdmH46/Transcriptomics/blob/main/Overview%20of%20Rheumatoid%20Arthritis%20and%20Scientific%20Understanding%20of%20the%20Disease.pdf)

Om RA beter te begrijpen, vroegtijdig te herkennen en ernstige schade aan gewrichten en botten te voorkomen, wordt gebruikgemaakt van transcriptomics. Deze techniek onderzoekt welke genen verhoogd of verlaagd tot expressie komen bij mensen met RA en welke metabole routes anders functioneren. [(Alivernini et al., 2022)](https://github.com/vdmH46/Transcriptomics/blob/main/Bronnen/The%20pathogenesis%20of%20rheumatoid%20arthritis.pdf)
# Methode
In dit onderzoek werd RNA-sequencing data geanalyseerd van acht synoviumbiopten: vier van gezonde individuen en vier van personen met reumatoïde artritis (RA). RA-patiënten hadden een bevestigde diagnose van meer dan 12 maanden en testten positief op ACPA. De gezonde individuen waren ACPA-negatief. De data is afkomstig uit een eerder gepubliceerd onderzoek [(Platzer et al., 2019)](https://github.com/vdmH46/Transcriptomics/blob/main/Dataset/Artikel%20Casus.pdf).

De ruwe sequencing data werd verwerkt in RStudio (versie 4.4.1). Eerst werd met het *Rsubread*-package (versie 2.20.0) een referentie-index opgebouwd op basis van het humane referentiegenoom GRCh38 (Homo_sapiens.GRCh38.dna.toplevel_1.fa). De reads werden uitgelijnd met de align()-functie van *Rsubread*. De verkregen BAM-bestanden werden gesorteerd met behulp van het *Rsamtools*-package (versie 2.22.0).

De reads werden geteld met featureCounts() op basis van een GTF-annotatiebestand dat informatie bevat over de ligging van genen (Homo_sapiens.GRCh38.114.gtf). De countmatrix werd geanalyseerd met *DESeq2* (versie 1.46.0) om differentieel tot expressie komende genen tussen RA en controle te identificeren (p-waarde < 0,05; log2 fold change ≥ 1). De resultaten werden gevisualiseerd in een *volcano plot* met het *EnhancedVolcano*-package (1.24.0).

Vervolgens werd een Gene Ontology (GO)-verrijkingsanalyse uitgevoerd met *goseq* (1.58.0). Significante termen werden gevisualiseerd met *ggplot2* (3.5.2). Tot slot werd KEGG-pathway analyse uitgevoerd met het *pathview*-package (1.46.0), waarbij onder andere pathway hsa04670 (“Leukocyte transendothelial migration”) werd gevisualiseerd.

[Flowchart](https://github.com/vdmH46/Transcriptomics/blob/main/Flowchart/Flowchart.png)
# Resultaten
Om inzicht te krijgen in welke genen significant verschillend tot expressie kwamen tussen de RA-groep en gezonde controles, werden 29.407 genen geanalyseerd. De resultaten zijn gevisualiseerd in een *volcano plot* ![(Figuur 1)](https://github.com/vdmH46/Transcriptomics/blob/main/Resultaten/VolcanoplotCasus.png) <sub><sup>Figuur: Top 10 verrijkte GO-termen op basis van differentieel tot expressie komende genen bij RA. De stipgrootte geeft het aantal genen aan, kleur geeft de -log₁₀(p-waarde) weer.</sup></sub>. Hierin staat de log₂ fold change (> 1) op de x-as, en de –log₁₀(p-waarde < 0,05) op de y-as.

Om beter te begrijpen welke biologische processen geassocieerd zijn met de differentieel tot expressie komende genen bij RA, werd een Gene Ontology (GO)-verrijkingsanalyse uitgevoerd met het *goseq*-package. Uit deze analyse kwamen meerdere verrijkte termen naar voren die sterk gerelateerd zijn aan het immuunsysteem [(Figuur 2)](https://github.com/vdmH46/Transcriptomics/blob/main/Resultaten/GO_Analyse.png). 

Ook werd een pathway-analyse uitgevoerd met behulp van het *pathview*-package. Hierbij werd de KEGG-pathway "Leukocyte transendothelial migration" (hsa04670) gevisualiseerd [(Figuur 3)](https://github.com/vdmH46/Transcriptomics/blob/main/Resultaten/hsa04670.pathview.png). In deze pathway werd een duidelijke opregulatie waargenomen van meerdere cell adhesion molecules (CAMs), waaronder VCAM1, ICAM1 en PECAM1. 
# Conclusie
Dit onderzoek maakte gebruik van transcriptomics om inzicht te krijgen in genexpressieveranderingen bij reumatoïde artritis (RA). Uit RNA-sequencing data van synoviumbiopten van RA-patiënten en gezonde controles bleek dat 1.690 genen significant verschillend tot expressie kwamen, met zowel verhoogde als verlaagde genactiviteit. De GO-analyse liet zien dat vooral immuun gerelateerde processen betrokken waren, waaronder de activatie van immuuncellen en migratie van leukocyten. De KEGG-pathwayanalyse toonde aan dat de route voor 'Leukocyte transendothelial migration' betrokken was. Hierbij spelen cell adhesion molecules (CAMs) een sleutelrol. Deze moleculen, gelegen op het endotheel, ondersteunden de hechting van leukocyten en maakten hun doorgang naar ontstoken weefsel mogelijk. De verhoogde expressie van CAMs bij RA wees op een actiever migratieproces, wat paste bij het chronische ontstekingsproces waarbij immuuncellen zich ophopen in het gewricht. [(Veale & Maple, 1996)](https://github.com/vdmH46/Transcriptomics/blob/main/Bronnen/Cell%20adhesion%20molecules%20in%20rheumatoid%20arthritis.pdf).

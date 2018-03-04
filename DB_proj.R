setwd('/Users/mitya/Desktop/RAN/')
library("org.Hs.eg.db")
library(plyr)
library(ggplot2)
library(ggrepel)
library(seqinr)
library(ape) 
library(dplyr)
library (rentrez)


db <- read_delim("~/Desktop/RAN/Human Full Genome (1).txt",  +     "\t", escape_double = FALSE, trim_ws = TRUE)
db <- Human_Full_Genome_1_
c <- db[db$GeneAttribute !='--',]
c1 <- c[c$GCount != 0, ]
c2 <- c1[c1$PatternSize >= 3,]
c3 <- c2[c2$PatternSize <= 6,]
#c3 <- c2[c2$CopyNumber >= 10,]
c4 <- c3[grep('G', c3$Pattern),]


c5 <- c4[c4$MatchPerc >= 90,] #norm
write.table(c5, file = "norm_database.txt", sep = "\t")


c6 <- c4[c4$MatchPerc == 100,] #hard
c6$res <- c6$CopyNumber %% 1 == 0 #hard
c7 <- c6[c6$res  == TRUE,] #hard
write.table(c7, file = "hard_database.txt", sep = "\t")


#Cider_____________________________________________

cider <- read.csv('cider_db.csv')
cider1 <- cider[cider$Organism == 'Homo sapiens',]
cider2 <- cider1[cider1$Subject.Resource == 'Entrez',]
cid1 <- cider2[, c(5,22)]
cid <- cid1[!duplicated(cid1),]
cid$Subject <- as.character(cid$Subject)
cid$Disease <- as.character(cid$Disease)

#ORF_+_G-quad_+_Kozak_____________________________________________


norm <- read.csv('norm_with_v3.csv')
scores <- read.csv('with_scores.csv', sep = ',')

read.table()





str(norm)
norm$GeneID <- as.character(norm$GeneID)
norm$Genome <- as.character(norm$Genome)
norm$Pattern <- as.character(norm$Pattern)

norm$Symbol <- mapIds(org.Hs.eg.db,
                            keys= norm$GeneID,
                            column="SYMBOL",
                            keytype="ACCNUM",
                            multiVals="first")


norm <- read.csv('annotated_norm_db.csv', sep = ';')

#Names_extraction_____________________________________________


name_norm <- norm[!duplicated(norm$Symbol),]
name_norm$Symbol

#_____________________________________________

head(cid)
str(fine_short)

cidx <- cid[!duplicated(cid$Subject),]
cidx1 <- cidx
names(cidx1) <- c('Symbol','Disease')
normx <- norm

finalll <- merge(x = normx, y = cidx1, by="Symbol", all.x = TRUE)
write.table(finalll, file = "annotated_norm_db.txt", sep = "\t")
write.csv(finalll, file = "annotated_norm_db.csv", sep = ";")





cidx2 <- ddply(cidx1, .(Symbol), summarize, Disease = sum(Disease))

paste(df$letters, df$numbers, sep=""))



#_____________________________________________

norm <- read.csv('annotated_norm_db.csv', sep = ';')
norm$GeneID <- as.character(norm$GeneID)
names <- read.table('names.txt')
norm$IDD <- names$V1
norm$Orf <- as.character(norm$Orf)
norm$IDD <- as.character(norm$IDD)
norm$Disease <- as.character(norm$Disease)
norm$Symbol <- as.character(norm$Symbol)
norm$Pattern <- as.character(norm$Pattern)
norm$Chr <- as.character(norm$Chr)
#length(which(duplicated(norm$Orf) == T))
#normx <- norm[!duplicated(norm),]
str(norm)




scores <- read.csv('with_scores.csv', sep = ',')
scores$Orf <- as.character(scores$Genome)
scores$Genome <- NULL

# sco <- scores
# sco[,1] <- NULL 
# sco[,1] <- NULL 
# sco[,1] <- NULL 
# sco[,2] <- NULL
# sco[,2] <- NULL
# sco[,4] <- NULL
# sco1 <- sco
# sco <- sco1
# sco[,2] <- NULL
# sco[,3] <- NULL
# #length(which(duplicated(sco$Orf) == T))
# scox <- sco[!duplicated(sco),]
# str(sco)
#normscox <- normsco[!duplicated(normsco),]
#normsco <- merge(x = norm, sco, by=c("Orf","Copy"))
#normsco <- merge(norm, sco, by=c("Orf","Copy"))
#normsco <- merge(normx, scox, by="Orf")
#normsco <- merge(x = normx, y = scox, by="Orf", all.x = TRUE)


normsco <- merge(norm, scores)
normsco$X <- NULL
normsco$Is_causac <- NULL
normsco$Is_gquad <- NULL
str(normsco)




nam <- read.table('names_filtred_short_last_VERYLAST_edited.txt')
nam$IDD <- as.character(nam$V1)
nam$V1 <- NULL
str(nam)


normscox <- merge(nam, normsco)
str(normscox)


boxplot(normscox$Score_Left)
boxplot(normscox$Score_Right)

median(normscox$Score_Left)
median(normscox$Score_Right)

mean(normscox$Score_Left)
mean(normscox$Score_Right)


normsco1 <- normscox[normscox$Score_Left >= median(normscox$Score_Left),]
normsco2 <- normsco1[normsco1$Score_Right >= median(normscox$Score_Right),]


normsco1 <- normscox[normscox$Score_Left >= 2,]
normsco2 <- normsco1[normsco1$Score_Right >= 3,]

sum(normsco2$Kozak)*100/nrow(normsco2)
sum(normsco2$Gquadruplex)*100/nrow(normsco2)

write.csv(normsco2, file = "final_db.csv", sep = ";")




boxplot(normsco2$Score_Left)
boxplot(normsco2$Score_Right)


normscoxL <- normsco2[normsco2$Score_Left >= 5.5,]
normscoxR <- normsco2[normsco2$Score_Right >= 6,]

normscoxLR <- rbind(normscoxL, normscoxR)
#normscoxLR$IDD <- NULL
normscoxLR <- normscoxLR[!duplicated(normscoxLR),]

normscoxLRD <- normscoxLR[!is.na(normscoxLR$Disease),]



normsco3 <- normsco2[normsco2$Score_Left >= 4,]
normsco4 <- normsco3[normsco3$Score_Left >= 6,]

dim(normsco4)
dim(normscoxLR)

finn <- rbind(normsco4, normscoxLR)
finnx <- finn[!duplicated(finn),]

finnxD <- finnx[!is.na(finnx$Disease),]

boxplot(finnxD$Score_Left)
boxplot(finnxD$Score_Right)

sum(finnxD$Kozak)*100/nrow(finnxD)
sum(finnxD$Gquadruplex)*100/nrow(finnxD)


finxd <- finnxD[finnxD$Score_Left >= 4,]
finxd <- finxd[finxd$Score_Right >= 4,]

write.csv(finxd, file = "Disease_db_new.csv", sep = ";")
_______________






normscoK <- normsco[normsco$Kozak == 1,]

boxplot(normscoK$Score_Left)
boxplot(normscoK$Score_Right)

median(normscoK$Score_Left)
median(normscoK$Score_Right)

mean(normscoK$Score_Left)
mean(normscoK$Score_Right)


normscoK1 <- normscoK[normscoK$Score_Left >= 3.79,]
normscoK2 <- normscoK1[normscoK1$Score_Right >= 4.47,]
sum(normscoK2$Gquadruplex)*100/nrow(normscoK2)


write.csv(normscoK2, file = "Kozak_db.csv", sep = ";")





normscoKD <- normscoK2[!is.na(normscoK2$Disease),]

boxplot(normscoKD$Score_Left)
boxplot(normscoKD$Score_Right)
median(normscoKD$Score_Left)
median(normscoKD$Score_Right)

mean(normscoKD$Score_Left)
mean(normscoKD$Score_Right)


normscoKD1 <- normscoKD[normscoKD$Score_Left >= 4.97,]
normscoKD2 <- normscoKD1[normscoKD1$Score_Right >= 5.3,]
normscoKD2$IDD <- NULL


diseas <- rbind(normscoxLRD,normscoKD2)
diseas <- diseas[!duplicated(diseas),]
diseasK <- diseas[diseas$Kozak == 1,]

write.csv(diseasK, file = "Disease_db.csv", sep = ";")
#___

nam <- read.table('names_filteres_edited.txt')
nam$IDD <- as.character(nam$V1)
nam$V1 <- NULL
str(nam)


normscox <- merge(nam, normsco)








nam$GeneID <- as.character(nam$V1)
str(nam)
nam$V1 <- NULL
norm1 <- norm[grep(nav, norm$GeneID),]



norm1 <- merge(x = norm, y = nam, by="GeneID", all.y = TRUE)










#_____________________________________________ 
y tc

reg_left <- 'G[C,A][G,A][G,A]CGGCG[A,G][T,C].G..G.T...G'

reg_right <- '[T,C][A,G][A,G][A,T,G][A,G]G..G[T,G,C].......G........G'

grep(reg_left, norm1$Orf)
grep(reg_right, norm1$Orf)


norm2 <- norm1[grep(reg_right, norm1$Orf),]
norm3 <- norm2[grep(reg_left, norm2$Orf),]



#_____________________________________________



norm <- read.csv('norm_with_v3.csv')
scores <- read.csv('with_scores.csv', sep = ',')

read.table()

















#_____________________________________________
#parce AA seqs
seq <- vector()
for (i in 1:length(c)){
  query <- paste(c[i],'[Gene] AND Homo sapiens[Orgn] AND refseq[filter]', sep='')
  prot <- entrez_search(db="protein", query)
  if (prot$count > 0){
    aa_seq <- entrez_fetch(db="protein", id=prot$ids[1], rettype="fasta")
    seq[i] <- as.vector(aa_seq)
  }
}



entrez_db_summary("sra")
entrez_db_summary("clinvar")

entrez_db_searchable("cdd")
# PLEN 	 Length of the PSSM or domain search model
#SD 	 The desription of functional sites in a domain 
#NS 	 The number of functional sites in a domain 

entrez_db_searchable("clinvar")
#DIS 	 Diseases or traits associated with this record 
#GID 	 Gene ID 




library (rentrez)
f_gene <- entrez_search(db="gene", term="NR_027946")
idg <- f_gene$ids #gene id
summ <- entrez_summary(db="gene", id=idg)
summ$name
#summ$summary


all_the_links <- entrez_link(dbfrom='gene', id=idg, db='all')
all_the_links$links
clin <- all_the_links$links$gene_clinvar
clin_summ <- entrez_summary(db="clinvar", id=clin)
clin_summ$`443532`$clinical_significance[1]
clin_summ$clinical_significance
str(clin_summ)
clin_summ1 <- entrez_summary(db="clinvar", id=443532)

cdd <- all_the_links$links$gene_cdd




hard <- read.csv('hard_with_v1.csv')
str(hard)
hard$GeneID <- as.character(hard$GeneID)
hard$Genome <- as.character(hard$Genome)
hard$Pattern <- as.character(hard$Pattern)

hard$Symbol <- mapIds(org.Hs.eg.db,
                      keys= hard$GeneID,
                      column="SYMBOL",
                      keytype="ACCNUM",
                      multiVals="first")

name_hard <- hard[!duplicated(hard$Symbol),]
name_hard$Symbol

head(cid)
str(fine_short)

cidx <- cid[!duplicated(cid$Subject),]
cidx1 <- cidx
names(cidx1) <- c('Symbol','Disease')
hardx <- hard

total1 <- merge(x = hardx, y = cidx1, by="Symbol", all.x = TRUE)







cidx2 <- ddply(cidx1, .(Symbol), summarize, Disease = sum(Disease))

paste(df$letters, df$numbers, sep=""))




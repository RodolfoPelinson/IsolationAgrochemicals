library(vegan)

data <- read.csv("Data/com_orig.csv", header = T, stringsAsFactors = F, row.names = 1)
Trait <- read.csv("Data/Trait.csv", header = T, stringsAsFactors = F)
coord <- read.csv("Data/coord.csv", header = T, stringsAsFactors = F, row.names = 1)
#chlorophyll <- read.csv("Data/chlorophyll.csv", header = T, row.names = 1, sep = ",")



SS <- as.factor(data$survey)

ID <- as.factor(data$ID)
ID_SS1 <- ID[SS == "1"]
ID_SS2_3_4 <- ID[SS == "2"]


isolation_all <- factor(data$isolation, levels = c("30","120","480"))
isolation_SS1 <- isolation_all[SS == "1"]
isolation_SS2_3_4 <- isolation_all[SS == "2"]

treatments_all <- factor(data$treatments, levels = c("control","pasture","sugar_cane"))
treatments_SS1 <- treatments_all[SS == "1"]
treatments_SS2_3_4 <- treatments_all[SS == "2"]

com_orig <- data[,-c(1:6)]



#Different combinations
SS_1_234 <- as.factor(c(rep("1", 45),rep("234", 3*180)))
SS_1_2_34 <- as.factor(c(rep("1", 45),rep("2", 180),rep("34", 2*180)))
SS_12_34 <- as.factor(c(rep("12", 45+180),rep("34", 2*180)))
SS_12_3_4 <- as.factor(c(rep("12", 45+180),rep("3", 180),rep("4", 180)))
SS_123_4 <- as.factor(c(rep("123", 45+180+180),rep("4", 180)))
SS_14_23 <- as.factor(c(rep("14", 45),rep("23", 2*180),rep("14", 180)))
SS_1_23_4 <- as.factor(c(rep("1", 45),rep("23", 2*180),rep("4", 180)))

isolation30120_480 <- factor(as.character(c(rep("30_120", 30),rep("480", 15))), levels = c("30_120","480"))
isolation30_120480 <- factor(as.character(c(rep("30", 15),rep("120_480", 30))), levels = c("30","120_480"))

treatments_contpast_sug <- as.factor(c(rep(c("control_pasture","control_pasture","sugar_cane"), 15)))
treatments_cont_pastsug <- as.factor(c(rep(c("control","pasture_sugar_cane","pasture_sugar_cane"), 15)))
treatments_contsug_past <- as.factor(c(rep(c("control_sugar_cane","pasture","control_sugar_cane"), 15)))

isolation30120_480_SS2_3_4 <- rep(isolation30120_480,4)
isolation30_120480_SS2_3_4 <- rep(isolation30_120480,4)

treatments_SS2_3_4_contpast_sug <- rep(treatments_contpast_sug,4)
treatments_SS2_3_4_cont_pastsug <- rep(treatments_cont_pastsug,4)
treatments_SS2_3_4_contsug_past <- rep(treatments_contsug_past,4)
#############################

com_SS1 <- com_orig[SS == "1",]
com_SS2 <- com_orig[SS == "2",]
com_SS3 <- com_orig[SS == "3",]
com_SS4 <- com_orig[SS == "4",]

Trait <- Trait[-(which(Trait$genus == "Staphylinidae")),]
Trait <- Trait[,c(1,2,3,5,6,7)]
Trait$trait <- as.factor(Trait$trait)
Trait <- Trait[order(Trait$genus),]
Trait_orig <- Trait

com_SS1_orig <- com_SS1
com_SS2_orig <- com_SS2
com_SS3_orig <- com_SS3
com_SS4_orig <- com_SS4


com_oc <- decostand(com_orig, method = "pa")
Trait_orig$total_ab <- colSums(com_orig)




com_SS1_oc <- decostand(com_SS1, method = "pa")
com_SS1 <- com_SS1[,colSums(com_SS1_oc) >= 4]
Trait_SS1 <- Trait[which(colSums(com_SS1_oc)>= 4),]
Trait_SS1$total_ab <- colSums(com_SS1)

com_SS2_oc <- decostand(com_SS2, method = "pa")
com_SS2 <- com_SS2[,colSums(com_SS2_oc) >= 4]
Trait_SS2 <- Trait[which(colSums(com_SS2_oc)>= 4),]
Trait_SS2$total_ab <- colSums(com_SS2)

com_SS3_oc <- decostand(com_SS3, method = "pa")
com_SS3 <- com_SS3[,colSums(com_SS3_oc) >= 4]
Trait_SS3 <- Trait[which(colSums(com_SS3_oc)>= 4),]
Trait_SS3$total_ab <- colSums(com_SS3)

com_SS4_oc <- decostand(com_SS4, method = "pa")
com_SS4 <- com_SS4[,colSums(com_SS4_oc) >= 4]
Trait_SS4 <- Trait[which(colSums(com_SS4_oc)>= 4),]
Trait_SS4$total_ab <- colSums(com_SS4)

com <- com_orig[,colSums(com_SS1_oc) >= 4 | colSums(com_SS2_oc) >= 4 | colSums(com_SS3_oc) >= 4 | colSums(com_SS4_oc) >= 4]
Trait <- Trait_orig[which(colSums(com_SS1_oc) >= 4 | colSums(com_SS2_oc) >= 4 | colSums(com_SS3_oc) >= 4 | colSums(com_SS4_oc) >= 4),]



##################################################################################
##################################################################################
##################################################################################
#Summed communities

summarize_community <- function(com, ID){
  ponds <- list()
  for (i in 1:length(levels(ID))){
    ponds[[i]] <- com[ID == ID[i],]
  }
  ponds <- lapply(ponds, colSums)
  summarized_com <- as.data.frame(do.call(rbind, ponds))
  rownames(summarized_com) <- as.character(unique(ID_SS2_3_4))
  return(summarized_com)
}

sum_com_SS1 <- summarize_community(com_SS1, ID_SS1)
sum_com_SS2 <- summarize_community(com_SS2, ID_SS2_3_4)
sum_com_SS3 <- summarize_community(com_SS3, ID_SS2_3_4)
sum_com_SS4 <- summarize_community(com_SS4, ID_SS2_3_4)

Trait_SS1_sum <- Trait_SS1[order(Trait_SS1$total_ab, decreasing = TRUE),]
sum_com_SS1 <- sum_com_SS1[,order(Trait_SS1$total_ab, decreasing = TRUE)]

Trait_SS2_sum <- Trait_SS2[order(Trait_SS2$total_ab, decreasing = TRUE),]
sum_com_SS2 <- sum_com_SS2[,order(Trait_SS2$total_ab, decreasing = TRUE)]

Trait_SS3_sum <- Trait_SS3[order(Trait_SS3$total_ab, decreasing = TRUE),]
sum_com_SS3 <- sum_com_SS3[,order(Trait_SS3$total_ab, decreasing = TRUE)]

Trait_SS4_sum <- Trait_SS4[order(Trait_SS4$total_ab, decreasing = TRUE),]
sum_com_SS4 <- sum_com_SS4[,order(Trait_SS4$total_ab, decreasing = TRUE)]







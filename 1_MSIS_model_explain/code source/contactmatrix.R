## read in matrix of effective contacts relevant to airborne infectious diseases
USAcontactnew <- read.csv("https://raw.githubusercontent.com/mobs-lab/mixing-patterns/main/data/contact_matrices/United_States_country_level_M_overall_contact_matrix_85.csv",header = F)

# visualize contact matrix
library(RColorBrewer)
n.cols=100
nice.cols <-  colorRampPalette(brewer.pal(9, "YlOrRd"))(n.cols)
heatmap(as.matrix(USAcontactnew)/sum(diag(as.matrix(USAcontactnew))), Rowv=NA, Colv=NA, scale='none', col=nice.cols)

## create a symmetric contact matrix
symmetricUSA <- 0.5 * (USAcontactnew + t(USAcontactnew))

## expand the contact under 1 year old
expendUSA_1 <- matrix(data = symmetricUSA[1,1],nrow = 5,ncol=5)
expendUSA_2 <- matrix(data = as.numeric(rep(symmetricUSA[1,],5)),byrow = T,nrow = 5,ncol=85)
expendUSA <- cbind(expendUSA_1,expendUSA_2)
expendUSA <- rbind(expendUSA,cbind(t(expendUSA_2),as.matrix(symmetricUSA)))


colnames(expendUSA)=c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y",rep("2-4Y",3),rep("5-9Y",5),
                             rep("10-19Y",10),rep("20-39Y",20),
                             rep("40-59Y",20),rep("60Y+",25))

# aggregate ages into 13 age groups; average the number of contacts within each age group
ave_mat_x = aggregate( expendUSA, by = list(colnames(expendUSA)), FUN='sum' )
ave_mat = aggregate( t(ave_mat_x[,-1]), by = list(colnames(expendUSA)), FUN='mean')
rownames(ave_mat) <- ave_mat[,1]
colnames(ave_mat)[2:length(ave_mat)] <- rownames(ave_mat)
library(dplyr)
ageorder <- c("Group.1","<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")
contactUSA <- select(ave_mat,ageorder)
contactUSA <- contactUSA %>% mutate(Group.1=  factor(Group.1, levels = c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+"))) %>%
  arrange(Group.1)
contactUSA <- as.matrix(contactUSA[,-1])
## create a symmetric contact matrix
contactUSA[lower.tri(contactUSA, diag = FALSE)] <- 0
contactUSA <- contactUSA + t(contactUSA)
diag(contactUSA) <- 0.5*diag(contactUSA)
rownames(contactUSA) <-  c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")

# visualize contact matrix
heatmapUSA <- as.matrix(contactUSA)
heatmap(heatmapUSA/sum(diag(heatmapUSA)), Rowv=NA, Colv=NA, scale='none', col=nice.cols)

saveRDS(heatmapUSA,"~/Box/aim3/RSVtransmissionmodel/data_and_parms/contactmatrix.rds")

setwd("/Volumes/mum2023/z6/Documents/workingdir/WES/e112")
library(dplyr)
library(tidyr)
#############################################read 1327 OA variants with matching ENSG
tab <- read.table("ensgOA_org1327.txt", header = TRUE)
org1327 <- tab %>% distinct() #remove duplicate lines


#############################################Variant Consequences
#get variants with only one result-unique, or more than one conseq
unqvep <- org1327 %>% group_by(Uploaded_variation) %>%
  filter(n()==1)
unqvep_table <- data.frame(table(number=unqvep$Consequence))
#check variants with more than one conseq
dupvep <- org1327 %>% group_by(Uploaded_variation) %>% 
  filter(n()>1)
all_vep <- rbind(unqvep, dupvep)
#write.table(all_vep, "conseq-ensgOA_org1327.txt", quote = FALSE, row.names = FALSE, sep = "\t")


#############################################Compile masterfile: vep, bcftools, hgdp
orgvep <- filter(org1327, CANONICAL =="YES")
orgvep_form <- orgvep %>% separate_wider_delim(Location, delim = "-", names = c("Location", "LocationSplit"))
orgvep_form$chrompos <- paste0("chr",orgvep_form$Location)

bcfsnp <- read.table("bcfann-ORG-annpass-OA_snp.txt", header = TRUE, stringsAsFactors=FALSE)
bcfindel <- read.table("bcfann-ORG-annpass-OA_indel.txt", header = TRUE, stringsAsFactors=FALSE)
all_bcf <- rbind(bcfsnp, bcfindel)
sapply(all_bcf, class)
i <- c(7,10, 13)
all_bcf[ , i] <- apply(all_bcf[ , i], 2, function(x) as.numeric(as.character(x)))
all_bcf$chrompos <- paste0(all_bcf$CHROM,":",all_bcf$POS)
ORGVEP_BCF <- merge(x=orgvep_form, y=all_bcf, by.x ="chrompos", by.y = "chrompos", all.y =  TRUE) %>%
  dplyr::rename("AF_1KGP"=AF.x,"AA_VEP.e112"=AA,"AF_WES"=AF.y)

#get HGDP-ORG
e112_hgdpauto <- read.table("bcfannfreq-hgdporg-auto",header=TRUE)
e112_hgdpchrx<- read.table("bcfannfreq-hgdporg_chrX",header=TRUE)
e112_hgdp <- rbind(e112_hgdpauto, e112_hgdpchrx)
e112_hgdp$chrompos <- paste0(e112_hgdp$CHROM,":",e112_hgdp$POS)
e112_hgdp <- select(e112_hgdp, c(3:7,11,14,16)) %>% rename("REF_HGDP"=REF,"ALT_HGDP"=ALT,"AN_HGDP"=AN, "ID_HGDP"=ID)
i <- c(6, 7)
e112_hgdp[ , i] <- apply(e112_hgdp[ , i], 2, function(x) as.numeric(as.character(x)))
sapply(e112_hgdp, class)
#Export MasterFile: ORGVEP_BCF_HGDP (1,327 unique ORG)
ORGVEP_BCF_HGDP <- merge(x=ORGVEP_BCF, y=e112_hgdp, by.x ="chrompos", by.y = "chrompos", all.x =  TRUE)
#write.table(ORGVEP_BCF_HGDP, "ORGVEP_BCF_HGDP.txt", quote = FALSE, row.names = FALSE, sep = "\t")


#############################################CHECKING 
ancestralstate <- select(ORGVEP_BCF_HGDP, c(1,2,24,43,46,47,62)) %>% distinct()
write.table(hgdpoa, "ancestralstate2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
check <- filter(hgdpoa, hgdpoa$ALT != hgdpoa$ALT_HGDP)
hgdpoa <- select(ORGVEP_BCF_HGDP, c(1,2,24,43,46:51,59:65)) %>% distinct()
missing <- missing %>% distinct() #remove duplicate lines
noDAF <- filter(ancestralSNPVEP_BCF_HGDP, is.na(ancestralSNPVEP_BCF_HGDP$DAF_WES))
noDAF <-select(noDAF, c(1,17,22,23,36,42))%>% distinct() 


#############################################CALCULATION DAF,FST
ORGVEP_BCF_HGDP <- read.table("ORGVEP_BCF_HGDP.txt", header = TRUE) #Import masterfile
calc_ORGVEP_BCF_HGDP<-  select(ORGVEP_BCF_HGDP, c(1:4,24,43,47:51,54,57,59:65))  %>% distinct()
#####dDAF Calculation
ancestralORGVEP_BCF_HGDP <- calc_ORGVEP_BCF_HGDP %>%
  mutate(DAF_JH = ifelse(REF == toupper(AA_VEP.e112), AF_JH, ifelse(ALT == toupper(AA_VEP.e112), 1 - AF_JH,
                                                                    ifelse(REF == AA_chimp, AF_JH, ifelse(ALT == AA_chimp, 1 - AF_JH, NA))))) %>%
  mutate(DAF_TM = ifelse(REF == toupper(AA_VEP.e112), AF_TM, ifelse(ALT == toupper(AA_VEP.e112), 1 - AF_TM,
                                                                    ifelse(REF == AA_chimp, AF_TM, ifelse(ALT == AA_chimp, 1 - AF_TM, NA))))) %>%
  mutate(DAF_WES = ifelse(REF == toupper(AA_VEP.e112), AF_WES, ifelse(ALT == toupper(AA_VEP.e112), 1 - AF_WES,
                                                                    ifelse(REF == AA_chimp, AF_WES, ifelse(ALT == AA_chimp, 1 - AF_WES, NA))))) 

ancestralORGVEP_BCF_HGDP <- ancestralORGVEP_BCF_HGDP %>%
  mutate(DAF_AFR = ifelse(REF_HGDP == toupper(AA_VEP.e112), AF_AFR , ifelse(ALT_HGDP == toupper(AA_VEP.e112), 1 - AF_AFR,
                                                                            ifelse(REF_HGDP == AA_chimp, AF_AFR, ifelse(ALT_HGDP == AA_chimp, 1 - AF_AFR, NA))))) %>%
  mutate(DAF_NONAFR = ifelse(REF_HGDP == toupper(AA_VEP.e112), AF_NONAFR, ifelse(ALT_HGDP == toupper(AA_VEP.e112), 1 - AF_NONAFR,
                                                                                 ifelse(REF_HGDP == AA_chimp, AF_NONAFR, ifelse(ALT_HGDP == AA_chimp, 1 - AF_NONAFR, NA))))) 

ancestralORGVEP_BCF_HGDP$dDAF_JHTM<- abs(ancestralORGVEP_BCF_HGDP$DAF_JH- ancestralORGVEP_BCF_HGDP$DAF_TM)
ancestralORGVEP_BCF_HGDP$dDAF_JHTM<- as.numeric(as.character(format(round(ancestralORGVEP_BCF_HGDP$dDAF_WES, 4), nsmall = 4)))

ancestralORGVEP_BCF_HGDP$dDAF_OA_AFR<- abs(ancestralORGVEP_BCF_HGDP$DAF_WES- ancestralORGVEP_BCF_HGDP$DAF_AFR)
ancestralORGVEP_BCF_HGDP$dDAF_OA_AFR<- as.numeric(as.character(format(round(ancestralORGVEP_BCF_HGDP$dDAF_OA_AFR, 4), nsmall = 4)))
ancestralORGVEP_BCF_HGDP$dDAF_OA_NON<- abs(ancestralORGVEP_BCF_HGDP$DAF_WES- ancestralORGVEP_BCF_HGDP$DAF_NONAFR)
ancestralORGVEP_BCF_HGDP$dDAF_OA_NON<- as.numeric(as.character(format(round(ancestralORGVEP_BCF_HGDP$dDAF_OA_NON, 4), nsmall = 4)))

#export to correct for indels
ancORGVEP_BCF_HGDP <- filter(ancestralORGVEP_BCF_HGDP, !is.na(ancestralORGVEP_BCF_HGDP$DAF_WES))
#correct_ancORGVEP_BCF_HGDP <- filter(ancestralORGVEP_BCF_HGDP, is.na(ancestralORGVEP_BCF_HGDP$DAF_WES))
#write.table(correct_ancORGVEP_BCF_HGDP, "correct_ancORGVEP_BCF_HGDP.txt", quote = FALSE, row.names = FALSE, sep = "\t")
indelcorrected <- read.table("correct_ancORGVEP_BCF_HGDP.txt", header=TRUE) 
corrANC_ORGVEPBCF_HGDP <- rbind(ancORGVEP_BCF_HGDP,indelcorrected)

#####FST Calculation
fst_snp <- read.table("ORG-genohweking_snp.JH.TM.fst.var",header = TRUE)
fstx<- read.table("ORG-genohweking_snp.x.JH.TM.fst.var",header = TRUE)
fstindel <- read.table("ORG-genohweking_indel.JH.TM.fst.var",header = TRUE)
fstall <- rbind(fst_snp,fstx,fstindel)
corrANCFST_ORGVEPBCF_HGDP <- merge(x=corrANC_ORGVEPBCF_HGDP, y=select(fstall, 3,5), by.x ="Location", by.y = "ID", all.x = TRUE)


#############################################Masterfile: ORGVEP_BCF_HGDP and DAF, FST
simpcorrANCFST_ORGVEPBCF_HGDP <- select(corrANCFST_ORGVEPBCF_HGDP, c(1,21:29))
ORGVEP_BCF_HGDP_CALCS <- merge(x=ORGVEP_BCF_HGDP, y=simpcorrANCFST_ORGVEPBCF_HGDP, by.x ="Location", by.y = "Location", all.x =  TRUE)
#write.table(ORGVEP_BCF_HGDP_CALCS, "ORGVEP_BCF_HGDP_CALCS.txt", quote = FALSE, row.names = FALSE, sep = "\t")


#####DAF FST: quantile and outliers
ddaf <- filter(corrANCFST_ORGVEPBCF_HGDP,corrANCFST_ORGVEPBCF_HGDP$AN >= 90 & !is.na(corrANCFST_ORGVEPBCF_HGDP$dDAF_JHTM )) #genotype rate
quantile(ddaf$dDAF_JHTM, c(.99), na.rm = TRUE)
upperdaf_JHTM <- filter(ddaf,ddaf$dDAF_JHTM> 0.3732864) 
breaks <- c(-0.1, 0.1, 0.2 ,0.3, 0.4,0.5)
pivotdaf <- cut(ddaf$dDAF_JHTM, breaks= breaks, labels = NULL)
daf <- data.frame(Range=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5"), 
                  Value="∆DAF_JH-TM", SNPs=c(summary(pivotdaf)))

hfst <- filter(corrANCFST_ORGVEPBCF_HGDP,corrANCFST_ORGVEPBCF_HGDP$AN >= 90 & !is.na(corrANCFST_ORGVEPBCF_HGDP$AF_WES))
cleanhfst <- hfst %>% filter(!grepl("NaN", HUDSON_FST)) 
quantile(cleanhfst$HUDSON_FST, c(.99), na.rm = TRUE)
upperfst_JHTM <- filter(cleanhfst,cleanhfst$HUDSON_FST> 0.2591648) 
breaks2 <- c(-0.1, 0,0.1, 0.2 ,0.3, 0.4, 0.5)
pivotfst <- cut(cleanhfst$HUDSON_FST, breaks= breaks2, labels = NULL)
fst <- data.frame(Range=c("< 0","0-0.1","0.1-0.2","0.2-0.3", "0.3-0.4","0.4-0.5"),
                  Value="FST_JH-TM", SNPs=c(summary(pivotfst)))          
#####DAF FST: plot
library(ggplot2)
library(ggpubr)
ggplot(rbind(daf,fst), aes(x =Range, y = SNPs, fill = Value, colour = Value)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9)+ theme_bw(base_size = 32)+ labs(y= "No. of Variants", x = "Bins")+
  geom_text(aes(label = SNPs), vjust = -0.5, size=7, colour = "black",position = position_dodge(.9))+
  facet_grid(cols = vars(Value), scales="free_x")+ 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#ggsave("DAF_FST.jpeg", width=50, height=30, units = "cm", bg = "white")

#HGDP-OA
ddaf_OAHGDP <- filter(corrANCFST_ORGVEPBCF_HGDP, corrANCFST_ORGVEPBCF_HGDP$AN >= 90 & !is.na(corrANCFST_ORGVEPBCF_HGDP$dDAF_OA_AFR) ) #genotype rate
quantile(ddaf_OAHGDP$dDAF_OA_AFR, c(.99), na.rm = TRUE)
upperdaf_OA_AFR <- filter(ddaf_OAHGDP,ddaf_OAHGDP$dDAF_OA_AFR> 0.6098) 
breaks3 <- c(-0.1, 0.1, 0.2 , 0.3, 0.4, 0.5, 0.6,0.9)
pivotdaf_OA_AFR <- cut(ddaf_OAHGDP$dDAF_OA_AFR, breaks= breaks3, labels = NULL)
daf_OA_AFR <- data.frame(Range=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.9*"),
                  Comparison="∆DAF_OA-AFR",       
                  Value="∆DAF", SNPs=c(summary(pivotdaf_OA_AFR)))

quantile(ddaf_OAHGDP$dDAF_OA_NON, c(.99), na.rm = TRUE)
upperdaf_OA_NON <- filter(ddaf_OAHGDP,ddaf_OAHGDP$dDAF_OA_NON> 0.28742) 
breaks4 <- c(-0.1, 0.1, 0.2 , 0.3, 0.4)
pivotdaf_OA_NON <- cut(ddaf_OAHGDP$dDAF_OA_NON, breaks= breaks4, labels = NULL)
daf_OA_NON <- data.frame(Range=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"),
                         Comparison="∆DAF_OA-NON",
                         Value="∆DAF", SNPs=c(summary(pivotdaf_OA_NON)))


ggplot(rbind(daf_OA_AFR,daf_OA_NON), aes(x =  Range, y = SNPs, fill = Value, colour = Value)) +geom_col (position = position_dodge2(preserve = "single"))+
  geom_bar(stat = "identity", position = "dodge", width = 0.6)+theme_bw(base_size = 32)+labs(y= "No. of Variants", x = "Bins")+
  geom_text(aes(label = SNPs), vjust = -0.5, size=7, colour = "black",position = position_dodge(.9))+
  facet_grid(cols = vars(Comparison), scales="free_x",space = "free" )+ 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line.y.right = )
#ggsave("DAF_FST-OAHGDP.jpeg", width=60, height=30, units = "cm", bg = "white")

allupper <- rbind(
  data.frame(upperdaf_JHTM,Comparison="∆DAF_JH-TM"),
  data.frame(upperfst_JHTM,Comparison="FST_JH-TM"),
  data.frame(upperdaf_OA_AFR,Comparison="∆DAF_OA_AFR"),
  data.frame(upperdaf_OA_NON,Comparison="∆DAF_OA-NON"))

simpallupper <- select(allupper, c(1,30))
all_upperORGVEP_BCF_HGDP_CALCS <- merge(x=simpallupper, y=ORGVEP_BCF_HGDP_CALCS, by.x ="Location", by.y = "Location", all.x =  TRUE)
filter(all_upperORGVEP_BCF_HGDP_CALCS,all_upperORGVEP_BCF_HGDP_CALCS$IMPACT == "MODERATE")

select(allupper, 1) %>% distinct()
hardysnp <- read.table("OA-homhet-snp.hardy", header = TRUE)
hardyindel <- read.table("OA-homhet-indel.hardy", header = TRUE)
hardyall <- rbind(hardysnp,hardyindel)
upperORGVEP_BCF_HGDP_CALCS_hardy <- merge(all_upperORGVEP_BCF_HGDP_CALCS, hardyall, by.x ="Location", by.y = "ID", all.x =  TRUE)
hardysnphgdp <- read.table("hgdporg-auto-homhet.hardy", header = TRUE)
upperORGVEP_BCF_HGDP_CALCS_allhardy <- merge(upperORGVEP_BCF_HGDP_CALCS_hardy, hardysnphgdp, by.x ="Location", by.y = "ID_HGDP", all.x =  TRUE)

export_upperORGVEPBCF_HGDP_CALCSallhardy <- select(upperORGVEP_BCF_HGDP_CALCS_allhardy, !c(16,17,28,29,32,34))
write.table(export_upperORGVEPBCF_HGDP_CALCSallhardy,"upperORGVEP_BCF_HGDP_CALCS_allhardy.txt", quote = FALSE, row.names = FALSE, sep = "\t")


41-9




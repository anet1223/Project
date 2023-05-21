library(readxl)
library(dplyr)
library(ggplot2)
library(multcomp)
library(plotrix)


#Извлечение данных из папки
dir1 <- this.dir()
setwd(dir1)

# чтение данных
raw_control1 <- read.csv('1T.thHB27_R23I_0.csv')
raw_control2 <- read.csv('2T.thHB27_R23I_0.csv')
raw_control3 <- read.csv('3T.thHB27_R23I_0.csv')
raw_control4 <- read.csv('4T.thHB27_R23I_0.csv')

raw_exp20_1 <- read.csv('8T.thHB27_R23I_20.csv')
raw_exp20_2 <- read.csv('9T.thHB27_R23I_20.csv')
raw_exp20_3 <- read.csv('10T.thHB27_R23I_20.csv')

raw_exp50_1 <- read.csv('11T.thHB27_R23I_50.csv')
raw_exp50_2 <- read.csv('12T.thHB27_R23I_50.csv')
raw_exp50_3 <- read.csv('13T.thHB27_R23I_50.csv')

raw_exp100_1 <- read.csv('5T.thHB27_R23I_100.csv')
raw_exp100_2 <- read.csv('6T.thHB27_R23I_100.csv')
raw_exp100_3 <- read.csv('7T.thHB27_R23I_100.csv')

# фильтрация и отбор столбцов
control1 <- raw_control1 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass) %>% mutate(Condition='control1')%>% mutate(Type='control')

control2 <- raw_control2 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='control2')%>% mutate(Type='control')

control3 <- raw_control3 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='control3')%>% mutate(Type='control')

control4 <- raw_control4 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='control4')%>% mutate(Type='control')


exp20_1 <- raw_exp20_1 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp20_1')%>% mutate(Type='drug')

exp20_2 <- raw_exp20_2 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp20_2')%>% mutate(Type='drug')

exp20_3 <- raw_exp20_3 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp20_3')%>% mutate(Type='drug')

exp50_1 <- raw_exp50_1 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp50_1')%>% mutate(Type='drug')

exp50_2 <- raw_exp50_2 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp50_2')%>% mutate(Type='drug')

exp50_3 <- raw_exp50_3 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp50_3')%>% mutate(Type='drug')

exp100_1 <- raw_exp100_1 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp100_1')%>% mutate(Type='drug')

exp100_2 <- raw_exp100_2 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp100_2')%>% mutate(Type='drug')

exp100_3 <- raw_exp100_3 %>% filter(!grepl('#DECOY#',Accession)) %>% select(Accession, Score...., Coverage...., X.Peptides, X.Unique,PTM,Avg..Mass,Description)%>%rename(ID = Accession, Score = Score...., Coverage = Coverage...., Peptides = X.Peptides, Unique = X.Unique, Avg_Mass = Avg..Mass)%>% mutate(Condition='exp100_3')%>% mutate(Type='drug')

# Объединение данных

df <- rbind(control1,control2,control3,control4,exp20_1, exp20_2, exp20_3,
            exp50_1, exp50_2, exp50_3, exp100_1, exp100_2, exp100_3)
df$ID <- as.factor(df$ID)
df$Condition <- as.factor(df$Condition)
df$Type <- as.factor(df$Type)



df$Protein_function <- unlist(strsplit(df$Description," OS="))[c(seq(from = 1,to = 12066, by = 2))]
df$Description <- NULL
df <- df %>% filter(Score>=90)


# создать таблицу для РСА

c1 <- subset(df,Condition =='control1',select = c(ID,Peptides))
colnames(c1)[2] <-'c1' 
c2 <- subset(df,Condition =='control2',select = c(ID,Peptides))
colnames(c2)[2] <-'c2' 
c3 <- subset(df,Condition =='control3',select = c(ID,Peptides))
colnames(c3)[2] <-'c3' 
c4 <- subset(df,Condition =='control4',select = c(ID,Peptides))
colnames(c4)[2] <-'c4' 

d20_1 <- subset(df,Condition =='exp20_1',select = c(ID,Peptides))
colnames(d20_1)[2] <-'d20_1' 
d20_2 <- subset(df,Condition =='exp20_2',select = c(ID,Peptides))
colnames(d20_2)[2] <-'d20_2' 
d20_3 <- subset(df,Condition =='exp20_2',select = c(ID,Peptides))
colnames(d20_3)[2] <-'d20_3' 

d50_1 <- subset(df,Condition =='exp50_1',select = c(ID,Peptides))
colnames(d50_1)[2] <-'d50_1' 
d50_2 <- subset(df,Condition =='exp50_2',select = c(ID,Peptides))
colnames(d50_2)[2] <-'d50_2' 
d50_3 <- subset(df,Condition =='exp50_2',select = c(ID,Peptides))
colnames(d50_3)[2] <-'d50_3' 

d100_1 <- subset(df,Condition =='exp100_1',select = c(ID,Peptides))
colnames(d100_1)[2] <-'d100_1' 
d100_2 <- subset(df,Condition =='exp100_2',select = c(ID,Peptides))
colnames(d100_2)[2] <-'d100_2' 
d100_3 <- subset(df,Condition =='exp100_2',select = c(ID,Peptides))
colnames(d100_3)[2] <-'d100_3' 

by <- join_by(ID)
new_table <- full_join(c1,c2,by)
new_table <- full_join(new_table,c3,by)
new_table <- full_join(new_table,c4,by)
new_table <- full_join(new_table,d20_1,by)
new_table <- full_join(new_table,d20_2,by)
new_table <- full_join(new_table,d20_3,by)
new_table <- full_join(new_table,d50_1,by)
new_table <- full_join(new_table,d50_2,by)
new_table <- full_join(new_table,d50_3,by)
new_table <- full_join(new_table,d100_1,by)
new_table <- full_join(new_table,d100_2,by)
new_table <- full_join(new_table,d100_3,by)

new_table_na_omit <- na.omit(new_table)
new_table_na_omit$result <- 0

new_table_na_omit$ID <- as.character(new_table_na_omit$ID)

chang_ids <- read.csv('112protein_change.csv')[,2]
chang_ids <- unique(chang_ids)
chang_ids <- chang_ids[chang_ids !='Q72IZ1|Q72IZ1_THET2']
chang_ids <- chang_ids[chang_ids !='sp|P61500|RECA_THET2']

for (i in 1:151){
  if (new_table_na_omit$ID[i]%in%chang_ids){
    new_table_na_omit$result[i] = 1
  }
}


write.csv(new_table_na_omit,'all_exp_na_omit.csv')


proba <- new_table[,-1]
row.names(proba) <- new_table$ID
proba_scaled <- as.data.frame(t(apply(proba,1, function(x) scale(x))))
colnames(proba_scaled) <- colnames(proba)
write.csv(proba_scaled,'proba_scaled.csv')

colnames(sorted.df)[1] <- 'ID'
row.names(sorted.df) <- NULL
df_fun <- unique(df[,c(1,10)])
s_id <- unique(sorted.df$ID)
a <- subset(df_fun, ID %in% s_id)
sorted_proba <- full_join(sorted.df,a,by)

# Удаление кто не во всех группах

df$Type <- as.character(df$Type)
df$Type[df$Condition=='exp20_1'] <- 'drug20'
df$Type[df$Condition=='exp20_2'] <- 'drug20'
df$Type[df$Condition=='exp20_3'] <- 'drug20'
df$Type[df$Condition=='exp50_1'] <- 'drug50'
df$Type[df$Condition=='exp50_2'] <- 'drug50'
df$Type[df$Condition=='exp50_3'] <- 'drug50'
df$Type[df$Condition=='exp100_1'] <- 'drug100'
df$Type[df$Condition=='exp100_2'] <- 'drug100'
df$Type[df$Condition=='exp100_3'] <- 'drug100'
df$Type <- as.factor(df$Type)
id <- unique(df$ID)

nou <- data.frame(ID=id, Control = NA, Drug20 = NA, Drug50 = NA, Drug100 = NA)

for (i in 1:444){
  nou$Control[i] <- nrow(df[df$ID==id[i]& df$Type=='control',])
  nou$Drug20[i] <- nrow(df[df$ID==id[i]& df$Type=='drug20',])
  nou$Drug50[i] <- nrow(df[df$ID==id[i]& df$Type=='drug50',])
  nou$Drug100[i] <- nrow(df[df$ID==id[i]& df$Type=='drug100',])
}


nul_id <- nou$ID[nou$Drug20 ==0 | nou$Control ==0 | nou$Drug50 ==0 | nou$Drug100 ==0]

new_full <- df[!df$ID%in%nul_id,]
new_full$ID <- as.character(new_full$ID)
new_full$ID <- as.factor(new_full$ID)

table_sum_full <- new_full %>%
  group_by(ID,Type)%>%
  summarise(mean_peptides = round(mean(Peptides),2),
            sd = round(sd(Peptides),2))

table_drug20_full <- table_sum_full[table_sum_full$Type == 'drug20',]

table_drug50_full <- table_sum_full[table_sum_full$Type == 'drug50',]

table_drug100_full <- table_sum_full[table_sum_full$Type == 'drug100',]

table_control_full <- table_sum_full[table_sum_full$Type == 'control',]
colnames(table_drug20_full) <- c("ID","Type","mean_peptides20",'sd20')
colnames(table_drug50_full) <- c("ID","Type","mean_peptides50",'sd50')
colnames(table_drug100_full) <- c("ID","Type","mean_peptides100",'sd100')

table_list3 <- list(table_control_full[,c(1,3,4)],table_drug20_full[,c(1,3,4)],
                    table_drug50_full[,c(1,3,4)],
                    table_drug100_full[,c(1,3,4)])
table_wide_full <- Reduce(function(x, y) merge(x, y, all=TRUE), table_list3)

table_wide_full[is.na(table_wide_full)] <- 0

table_wide_full$step1 <- NA
table_wide_full$step2 <- NA
table_wide_full$step3 <- NA



for (i in 1:202){
  if ((table_wide_full$mean_peptides[i]-3*table_wide_full$sd[i])
      >(table_wide_full$mean_peptides20[i]+3*table_wide_full$sd20[i])){
    table_wide_full$step1[i] <- 'down'
  }
  else if ((table_wide_full$mean_peptides[i]+3*table_wide_full$sd[i])
           <(table_wide_full$mean_peptides20[i]-3*table_wide_full$sd20[i])){
    table_wide_full$step1[i] <- 'up'
  }
  else (table_wide_full$step1[i] <- 'level')
}

for (i in 1:202){
  if ((table_wide_full$mean_peptides20[i]-3*table_wide_full$sd20[i])
      >(table_wide_full$mean_peptides50[i]+3*table_wide_full$sd50[i])){
    table_wide_full$step2[i] <- 'down'
  }
  else if ((table_wide_full$mean_peptides20[i]+3*table_wide_full$sd20[i])
           <(table_wide_full$mean_peptides50[i]-3*table_wide_full$sd50[i])){
    table_wide_full$step2[i] <- 'up'
  }
  else (table_wide_full$step2[i] <- 'level')
}

for (i in 1:202){
  if ((table_wide_full$mean_peptides50[i]-3*table_wide_full$sd50[i])
      >(table_wide_full$mean_peptides100[i]+3*table_wide_full$sd100[i])){
    table_wide_full$step3[i] <- 'down'
  }
  else if ((table_wide_full$mean_peptides50[i]+3*table_wide_full$sd50[i])
           <(table_wide_full$mean_peptides100[i]-3*table_wide_full$sd100[i])){
    table_wide_full$step3[i] <- 'up'
  }
  else (table_wide_full$step3[i] <- 'level')
}

table_wide_full$change_by_Anova <- 'no'

for (i in 1:202){
  if (table_wide_full$ID[i]%in%chang_ids){
    table_wide_full$change_by_Anova[i] = 'yes'
  }
}

table_wide_full$flow <- paste(table_wide_full$step1,table_wide_full$step2,
                              table_wide_full$step3)

table_wide_full$Protein_function <- NA
table_wide_full$ID <- as.character(table_wide_full$ID)
df$ID <- as.character(df$ID)

for (i in 1:202){
  name = table_wide_full$ID[i]
  fun = df$Protein_function[df$ID==name]
  table_wide_full$Protein_function[i] = fun
}

table_wide_full$ID <- as.factor(table_wide_full$ID)
df$ID <- as.factor(df$ID)

flow <- table_wide_full %>%
  group_by(flow)%>%
  summarise(count = n())

write.csv(file = 'table_wide_full.csv', table_wide_full)
write.csv(file = 'flow.csv', flow)

# Удаление лишних значений

id <- unique(df$ID)

nou <- data.frame(ID=id, Control = 0, Drug = 0)

for (i in 1:444){
  nou$Control[i] <- nrow(df[df$ID==id[i]& df$Type=='control',])
  nou$Drug[i] <- nrow(df[df$ID==id[i]& df$Type=='drug',])
}


nul_id <- nou$ID[nou$Drug ==0 | nou$Control ==0]

new <- df[!df$ID%in%nul_id,]
new$ID <- as.factor(new$ID)
new$Type <- as.character(new$Type)
new$Type[new$Condition=='exp20_1'] <- 'drug20'
new$Type[new$Condition=='exp20_2'] <- 'drug20'
new$Type[new$Condition=='exp20_3'] <- 'drug20'
new$Type[new$Condition=='exp50_1'] <- 'drug50'
new$Type[new$Condition=='exp50_2'] <- 'drug50'
new$Type[new$Condition=='exp50_3'] <- 'drug50'
new$Type[new$Condition=='exp100_1'] <- 'drug100'
new$Type[new$Condition=='exp100_2'] <- 'drug100'
new$Type[new$Condition=='exp100_3'] <- 'drug100'

new$Type <- as.factor(new$Type)
new$Type <- factor(new$Type, levels =c('control','drug20','drug50','drug100'),  labels=c('control','drug20','drug50','drug100'))
write.csv(file = '278common_protein.csv', new)

# Сравнение ANOVA

new_unique <- unique(new$ID)
comp_anova <- as.data.frame(new_unique)
comp_anova$p_value <- NA

for (i in 1:278){
  name <- new_unique[i]
  test <- aov(Peptides~Type,new[new$ID==name,])
  p_value <- summary(test)[[1]][["Pr(>F)"]][1]
  if (!is.null(p_value)){
  comp_anova$p_value[comp_anova$new_unique==name] <- p_value}
}

comp_anova <- na.omit(comp_anova)
anova_sign <- comp_anova[comp_anova$p_value<0.05,]
no_change <- comp_anova[comp_anova$p_value>=0.05,]

sign_id <- as.vector(anova_sign$new_unique)
no_change_id <- as.vector(no_change$new_unique)

# Фильтр по значимым

sign_df <- new[new$ID%in%sign_id,]
sign_df <- rbind(sign_df,new[new$ID =='Q72J50|Q72J50_THET2',])
sign_df <- rbind(sign_df,new[new$ID =='Q72GV7|Q72GV7_THET2',])

sign_df$Type <- factor(sign_df$Type, levels =c('control','drug20','drug50','drug100'),  labels=c('control','drug20','drug50','drug100'))

no_change_df <- new[new$ID%in%no_change_id,]
no_change_df$Type <- factor(no_change_df$Type, levels =c('control','drug20','drug50','drug100'),  labels=c('control','drug20','drug50','drug100'))

write.csv(file = '147protein_no_change.csv', no_change_df)
write.csv(file = '112protein_change.csv', sign_df)

# суммарная таблица по различиям
table_sum <- sign_df %>%
  group_by(ID,Type)%>%
  summarise(p_value = 0, mean_peptides = round(mean(Peptides),2), 
            sd_peptides = round(sd(Peptides),2),
            mean_unique = round(mean(Unique),2),
            mean_coverage = round(mean(Coverage),2), mean_score = round(mean(Score),2),
            mass = mean(Avg_Mass))

table_drug20 <- table_sum[table_sum$Type == 'drug20',]
table_drug20$p_value <- sapply(table_drug20$ID,function(x) TukeyHSD(
  aov(Peptides~Type, data = sign_df[sign_df$ID==x,]),
  conf.level=.95)$Type[19])

table_drug50 <- table_sum[table_sum$Type == 'drug50',]
table_drug50$p_value <- sapply(table_drug50$ID,function(x) TukeyHSD(
  aov(Peptides~Type, data = sign_df[sign_df$ID==x,]),
  conf.level=.95)$Type[20])

table_drug100 <- table_sum[table_sum$Type == 'drug100',]
table_drug100$p_value <- sapply(table_drug100$ID,function(x) TukeyHSD(
  aov(Peptides~Type, data = sign_df[sign_df$ID==x,]),
  conf.level=.95)$Type[21])

table_control <- table_sum[table_sum$Type == 'control',]

table_sum <- rbind(table_control,table_drug100,table_drug50,table_drug20)
write.csv(file = 'table_sum.csv', table_sum)

colnames(table_drug20) <- c("ID","Type","p_value20","mean_peptides20","sd_peptides20",
                            "mean_unique20","mean_coverage20","mean_score20","mass")
colnames(table_drug50) <- c("ID","Type","p_value50","mean_peptides50","sd_peptides50",
                            "mean_unique50","mean_coverage50","mean_score50","mass")
colnames(table_drug100) <- c("ID","Type","p_value100","mean_peptides100","sd_peptides100",
                            "mean_unique100","mean_coverage100","mean_score100","mass")


table_list <- list(table_control[,c(1,3,4,5,6,7,8,9)],table_drug20[,c(1,3,4,5,6,7,8)],
                        table_drug50[,c(1,3,4,5,6,7,8)],
                        table_drug100[,c(1,3,4,5,6,7,8)])
table_wide_diff <- Reduce(function(x, y) merge(x, y, all=TRUE), table_list)


# суммарная таблица без изменений
table_no_change <- no_change_df %>%
  group_by(ID,Type)%>%
  summarise(mean_peptides = round(mean(Peptides),2), 
            sd_peptides = round(sd(Peptides),2),
            mean_unique = round(mean(Unique),2),
            mean_coverage = round(mean(Coverage),2), mean_score = round(mean(Score),2),
            mass = mean(Avg_Mass))

table_drug20_no <- table_no_change[table_no_change$Type == 'drug20',]

table_drug50_no <- table_no_change[table_no_change$Type == 'drug50',]

table_drug100_no <- table_no_change[table_no_change$Type == 'drug100',]

table_control_no <- table_no_change[table_no_change$Type == 'control',]
colnames(table_drug20_no) <- c("ID","Type","mean_peptides20","sd_peptides20",
                            "mean_unique20","mean_coverage20","mean_score20","mass")
colnames(table_drug50_no) <- c("ID","Type","mean_peptides50","sd_peptides50",
                            "mean_unique50","mean_coverage50","mean_score50","mass")
colnames(table_drug100_no) <- c("ID","Type","mean_peptides100","sd_peptides100",
                             "mean_unique100","mean_coverage100","mean_score100","mass")

table_list2 <- list(table_control_no[,c(1,3,4,5,6,7,8)],table_drug20_no[,c(1,3,4,5,6,7)],
                   table_drug50_no[,c(1,3,4,5,6,7)],
                   table_drug100_no[,c(1,3,4,5,6,7)])
table_wide_no_change <- Reduce(function(x, y) merge(x, y, all=TRUE), table_list2)


table_wide_final_diff <- inner_join(table_wide_diff,sign_df[,c(1,10)],keep = F) 
table_wide_final_diff<- table_wide_final_diff[!duplicated(table_wide_final_diff),]
write.csv(file = 'table_wide_diff.csv', table_wide_final_diff)

table_wide_final_no <- inner_join(table_wide_no_change,no_change_df[,c(1,10)],keep = F) 
table_wide_final_no<- table_wide_final_no[!duplicated(table_wide_final_no),]
write.csv(file = 'table_wide_no_change.csv', table_wide_final_no)


# передислокация

table_wide_final_no <- rbind(table_wide_final_no,
                              table_wide_final_diff[table_wide_final_diff$Protein_function == 'DNA gyrase subunit A',-c(2,9,15,21)])

table_wide_final_no <- rbind(table_wide_final_no,
                             table_wide_final_diff[table_wide_final_diff$Protein_function == 'Protein RecA',-c(2,9,15,21)])

table_wide_final_diff <- subset(table_wide_final_diff,table_wide_final_diff$Protein_function != 'Protein RecA')

table_wide_final_diff <- subset(table_wide_final_diff,table_wide_final_diff$Protein_function != 'DNA gyrase subunit A')

table_wide_final_no <- subset(table_wide_final_no,
                             table_wide_final_no$ID != 'Q72J50|Q72J50_THET2')
table_wide_final_no <- subset(table_wide_final_no,
                              table_wide_final_no$ID != 'Q72GV7|Q72GV7_THET2')



write.csv(file = 'table_wide_diff.csv', table_wide_final_diff)
write.csv(file = 'table_wide_no_change.csv', table_wide_final_no)




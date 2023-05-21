library(readxl)
library(dplyr)
library(ggplot2)
library(multcomp)
library(plotrix)


#Извлечение данных из папки
dir1 <- this.dir()
setwd(dir1)

# чтение данных
raw_control1 <- read.table('1THET2_R23I_0.tsv',header = T)
raw_control2 <- read.table('2THET2_R23I_0.tsv',header = T)
raw_control3 <- read.table('3THET2_R23I_0.tsv',header = T)
raw_control4 <- read.table('4THET2_R23I_0.tsv',header = T)

raw_exp20_1 <- read.table('8THET2_R23I_20.tsv',header = T)
raw_exp20_2 <- read.table('9THET2_R23I_20.tsv',header = T)
raw_exp20_3 <- read.table('10THET2_R23I_20.tsv',header = T)

raw_exp50_1 <- read.table('11THET2_R23I_50.tsv',header = T)
raw_exp50_2 <- read.table('12THET2_R23I_50.tsv',header = T)
raw_exp50_3 <- read.table('13THET2_R23I_50.tsv',header = T)

raw_exp100_1 <- read.table('5THET2_R23I_100.tsv',header = T)
raw_exp100_2 <- read.table('6THET2_R23I_100.tsv',header = T)
raw_exp100_3 <- read.table('7THET2_R23I_100.tsv',header = T)

# фильтрация и отбор столбцов
control1 <- raw_control1 %>% mutate(Condition='control1')%>% mutate(Type='control')
control2 <- raw_control2 %>% mutate(Condition='control2')%>% mutate(Type='control')
control3 <- raw_control3 %>% mutate(Condition='control3')%>% mutate(Type='control')
control4 <- raw_control4 %>% mutate(Condition='control4')%>% mutate(Type='control')

exp20_1 <- raw_exp20_1 %>% mutate(Condition='exp20_1')%>% mutate(Type='drug')
exp20_2 <- raw_exp20_2 %>% mutate(Condition='exp20_2')%>% mutate(Type='drug')
exp20_3 <- raw_exp20_3 %>% mutate(Condition='exp20_3')%>% mutate(Type='drug')

exp50_1 <- raw_exp50_1 %>% mutate(Condition='exp50_1')%>% mutate(Type='drug')
exp50_2 <- raw_exp50_2 %>% mutate(Condition='exp50_2')%>% mutate(Type='drug')
exp50_3 <- raw_exp50_3 %>% mutate(Condition='exp50_3')%>% mutate(Type='drug')

exp100_1 <- raw_exp100_1 %>% mutate(Condition='exp100_1')%>% mutate(Type='drug')
exp100_2 <- raw_exp100_2 %>% mutate(Condition='exp100_2')%>% mutate(Type='drug')
exp100_3 <- raw_exp100_3 %>% mutate(Condition='exp100_3')%>% mutate(Type='drug')

# Объединение данных

df <- rbind(control1,control2,control3,control4,exp20_1, exp20_2, exp20_3,
            exp50_1, exp50_2, exp50_3, exp100_1, exp100_2, exp100_3)
colnames(df)[8] <- 'ID'
colnames(df)[4] <- 'Peptides'
df$ID <- as.factor(df$ID)
df$Condition <- as.factor(df$Condition)
df$Type <- as.factor(df$Type)

df <- df %>% filter(Peptides>1)


# Удаление лишних значений

id <- unique(df$ID)

nou <- data.frame(ID=id, Control = 0, Drug = 0)

for (i in 1:384){
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
write.csv(file = '273common_protein.csv', new)

# Сравнение ANOVA

new_unique <- unique(new$ID)
comp_anova <- as.data.frame(new_unique)
comp_anova$p_value <- NA

for (i in 1:273){
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

sign_df$Type <- factor(sign_df$Type, levels =c('control','drug20','drug50','drug100'),  labels=c('control','drug20','drug50','drug100'))

no_change_df <- new[new$ID%in%no_change_id,]
no_change_df$Type <- factor(no_change_df$Type, levels =c('control','drug20','drug50','drug100'),  labels=c('control','drug20','drug50','drug100'))

write.csv(file = '166protein_no_change.csv', no_change_df)
write.csv(file = '91protein_change.csv', sign_df)

# суммарная таблица по различиям
table_sum <- sign_df %>%
  group_by(ID,Type)%>%
  summarise(p_value = 0, mean_peptides = round(mean(Peptides),2), 
            sd_peptides = round(sd(Peptides),2))

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

colnames(table_drug20) <- c("ID","Type","p_value20","mean_peptides20","sd_peptides20")
colnames(table_drug50) <- c("ID","Type","p_value50","mean_peptides50","sd_peptides50")
colnames(table_drug100) <- c("ID","Type","p_value100","mean_peptides100","sd_peptides100")


table_list <- list(table_control[,c(1,3,4,5)],table_drug20[,c(1,3,4,5)],
                        table_drug50[,c(1,3,4,5)],
                        table_drug100[,c(1,3,4,5)])
table_wide_diff <- Reduce(function(x, y) merge(x, y, all=TRUE), table_list)

peaks_dif <- read.csv('C:/Users/anja1/OneDrive/Documents/Биоинформатика/Проект/T. thermophilus/Peaks/table_wide_diff.csv')

peaks_no_change <- read.csv('C:/Users/anja1/OneDrive/Documents/Биоинформатика/Проект/T. thermophilus/Peaks/table_wide_no_change.csv')

table_wide_diff$Peaks_dif <- NA
table_wide_diff$Peaks_no_change <- NA
table_wide_diff$ID <- sub("tr\\|",'',table_wide_diff$ID)

for (i in 1:91){
  if (table_wide_diff$ID[i]%in%peaks_dif$ID){
    table_wide_diff$Peaks_dif[i] <- 'yes' 
  }
}
for (i in 1:91){
  if (table_wide_diff$ID[i]%in%peaks_no_change$ID){
    table_wide_diff$Peaks_no_change[i] <- 'yes' 
  }
}

sum(table_wide_diff$Peaks_dif == 'yes',na.rm = T)
#64
sum(table_wide_diff$Peaks_no_change == 'yes',na.rm = T)
#21
write.csv(file = 'table_wide_diff.csv', table_wide_diff)

# суммарная таблица без изменений
table_no_change <- no_change_df %>%
  group_by(ID,Type)%>%
  summarise(mean_peptides = round(mean(Peptides),2), 
            sd_peptides = round(sd(Peptides),2))

table_drug20_no <- table_no_change[table_no_change$Type == 'drug20',]

table_drug50_no <- table_no_change[table_no_change$Type == 'drug50',]

table_drug100_no <- table_no_change[table_no_change$Type == 'drug100',]

table_control_no <- table_no_change[table_no_change$Type == 'control',]
colnames(table_drug20_no) <- c("ID","Type","mean_peptides20","sd_peptides20")
colnames(table_drug50_no) <- c("ID","Type","mean_peptides50","sd_peptides50")
colnames(table_drug100_no) <- c("ID","Type","mean_peptides100","sd_peptides100")

table_list2 <- list(table_control_no[,c(1,3,4)],table_drug20_no[,c(1,3,4)],
                   table_drug50_no[,c(1,3,4)],
                   table_drug100_no[,c(1,3,4)])
table_wide_no_change <- Reduce(function(x, y) merge(x, y, all=TRUE), table_list2)



table_wide_no_change$Peaks_dif <- NA
table_wide_no_change$Peaks_no_change <- NA
table_wide_no_change$ID <- sub("tr\\|",'',table_wide_no_change$ID)



for (i in 1:166){
  if (table_wide_no_change$ID[i]%in%peaks_no_change$ID){
    table_wide_no_change$Peaks_no_change[i] <- 'yes' 
  }
}

for (i in 1:166){
  if (table_wide_no_change$ID[i]%in%peaks_dif$ID){
    table_wide_no_change$Peaks_dif[i] <- 'yes' 
  }
}

sum(table_wide_no_change$Peaks_no_change == 'yes',na.rm = T)
#90
sum(table_wide_no_change$Peaks_dif == 'yes',na.rm = T)
#43

write.csv(file = 'table_wide_no_change.csv', table_wide_no_change)

# отбор общих для peaks и identypy белков

big_dif <- table_wide_diff %>% filter(Peaks_dif=='yes')
big_dif$impact <- 'difference'
big_peaks <- peaks_dif %>% filter(ID%in%big_dif$ID)
big_dif_peaks <- 
inner_join(big_dif[,c(1,3,4,6,7,9,10,12,13,16)],
           big_peaks[,c(2,4,5,11,12,17,18,23,24,28)], by='ID')


big_no <- table_wide_no_change %>% filter(Peaks_no_change=='yes')
big_no$impact <- 'no_change'
big_no_peaks <- peaks_no_change %>% filter(ID%in%big_no$ID)
big_no_change <- inner_join(big_no[,c(1:9,12)],
                            big_no_peaks[,c(2,3,4,9,10,14,15,19,20,24)], by='ID')

big_table <- rbind(big_dif_peaks,big_no_change)

write.csv(file = 'combine_identiry_peaks_proteins.csv', big_table)


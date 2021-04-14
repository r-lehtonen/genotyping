# genotipagem e equilíbrio de HW com o pacote genetics
# arquivo de entrada: as colunas com os SNPs devem 
# ser sequenciais
# os genótipos devem estar com letras e separados 
# por barra: ex: G/T

library(genetics)
library(dplyr)

dados<-read.table("genotipagem.txt", header = T)
caso<-filter(dados, grupo==1)   
controle<-filter(dados, grupo==2) 


#cria a tabela de saida
tabela.saida<-data.frame(grupo=as.character(NA),SNP=as.character(NA),genot1=as.character(NA), N1=NA, F1=NA,genot2=as.character(NA), N2=NA, F2=NA,genot3=as.character(NA),N3=NA,F3=NA,alelo1=as.character(NA),N4=NA,F4=NA,alelo2=as.character(NA),N5=NA,F5=NA,hw=NA,qui_genot=NA,qui_alel=NA,stringsAsFactors = F)
tabela.saida

head(dados)
nomes<-colnames(dados) #nomes dos SNPs
i<-2 #número da coluna inicial para genotipar
snp<-7 # número da coluna final para genotipar
j<-1

for (i in i:snp) {
  nome_snp<-nomes[i]
  genotipos1<-na.omit(genotype(caso[,i], sep="/"))#genotipagem
  resultado1<-summary(genotipos1)#sumariza resultados
  eqhw1<-HWE.chisq(genotipos1)#equilíbrio de HW
  nome_genot<-rownames(resultado1$genotype.freq)
  nome_alelos<-rownames(resultado1$allele.freq)
  
  
  tabela.saida[j,]<-data.frame(grupo=as.character("caso"),SNP=as.character(nome_snp),genot1=as.character(nome_genot[1]), N1=resultado1$genotype.freq[1,1], F1=resultado1$genotype.freq[1,2],genot2=as.character(nome_genot[2]), N2=resultado1$genotype.freq[2,1],F2=resultado1$genotype.freq[2,2], genot3=as.character(nome_genot[3]),N3=resultado1$genotype.freq[3,1],F3=resultado1$genotype.freq[3,2],alelo1=as.character(nome_alelos[1]),N4=resultado1$allele.freq[1,1],F4=resultado1$allele.freq[1,2],alelo2=as.character(nome_alelos[2]),N5=resultado1$allele.freq[2,1],F5=resultado1$allele.freq[2,2],hw=eqhw1$p.value,qui_genot=NA,qui_alel=NA,stringsAsFactors = F)
  j<-j+1
  
  genotipos2<-na.omit(genotype(controle[,i], sep="/"))#genotipagem
  resultado2<-summary(genotipos2)#sumariza resultados
  eqhw2<-HWE.chisq(genotipos2)#equilíbrio de HW
  
  #qui-quadrado alelo
  a<-resultado1$allele.freq[1,1]
  b<-resultado1$allele.freq[2,1]
  c<-resultado2$allele.freq[1,1]
  d<-resultado2$allele.freq[2,1]
  quadrado<-rbind(c(a,c),c(b,d))
  quadrado
  res_qui<-chisq.test(quadrado,correct=F)
  
  #qui-quadrado genotipo
  g<-resultado1$genotype.freq[1,1]
  h<-resultado1$genotype.freq[2,1]
  l<-resultado1$genotype.freq[3,1]
  k<-resultado2$genotype.freq[1,1]
  m<-resultado2$genotype.freq[2,1]
  o<-resultado2$genotype.freq[3,1]
  quadrado2<-rbind(c(g,k),c(h,m),c(l,o))
  res_qui2<-chisq.test(quadrado2,correct=F)
  
  tabela.saida[j,]<-data.frame(grupo=as.character("controle"),SNP=as.character(nome_snp),genot1=as.character(nome_genot[1]), N1=resultado2$genotype.freq[1,1], F1=resultado2$genotype.freq[1,2],genot2=as.character(nome_genot[2]), N2=resultado2$genotype.freq[2,1],F2=resultado2$genotype.freq[2,2], genot3=as.character(nome_genot[3]),N3=resultado2$genotype.freq[3,1],F3=resultado2$genotype.freq[3,2],alelo1=as.character(nome_alelos[1]),N4=resultado2$allele.freq[1,1],F4=resultado2$allele.freq[1,2],alelo2=as.character(nome_alelos[2]),N5=resultado2$allele.freq[2,1],F5=resultado2$allele.freq[2,2],hw=eqhw2$p.value,qui_genot=res_qui2$p.value,qui_alel=res_qui$p.value,stringsAsFactors = F)
  
  
  
  i<-(i+1)
  j<-j+1
}

head(tabela.saida)

library(openxlsx)
write.xlsx(tabela.saida, file="genotipagem.xlsx")
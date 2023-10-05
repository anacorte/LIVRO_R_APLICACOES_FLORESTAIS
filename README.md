# LIVRO_R_APLICACOES_FLORESTAIS
# 11 Teste de hipóteses
## 11.1 Teste t-Student

# A seguir, considere os diâmetros de árvores oriundos de 
# diferentes amostras:
dap1 <- c(30.5,35.3,33.2,40.8,42.3,41.5,36.3,43.2,34.6,38.5)
dap2 <- c(28.2,35.1,33.2,35.6,40.2,37.4,34.2,42.1,30.5,38.4)

## 11.1.1 Teste t para uma média

t.test(dap1, # amostra a ser testada
       mu=35, # hipótese de nulidade
       alternativa = "greater") # teste unilateral à direita

## 11.1.2 Teste t para as médias de duas amostras independentes

t.test(dap1, dap2, # amostras a serem comparadas
       conf.level = 0.99) # nível de significância

## 11.1.3 Teste t para médias de duas amostras dependentes
t.test(dap1, dap2,
       conf.level = 0.99, # nível de significância
       paired = TRUE) # afirma dependência entre as amostras

# Fim-----------


# 12. Aplicações Florestais

# Instalar pacotes necessários
#install.packages("data.table")
#install.packages("car")
#install.packages("lmtest")
#install.packages("faraway")

# Carregar pacotes
library(data.table)
library(car)
library(lmtest)
library(faraway)

## 12.1. Análise de Variância (ANOVA)
### 12.1.1. Delineamento Inteiramente Casualizado (DIC)
parica <- fread(file = "Data/parica.txt", stringsAsFactors = T)
print(parica)

# Resumo
summary(parica)

# média de altura em cada tratamento?
media <- tapply(parica$Rep, parica$Trat, mean) 
print(media)

# variância da altura em cada tratamento?
var <- tapply(parica$Rep, parica$Trat, var)
print(var)

bartlett.test(parica$Rep, parica$Trat)  # homogeneidade de variâncias

# Gráfico: Tratamentos versus Alturas das plântulas
par(mar = c(4.5,3,1.5,1), mgp = c(2,1,0), mfrow = c(1,2))

# plot: Rep x Trat
plot(parica$Trat, parica$Rep, type="o", 
     main="Schizolobium parahyba",
     xlab="Tratamentos", ylab="Altura (cm)")

points(media, pch=20, col=2, cex=1.5)      # Adiciona médias/tratamento

# Usando boxplot
boxplot(parica$Rep ~ parica$Trat, main="Schizolobium parahyba",
        xlab="Tratamentos", ylab="Altura (cm)")
points(media, pch=20, col=2, cex=1.5)

# ANOVA - Delineamento Inteiramente Casualizado (DIC)
anova.DIC <- aov(Rep~Trat,data=parica)
anova(anova.DIC)

# Resultados da ANOVA - DIC
par(mfrow=c(2,2))
plot(anova.DIC)       # análise dos resíduos

# Teste de Tukey
par(mfrow=c(1,1))
Tukey <- TukeyHSD(anova.DIC)
plot(Tukey)

# Pressupostos da ANOVA
## Teste de Normalidade - shapiro-wilk
shapiro.test(resid(anova.DIC))

# Normal Q-Q plot
qqnorm (resid(anova.DIC))       # obtendo o papel de probabilidade normal.
qqline(resid(anova.DIC))        # inserindo uma linha auxiliar (linear).

# Homogeneidade de variâncias - Teste de Levene
leveneTest(Rep~Trat,parica)

#----------------------------------------------------------------------------
# 12.2 Regressão Linear
teca <- fread("Data/Tectona.csv", stringsAsFactors = T)

## Gráficos de dispersão
par(mar = c(4.5,3.5,1.5,1), mgp = c(2,1,0), mfrow = c(1,2))

plot(DAP,Volume, type = "p", main=NULL, font.main=NULL, col.main=NULL, 
     xlab="DAP (cm)", ylab=expression(Volume~(m^3)), font.lab=1, col.lab="black",
     font.axis=1, col.axis = "black")

plot(H,Volume, type = "p", main=NULL, font.main=NULL, col.main=NULL, 
     xlab="Altura (m)", ylab=expression(Volume~(m^3)), font.lab=1, col.lab="black",
     font.axis=1, col.axis = "black")

## Ajuste de modelos
Berkhout <- lm(Volume ~ DAP, data=teca)                 # Berkhout
KGehrardt <- lm(Volume ~ I(DAP^2), data=teca)           # Kopezky-Gehrardt
SHall <- lm(log(Volume) ~ log(DAP) + log(H), data=teca) # Shumacher-Hall
print(Berkhout); print(KGehrardt); print(SHall)

### Resumo dos ajustes
summary(Berkhout)
summary(KGehrardt)
summary(SHall)

### ANOVA da Regressão
anova(Berkhout)
anova(KGehrardt)
anova(SHall)

### Reta ajustada - Berkhout
par(mar = c(4.5,3.5,1.5,1), mgp = c(2,1,0), mfrow = c(1,1))

plot(teca$Volume~teca$DAP, main="Berkhout", 
     xlab="DAP (cm)", ylab = expression(Volume~(m^3)))
abline(Berkhout,lty=2, col="red")

# Predições e Resíduos - Berkhout
predict(Berkhout)
residuals(Berkhout)

### Análise de Resíduos
#### Normalidade dos resíduos
shapiro.test(Berkhout$residuals)
shapiro.test(KGehrardt$residuals)
shapiro.test(SHall$residuals)

#### Autocorrelação de resíduos
durbinWatsonTest(Berkhout)
durbinWatsonTest(KGehrardt)
durbinWatsonTest(SHall)

#### Heterocedasticidade dos resíduos
bptest(Berkhout)
bptest(KGehrardt)
bptest(SHall)

#### Gráficos de resíduos
# Berkhout
par(mar = c(4.5,3,1.5,1), mgp = c(2,1,0), mfrow = c(2,2))
plot(Berkhout)

# KGehrardt
par(mar = c(4.5,3,1.5,1), mgp = c(2,1,0), mfrow = c(2,2))
plot(KGehrardt)

# Schumacher-Hall
par(mar = c(4.5,3,1.5,1), mgp = c(2,1,0), mfrow = c(2,2))
plot(SHall)

### Multicolinearidade - Regressão múltipla
vif(SHall)

# Fim-----------


13. Amostragem Aleatória Simples (AAS)

# Instalar pacotes necessários
#install.packages("data.table")
#install.packages("FinCal")

# Carregar pacotes
library(data.table)
library(FinCal)

# Carregando o conjunto de dados
AAS <- fread("Data/AAS.csv", stringsAsFactors = T)
AAS

# Estatísticas básicas para as unidades amostrais
## 1. Média aritmética
mean(AAS$Volume)

## 2. Variância
var(AAS$Volume)

## 3. Desvio Padrão
sd(AAS$Volume)

## 4. Coeficiente de variação
coefficient.variation(sd=sd(AAS$Volume),
                      avg=mean(AAS$Volume))*100

## 5. Intensidade amostral
E <- function(x){
  media = mean(x)
  E = signif(0.1*media, 4)
  return(E)
}

E(AAS$Volume)

## Número de unidades amostrais possíveis na população
N <- function(A,a){
  N <- A/a
  return(N)
}

N <- N(400000,600)
N

## Determinar se a população é finita ou infinita
FC <- function(x,A,a){
  n <- length(x)
  N <- ceiling(A/a)
  f <- n/N
  FC <- 1-f
  
  if(FC >= 0.98){
    cat("A população é Infinita. Portanto, despreze o FC na fórmula da n.\n")
  }else{
    cat("A população é Finita. Portanto, use o FC para corrigir n.\n")
  }
  return(list(f=f,FC=FC))
}

FC <- FC(x=AAS$Volume, A=400000, a = 600)

# Intensidade de amostragem ideal em função da variância
n <- function(x,A,a){
  N <- ceiling(A/a)
  E = 0.1*mean(x)
  t = qt(1-.05/2, df=length(x)-1)
  n <- ceiling((N*t^2*var(x))/(N*E^2 + t^2*var(x)))
  cat(paste("Para atender ao erro estabelecido você deve amostrar", n, 
            "parcelas.\n"))
  
  if(n <= length(x)){
    cat("Esforço amostral satisfatório. O IF é definitivo!")
  }else{
    cat(paste("Retorne a campo e meça mais", abs(length(x)-n), "parcelas."))
  }
}

n(x = AAS$Volume, A = 400000, a = 600)

## 6. Variância da média
var(AAS$Volume)/length(AAS$Volume)*(FC$FC)

## 7. Erro padrão da média
sbarx <- sd(AAS$Volume)/sqrt(length(AAS$Volume))*(sqrt(FC$FC))
sbarx

## 8. Erro de amostragem
# a) Erro de amostragem absoluto
Ea <- qt(1-.05/2, df=length(AAS$Volume)-1)*sbarx
Ea

# b) Erro de amostragem relativo
Er <- Ea/mean(AAS$Volume)*100
Er

## 9. Intervalo de confiança para média
# Limite inferior para média (LI)
LIbarx <- mean(AAS$Volume)-Ea
LIbarx

# Limite superior para média (LS)
LSbarx <- mean(AAS$Volume)+Ea
LSbarx

## 10. Total da população
hatX <- N*mean(AAS$Volume)
hatX

## 11. Intervalo de confiança para o total
# Limite inferior para total da população (LI)
LIhatx <- hatX - N*Ea
LIhatx

# Limite superior para total da população (LS)
LShatx <- hatX + N*Ea
LShatx

# Usando o pacote "forester" (em construção)
# O pacote "forester" ainda não está disponível no CRAN. 
# Mas, as funções existentes podem ser usadas a partir da 
# versão em desenvolvimento disponível no GitHub.

# Instalar a versão em desenvolvimento...
# remotes::install_github("DeivisonSouza/forester")
library(forester)
RS(x = AAS$Volume, A = 400000, a = 600, LE = 0.1)

# Resultados no painel de visualização
RS(x = AAS$Volume, A = 400000, a = 600, 
   LE = 0.1, DT = TRUE)

# Fim-----------

## 14. Análise fitossociológica

### 14.1 Estrutura Horizontal

# Instalar pacotes necessários
#install.packages("data.table")
#install.packages("ggplot2")
#install.packages("plyr")

# Carregar pacotes
library(data.table)
library(ggplot2)
library(plyr)

# Carregando o conjunto de dados
FOM <- fread("Data/Fito.csv", stringsAsFactors = T)
FOM

# Inspeção dos dados
nrow(FOM)
names(FOM)
dim(FOM)

# n = Número total de indivíduos amostrados na j-ésima parcela
FOM[, .(n=.N), by=Parcela][]

# Número de indivíduos amostrados da i-ésima espécie na j-ésima parcela
FOM[, .(ni=.N), by=c("Parcela", "Especie")]

# Uma visualização gráfica
ggplot(FOM[, .(ni=.N), by=c("Parcela", "Especie")], 
       aes(x=Especie, y=ni, fill=Especie)) + 
  geom_bar(stat="identity",position="dodge",width = 1,colour="black")+
  geom_text(aes(label=ni,hjust=-.3, vjust=0.5),
            position=position_dodge(width = 0.7))+
  facet_grid(~ Parcela, labeller=labeller(
    Parcela = Parcela<-as_labeller(
      c(`1`="Parcela 1",`2`="Parcela 2",`3`="Parcela 3"))))+
  coord_flip()+
  geom_text(data=ddply(.data=FOM, .(Parcela), summarize, 
                       n=paste("n =", length(Especie))), 
            aes(x=23, y=7, label=n), colour="black", 
            inherit.aes=FALSE, parse=FALSE)+
  theme_bw()+
  theme(axis.line.x=element_line(size=0.5,colour="black"),
        axis.line.y=element_line(size=0.5,colour="black"),
        axis.line=element_line(size=1,colour="black"),
        strip.text.x=element_text(colour=1,size=12,family="serif",face="bold"),
        strip.background = element_rect(colour="black", fill="snow2"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black"),
        panel.background=element_blank(),
        axis.text.x=element_text(colour="black",size=12,family="serif",angle=0),
        axis.text.y=element_text(colour=1,size=12,family="serif",face="italic"),
        legend.position="none")+
  scale_x_discrete(name="Espécie")+
  scale_y_continuous(name="Número de indivíduos",
                     limits=c(0,8))

# Parâmetros da estrutura horizontal - Araucaria angustifolia

## 1. Densidade absoluta
DAi <- function(x, A){
  ni <- nrow(subset(FOM, Especie=="Araucaria angustifolia"))
  DAi <- ni/A
  return(DAi)
}

DAi(x = FOM$Especie, A = 0.3)

## 2. Densidade relativa
DRi <- function(x, A){
  ni <- nrow(subset(FOM, Especie=="Araucaria angustifolia"))
  DAi <- ni/A
  DTA <- length(x)/A
  DRi <- (DAi/DTA)*100
  return(DRi)
}

DRi(x = FOM$Especie, A = 0.3)

## 3. Dominância absoluta
DoAi <- function(data, A, ...){
  data <- data[Especie=="Araucaria angustifolia"]
  gi <- data[, .(gi=pi*DAP^2/40000)]
  Gi <- sum(gi)
  DoAi <- Gi/A
  return(DoAi)
}

DoAi(data=FOM, A=0.3)

## 4. Dominância relativa
DoRi <- function(data, A, ...){
  Gt <- data[, .(gi=pi*DAP^2/40000)]
  data <- data[Especie=="Araucaria angustifolia"]
  gi <- data[, .(gi=pi*DAP^2/40000)]
  Gi <- sum(gi)
  DoAi <- Gi/A
  DoRi <- (Gi/sum(Gt))*100
  return(DoRi)
}

DoRi(data=FOM, A=0.3)

## 5. Frequência absoluta
FAi <- function(data, ...){
  Ut <- length(unique(data$Parcela))
  Ui <- unique(data, by=c("Especie", "Parcela"))[, .(Ui=.N), by="Especie"]
  Ui <- Ui[Especie=="Araucaria angustifolia", Ui]
  FAi <- (Ui/Ut)*100
  return(FAi)
}

FAi(data=FOM)

## 6. Frequência relativa
FRi <- function(data, ...){
  Ut <- length(unique(FOM$Parcela))
  Ui <- unique(FOM, by=c("Especie", "Parcela"))[, .(Ui=.N), by="Especie"]
  FAi <- Ui[, .(FAi=(Ui/length(unique(FOM$Parcela)))*100)]
  Ui_AA <- Ui[Especie=="Araucaria angustifolia", Ui]
  FAi_AA <- (Ui_AA/Ut)*100
  FRi <- (FAi_AA /sum(FAi))*100
  return(FRi)
}

FRi(data=FOM)

## 7. Valor de cobertura
VCi <- DRi(x = FOM$Especie, A = 0.3) + DoRi(data=FOM, A=0.3)
VCi

## 8. Porcentagem de cobertura
PCi <- (DRi(x = FOM$Especie, A = 0.3) + DoRi(data=FOM, A=0.3))/2
PCi

## 9. Valor de importância
VIi <- DRi(x = FOM$Especie, A = 0.3) + DoRi(data=FOM, A=0.3) + FRi(data=FOM)
VIi

## 10. Porcentagem de importância
PIi <- (DRi(x = FOM$Especie, A = 0.3) + DoRi(data=FOM, A=0.3) + FRi(data=FOM))/3
PIi

# Uma função genérica
EH <- function(species, sample, d, A,...){
  DT <- data.table(species=species,sample=sample,d=d)
  DT <- DT[,`:=`(gi=pi*d^2/40000)]
  Ui <- unique(DT, by=c("species", "sample"))[, .(Ui=.N), by="species"][order(species)]
  ni <- DT[, .(ni=.N, Gi = sum(gi)), by="species"]
  ni <- ni[Ui,on="species"]
  EH <- ni[,DAi := ni/A,
  ][,DRi := (DAi/sum(DAi))*100,
  ][,DoAi := Gi/A,
  ][,DoRi := (DoAi/sum(DoAi))*100,
  ][,VC := DRi + DoRi,
  ][,PC := VC/2,
  ][,FAi := (Ui/length(unique(DT$sample)))*100,
  ][,FRi := (FAi/sum(FAi))*100,
  ][,VIi := DRi + DoRi + FRi,
  ][,PIi := VIi/3][order(-VIi)]
  return(EH)
}

EH <- EH(species=FOM$Especie, sample=FOM$Parcela, 
         d=FOM$DAP, A=0.3)

# Fim-----------


---
output:
  html_document: default
  pdf_document: default
---

<< Analyze dataset from hospital w. statistic & machine learning >>

Scenario:
- diagnosing phenotypically distinct gall bladder tumor type currently only diagnosible by histopathology. 
- suspecting cancer to be hereditary 
- cohort: unrelated patients suspected to have same gall bladder cancer phenotype
- pilot study: 2 families -> identified 2 variants in genes BMX & BRCA2 (may underlie heritable tumor phenotype)
- AIM: confirm link between genes & heritable tumor to see if they are suitable predictor of phenotype development 

- dataset: clinical genotyping & expression data from patients 
    - tsv file (phenotyping, genotyping, expression data)
    "data_blitz.tsv"
      columns: patient disease phenotype [categorical] unaffected / familial cancer
               karyotype of sex chromosomes [categorical] XX / XY
               genotype of variants rs# in genes BMX, BRCA2 [categorical]
               RNA-seq normalised expression values of BMX, BRCA2 [continuous] 
               
#               

SECTION 1 _ STATISTICS

1. import data
```{r}
setwd("/Users/ellialee/Desktop/Coding/R/Supervision")

library(tidyverse)
data <- read_tsv("data_blitz.tsv")
data
```
2. Report data types of columns (continuous / categorical)
```{r}

data.class <- sapply(data, class)
classOfCol <- data.frame(data.class)
classOfCol
classList <- c(classOfCol[ ,1])
classList

dataType <- c()
for (i in classList) {
  if (i == "character") {dataType <- append(dataType, 'categorical')}
  else {dataType <- append(dataType, 'continuous')}
}
dataType

typeReport <- cbind(classOfCol, dataType)
typeReport

```

3. Perform statistical test to identify if variants in BMX/BRCA2 are associated with patient phenotype
- predictor: expression of BMX/BRCA2 (continuous)
- outcome: phenotype (category)
- Regression model

Data cleaning 
- Phenotype data: 0=unaffected, 1=familial cancer
```{r}
library(tidyverse)
data <- data %>%
  mutate(Phenotype = recode_factor(Phenotype, 'unaffected' = '0', 'familial cancer' = '1'))
data

```
Variants in BMX/BRCA2 to phenotype (both categorical)
- chi-square test (test of frequency)
```{r}
chisq.test(data$rs886040801, data$Phenotype)
chisq.test(data$rs139052738, data$Phenotype)
fisher.test(data$rs139052738, data$Phenotype)
```

++ Expression level to phenotype
Continuous predictor (x) and categorical outcome (y)

```{r}
library(ggplot2)

y <- data$Phenotype

plot(y=y, x=data$`Expression BMX`, xlab='expression of BMX', ylab='phenotype')
plot(y=y, x=data$`Expression BRCA2`, xlab='expression of BRCA2', ylab='phenotype')

```
- expression level of BMX is not correlated with familial cancer
- expression level of BRCA2 is correlated with familial cancer (low expression of BRCA2 lead to familial cancer)





 Model for BRCA2 & phenotype
```{r}

model_regression <- glm(y~data$`Expression BRCA2`, data=data, family='binomial')
summary(model_regression)
  

```
p values significant! reject null hypothesis _ it didn't happen by chance 
BRCA2 expression level associated with familial cancer



4. identify if expression of BMX is affected by rs139052738 (C -> T mutation)

Understanding data:
- The variant is on X-chromosome -> XX have 2 alleles & XY have 1 allele 
- Plan: check dose dependency of BMX in male & female separately and then together. 
    - data cleaning: one copy of variant accounts for one unit '1'



Observe difference in [BMX] between XY w. normal(C) and mutated(T) in Males 
2 samples (2 possible genotypes in male) _ u-test (non-parametric)

```{r}
XY <- data %>%
  filter(Sex_chromosomes == 'XY') %>%
  select(rs139052738, 'Expression BMX')

xBMX <- XY$rs139052738
yBMX <- XY$`Expression BMX`

ggplot(XY, aes(x=xBMX, y=yBMX)) +
  geom_point() +
  xlab("Base at rs139052738") +
  ylab("expression level of BMX")

ggplot(XY, aes(x=xBMX, y=yBMX)) +
  geom_boxplot() +
  xlab("Base at rs139052738") +
  ylab("expression level of BMX")

wilcox.test(yBMX~xBMX)

```
Insignificant; P-value = 0.095 (≥ 0.05) --> Don't reject Null Hypothesis

Expression of BMX is not affected by variant rs139052738 in Males
BMX is not dose dependent.




Check if it's true for females.
more than 2 samples in females _ Kruskal-Wallis test (non-parametric)

Data cleaning needed: 
    Phased data. C/T and T/C are same 
    Possible variables: C/C, C/T, T/C, T/T
    
    C/C = 0
    C/T, T/C = 1
    T/T = 2
    
```{r}
XX <- data %>%
  filter(Sex_chromosomes == 'XX') %>%
  select(rs139052738, 'Expression BMX')

XX <- XX %>%
  mutate(
    rs139052738 = case_when(
      rs139052738 == 'C/C' ~ '0',
      rs139052738 == 'T/T' ~ '2',
      TRUE ~ '1'
    )
  )

```

Checking dose dependency in female 
```{r}
xBMX2 <- XX$rs139052738
yBMX2 <- XX$`Expression BMX`

ggplot(XX, aes(x=xBMX2, y=yBMX2)) +
  geom_point() +
  xlab("number of variant allele") +
  ylab("expression level of BMX")

ggplot(XX, aes(x=xBMX2, y=yBMX2)) +
  geom_boxplot() +
  xlab("number of variant allele") +
  ylab("expression level of BMX")

kruskal.test(yBMX2 ~ xBMX2, data=XX)

```
Insignificant; P-value = 0.58 (≥ 0.05) --> Don't reject Null Hypothesis

Expression of BMX is not affected by variant rs139052738 in Females
BMX is not dose dependent.






Check if it is true when female & male data are integrated:
data cleaning + mutating rs139052738 categorical values 
  0 = homozygous normal
  1 = heterozygou (normal & variant)
  2 = homozygous variant 
  
  for males, 0 = normal, 2=variant as it BMX is not dose dependent 
```{r}
XY2 <- XY %>%
  mutate(
    rs139052738 = case_when(
      rs139052738 == 'C' ~ '0',
      rs139052738 == 'T' ~ '2',
    )
  )

XXandXY_BMX <- rbind(XX, XY2)
yforAll <- XXandXY_BMX$`Expression BMX`
xforAll <- XXandXY_BMX$rs139052738

ggplot(XXandXY_BMX, aes(x=xforAll, y=yforAll)) +
  geom_point() +
  xlab("number of variant allele") +
  ylab("expression level of BMX")

kruskal.test(yforAll ~ xforAll, data=XXandXY_BMX)

```
insignificant; p-value = 0.23 (>0.05) ; don't reject null hypothesis
average level of BMX expression not affected by rs139052738 variant




5. Identify if expression of BRCA2 if affected by rs886040801 (T>G mutation)

Understanding data: 
- All individuals have two copies of allele (somatic)
- one copy of mutation = unit '1' 
    T/T = 0 
    T/G = 1 : data cleaning needed (phased)
    G/G = 2
```{r}
BRCA2_data <- data %>%
  select(rs886040801, 'Expression BRCA2')

BRCA2 <- BRCA2_data %>%
  mutate(
    rs886040801 = case_when(
      rs886040801 == 'T/T' ~ '0',
      rs886040801 == 'G/G' ~ '2',
      TRUE ~ '1'
    )
  )

x_BRCA2 <- BRCA2$rs886040801
y_BRCA2 <- BRCA2$`Expression BRCA2`

ggplot(BRCA2, aes(x=x_BRCA2, y=y_BRCA2)) +
  geom_point() +
  xlab("number of variant allele") +
  ylab("expression level of BRCA2")

kruskal.test(y_BRCA2 ~ x_BRCA2, data=BRCA2)


```
Significant P-Value! rs886040801 affects BRCA2 expression!
BRCA2 is dose dependent; the number of normal allele is proportional to BRCA2 expression.
    one copy of normal = 1/2 BRCA2 expression
    
    


see if phenotypes are classified according to BRCA2 expression:
```{r}
data2 <- read_tsv("data_blitz.tsv")
data2_BRCA <- data2 %>%
  select(Phenotype,rs886040801,'Expression BRCA2')

data2_BRCA <- data2_BRCA %>%
  mutate(
    rs886040801 = case_when(
      rs886040801 == 'T/T' ~ '0',
      rs886040801 == 'G/G' ~ '2',
      TRUE ~ '1'
    )
  )

ggplot(data2_BRCA, aes(x=x_BRCA2, y=y_BRCA2, color = Phenotype)) +
  geom_point() +
  xlab("number of variant allele") +
  ylab("expression level of BRCA2") 

```
Individuals with only variant alleles display familial cancer. 
1/2 expression of variant allele is insufficient for tumor development.





6. Identify if expression of BMX affects BRCA2 expression 
Understanding data: Tests of relationship of two scale data 
```{r}
GeneExpressionData <- data %>%
  select('Expression BMX', 'Expression BRCA2')

BRCA2Level <- GeneExpressionData$`Expression BRCA2`
MBXLevel <- GeneExpressionData$`Expression BMX`

ggplot(GeneExpressionData, aes(x=MBXLevel, y=BRCA2Level))+
  geom_point() +
  geom_smooth(method='lm')

cor.test(BRCA2Level, MBXLevel, method="spearman")

```
insignificant p-value = 0.90 (>0.05) ; do not reject null hypothesis of no correlation 
correlation coefficient close to 0 -> do no reject null hypothesis 
    (H0 rejected if correlation coefficient is extreme (1 or -1))

Expression of BRCA2 and BMX are not correlated. 
  Expression of BMX does not affect expression of BRCA2.









SECTION 2: Machine learning 
use dataset to predict patient phenotype.


1. machine learning approach to visualize data
      
to use 'pair' function, data types have the be numeric 
--> data cleaning 
```{r}
head(data)  #Phenotype already numeric 

# change Sex chromosome to numeric// Number of X-chrom (XY=1, XX=2)
data_sex <- data %>%
  mutate(
    Sex_chromosomes = case_when(
      Sex_chromosomes == 'XX' ~ '2',
      Sex_chromosomes == 'XY'~ '1'
    )
  )

head(data_sex)

# rs139052738 variant _for males, C=0, T=2
data_rs13 <- data_sex %>%
  mutate(
    rs139052738 = case_when(
      rs139052738 == 'C' ~ '0',
      rs139052738 == 'C/C' ~ '0',
      rs139052738 == 'T/T' ~ '2',
      rs139052738 == 'T' ~ '2',
      TRUE ~ '1'
    )
  )
head(data_rs13)

#rs886040801 variant
data_final <- data_rs13 %>%
  mutate(
    rs886040801 = case_when(
      rs886040801 == 'T/T' ~ '0',
      rs886040801 == 'G/G' ~ '2',
      TRUE ~ '1'
    )
  )

head(data_final)

# need to be numeric for machine learning 
for (i in colnames(data_final)) {
  data_final[[i]] <- as.numeric(data_final[[i]])
}

data_final$Phenotype <- data_final$Phenotype-1
head(data_final)

```
Try to visualize correlation:

```{r}
pairs(data_final[1:6], upper.panel = NULL)
```
From the pairs graphs, few possible correlations are present:
I. Change in BRCA2 expression level in unaffected patients & patients with familial cancer
II. Difference in BRCA2 expression level depending on genotypes of variant rs886040801 

rs886040801 genoytpe -> BRCA2 expression -> Phenotype


(((better to have 1,2,3 than from 0,1,2)))


2. deploy a single machine learning method to classify the data based on patient disease phenotype 
      - supervised approaches: Linear regression, Logistic regression, k means clustering, decision trees, neural networks 
      - provide evidence-based explanation why you have chosen your approach.
      - use your approach to predict disease state of individuals in file "predict_data.tsv" (record predictions of these individuals)
      

Machine learning _ supervised _ training the data 
- Making training & test sets
```{r}
split_size = 0.8
sample_size = floor(split_size * nrow(data_final))

set.seed(123)
train_indices <- sample(seq_len(nrow(data_final)), size=sample_size)

train <- data_final[train_indices,]
test <- data_final[-train_indices,]

```


- logistic regression
    - binary output (phenotype affected/familial cancer) _classification
    - data measurements: discrete (genotype) and continuous (expression level)
    - predict likelihood of an event 

```{r}
library(caTools)
library(e1071)

Label.train <- train[,1]
Data.train <- train[,-1]

Label.train <- as.matrix(Label.train)
Data.train <- as.matrix(Data.train)


?LogitBoost
model <-  LogitBoost(Data.train, Label.train)
Data.test <- test
Lab <- predict(model, Data.test, type = 'raw')
data.frame(row.names(test), test$Phenotype, Lab)



```
it is not doing well with predicting phenotypes correctly,,


Decision tree_ y = Expression BRCA2 (continuous)
Good for visualising gene expression 


neural network 

```{r}
library(nnet)
cancer.nnet <- nnet(Phenotype ~ ., data=train, size=4, decay=0.0001, maxit=500, trace=FALSE)
predictions <- predict(cancer.nnet, test[, 1:6])
table(round(predictions), test$Phenotype)
```
Accurate prediction!
- long computational time


using the model for new dataset:
1. getting new dataset cleaned
```{r}
library(tidyverse)
setwd("/Users/ellialee/Desktop/Coding/R/Supervision")
newData <- read_tsv("predict_data.tsv")
head(newData)

# change Sex chromosome to numeric// Number of X-chrom (XY=1, XX=2)
newdata_sex <- newData %>%
  mutate(
    Sex_chromosomes = case_when(
      Sex_chromosomes == 'XX' ~ '2',
      Sex_chromosomes == 'XY'~ '1'
    )
  )


# rs139052738 variant _for males, C=0, T=2
newdata_rs13 <- newdata_sex %>%
  mutate(
    rs139052738 = case_when(
      rs139052738 == 'C' ~ '0',
      rs139052738 == 'C/C' ~ '0',
      rs139052738 == 'T/T' ~ '2',
      rs139052738 == 'T' ~ '2',
      TRUE ~ '1'
    )
  )
head(data_rs13)

#rs886040801 variant
newData_final <- newdata_rs13 %>%
  mutate(
    rs886040801 = case_when(
      rs886040801 == 'T/T' ~ '0',
      rs886040801 == 'G/G' ~ '2',
      TRUE ~ '1'
    )
  )

head(newData_final)

# need to be numeric for machine learning 
for (i in colnames(newData_final)) {
  newData_final[[i]] <- as.numeric(newData_final[[i]])
}

newData_final

```


```{r}
Phenotypes01 <- predict(cancer.nnet, newdata = newData_final)
PredictedPheno <- round(Phenotypes01)
PredictedPheno

BinaryPheno <- as.vector(PredictedPheno)

Phenotype <- c()
for (i in BinaryPheno) {
  if (i == 0) {Phenotype <- append(Phenotype, 'unaffected')}
  else {Phenotype <- append(Phenotype, 'familial cancer')}
}

Phenotype

Final_report <- cbind(Phenotype, newData)
notNeededCol <- c(3,4)
Final_report[, -notNeededCol]




```


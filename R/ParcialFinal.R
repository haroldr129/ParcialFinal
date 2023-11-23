library(devtools)
library(tidyverse)
library(dplyr)
library(rpart)
library(caret)
library(glmnet)

PPI <- DynamicCancerDriverKM::PPI

Matriz_Normal_Tumor <- merge(DynamicCancerDriverKM::BRCA_normal,DynamicCancerDriverKM::BRCA_PT,all = TRUE)
Matriz_Normal_Tumor <- Matriz_Normal_Tumor[,-c(1:3)]
Matriz_Normal_Tumor <- Matriz_Normal_Tumor[,-c(2:4)]

## En esta linea, hallamos el valor maximo para toda la tabla y sacamos el 0,1%
##para posteriormente ver si el gen esta expresado
Valor_Maximo <- max(Matriz_Normal_Tumor[,2:23688], na.rm=FALSE)*0.0002


### En las siguientes lineas de codigo vamos a binarizar la tabla, donde los valores
## que superen el umbral o valor maximo propuesto, se conviertan en 1 y los que no en 0


Matriz_Normal_Tumor_binario <- Matriz_Normal_Tumor
Matriz_Normal_Tumor_binario[, 2:23688] <- apply(Matriz_Normal_Tumor[, 2:23688], 2, function(x) ifelse(x > Valor_Maximo, 1, 0))

## En el siguiente paso, vamos a eliminar los genes que no estan expresados como minimo en el
## 20% de las observaciones

Porcentajes <- colMeans(Matriz_Normal_Tumor_binario[, 2:ncol(Matriz_Normal_Tumor_binario)] == 1)

Matriz_Normal_Tumor_binario_filtrado <- Matriz_Normal_Tumor_binario[, Porcentajes >= 0.2]

porcentaje_columna <- colMeans(Matriz_Normal_Tumor_binario_filtrado[, 2:ncol(Matriz_Normal_Tumor_binario_filtrado)]==1)
porcentaje_columna

###### A continuación renombramos las variables de acuerdo con el HGNC.symbol para por cruzarlas con el df PPI

PPIs_Renombrados <- AMCBGeneUtils::changeGeneId(PPI$`Input-node Gene Symbol`,from="HGNC.symbol")
colnames(Matriz_Normal_Tumor_binario_filtrado)[c(2:ncol(Matriz_Normal_Tumor_binario_filtrado))] <- PPIs_Renombrados$HGNC.symbol


####### A continuación contaremos el numero de veces que se repite cada gen, tanto en INPUT como en Output

Agrupacion <- PPI %>% group_by(`Input-node Gene Symbol`) %>%
  summarise(NN = n()) %>%
  rename(`Gen`=`Input-node Gene Symbol`)

Agrupacion2 <- PPI %>% group_by(`Output-node Gene Symbol`) %>%
  summarise(NN = n()) %>%
  rename(`Gen`=`Output-node Gene Symbol`)

### Combinamos las dos tablas que hacen el conteo tanto de Input como de Output
tabla_combinada <- bind_rows(Agrupacion, Agrupacion2) %>%
  group_by(`Gen`) %>%
  summarise(NN = sum(NN)) %>%
  arrange(desc(NN)) %>%
  top_n(101, NN)

Cambiar_nombres_tabla_combinada <- AMCBGeneUtils::changeGeneId(tabla_combinada$Gen,from="HGNC.symbol")
rownames(tabla_combinada) <- Cambiar_nombres_tabla_combinada$HGNC.symbol

nombres_variables <- colnames(Matriz_Normal_Tumor_binario_filtrado)
# Convertir los nombres de las variables en una matriz
matriz_nombres_variables <- as.matrix(nombres_variables)

valores_comunes <- intersect(tabla_combinada$Gen, matriz_nombres_variables)
valores_comunes <- as.matrix(valores_comunes)

##### a Continuación la matriz final con la que trabajaremos los modelos
Matriz_Final <- Matriz_Normal_Tumor_binario_filtrado[,c("sample_type",valores_comunes)]
Matriz_Final <- Matriz_Final[, !names(Matriz_Final) == "APP"]


##################################################################################
#A partir de aca se implementaran los modelos#####################
##################################################################################

######################## KNN #######################################################

MatrizFinal_df <- Matriz_Final
MatrizFinal_df$sample_type <- as.factor(MatrizFinal_df$sample_type)

set.seed(1)
sample.index <- sample(1:nrow(MatrizFinal_df)
                       ,nrow(MatrizFinal_df)*0.7
                       ,replace = F)

predictor <- colnames(MatrizFinal_df)[-1]

train.data <- MatrizFinal_df[sample.index,c(predictor,"sample_type"),drop=F]
test.data <- MatrizFinal_df[-sample.index,c(predictor,"sample_type"),drop=F]

ctrl <- trainControl(method="cv", p = 0.7)
knnTrain <- train( sample_type ~ .
                   , data = train.data
                   , method = "knn", trControl = ctrl
                   , preProcess = c("range") #c("center", "scale")
                   , tuneLength = 25)

knnTrain
plot(knnTrain)

knnPredict <- predict(knnTrain, newdata = test.data)
knnPredict

confusionMatrix(knnPredict, test.data$sample_type)

####################### REGRESIÓN LINEAL #######################################
ins_model_1 <- lm(as.numeric(sample_type) ~ ., data = train.data)
ins_model_1

predict(ins_model_1, newdata = test.data)
summary(ins_model_1)

excl_v <- c("TP53", "CREBBP", "YWHAG", "SMAD3", "GRB2", "SRC", "SMAD2", "CDKN1A",
            "FYN", "TK1", "SMAD4", "JUN", "CCDC85B", "MAPK6", "GSK3B","PIK3R1",
            "SHC1", "TRAF2", "YWHAZ", "CASP3", "UBE2I", "SP1", "VIM", "ATXN1",
            "SMN1", "UBQLN4", "PRKACA", "TGFBR1", "CALM1", "SETDB1", "BRCA1",
            "CTNNB1", "LCK", "RXRA", "EEF1A1", "AKT1", "STAT3", "PTPN11", "NCOA1",
            "PLCG1", "ACTB", "MDFI",  "EWSR1", "RAC1", "NFKB1", "NR3C1", "UNC119",
            "ABL1", "DLG4", "ATN1","NCOR2", "CDK2", "CHD3", "PRKCD", "MAPK14", "TLE1",
            "XRCC6", "CBL", "INSR", "PTN", "ZBTB16", "KAT5", "CAV1", "RAF1", "STAT1",
            "COPS6", "KAT2B", "PTPN6", "MAPK8", "PXN", "ACTA1", "NCOR1", "PIN1", "TRAF6",
            "ANXA7", "LYN"
)

train.data_Rl <- train.data[, !(names(train.data) %in% excl_v)]
test.data_Rl <- test.data[, !(names(test.data) %in% excl_v)]

model_excl_Rl <- lm(as.numeric(sample_type) ~ ., data = train.data_Rl)
summary(model_excl_Rl)

###################### REGRESIÓN LOGISTICA #####################################

MF_df <- Matriz_Final
MF_df$sample_type <- MF_df$sample_type=="Primary Tumor"
MF_df$sample_type <- 1*MF_df$sample_type
MF_df$sample_type <- as.factor(MF_df$sample_type)

sample.index2 <- sample(1:nrow(MF_df)
                        ,nrow(MF_df)*0.7
                        ,replace = F)

predictor2 <- colnames(MF_df)[-1]

train.data2 <- MF_df[sample.index2,c(predictor2,"sample_type"),drop=F]
test.data2 <- MF_df[-sample.index2,c(predictor2,"sample_type"),drop=F]

LogisticTrain <- glm(formula = sample_type ~ .
                     ,data = train.data2
                     ,family = "binomial")
LogisticTrain

Lpredict <- predict(LogisticTrain, newdata = test.data2
                    , type = "response")
Lpredict

predicted_classes <- ifelse(Lpredict > 0.5, 1, 0)
predicted_classes <- factor(predicted_classes, levels = levels(test.data2$sample_type))

confusion_matrix <- confusionMatrix(data = predicted_classes, reference = test.data2$sample_type)
print(confusion_matrix)

################### ÁRBOL DE DECISION ##########################################


head (Matriz_Final)

Matriz_Final[, "sample_type"] <- as.factor(Matriz_Final[, "sample_type"])

modelo_arbol <- rpart(sample_type ~ ., data = Matriz_Final, na.action = na.omit)

plot(modelo_arbol, uniform=T, margin=0.1)
text(modelo_arbol, use.n=T, all=T, cex=0.8)

rpart.plot::rpart.plot(modelo_arbol, tweak=1.5)


modelo_arbol$cptable

Etiquetas_Modelo_Arbol <- predict(modelo_arbol, Matriz_Final[,-1], type="class")
confusionMatrix(Etiquetas_Modelo_Arbol, Matriz_Final[,1])



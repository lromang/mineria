################################################
## Luis Manuel Román García
## Minería de Datos
## Proyecto 1
################################################

################################################
## Librerías utilizadas
################################################
library(caret)
library(neuralnet)
library(plyr)
library(GGally)
################################################
## Descripción de variables
################################################
##1. Sequence Name: Accession number for the SWISS-PROT database
##2. mcg: McGeoch's method for signal sequence recognition.
##3. gvh: von Heijne's method for signal sequence recognition.
##4. alm: Score of the ALOM membrane spanning region prediction program.
##5. mit: Score of discriminant analysis of the amino acid content
##of the N-terminal region (20 residues long) of mitochondrial
##and non-mitochondrial proteins.
##6. erl: Presence of "HDEL" substring (thought to act as a signal
##for retention in the endoplasmic reticulum lumen). Binary attribute.
##7. pox: Peroxisomal targeting signal in the C-terminus.
##8. vac: Score of discriminant analysis of the amino acid content
##of vacuolar and extracellular proteins.
##9. nuc: Score of discriminant analysis of nuclear localization signals of nuclear and non-nuclear proteins.

################################################
## Lectura de datos
################################################
data <- read.table("./data/yeast.data",
                stringsAsFactors = FALSE,
                header = FALSE
                )
## La primera columna continene los nombres de registro
## en la base de datos por tanto, no tiene ningún significado
## para la predicción.
data <- data[,-1]

################################################
## Dado que la variable 6 es binaria, la cambiaremos
## a tipo factor.
data[,5] <- as.factor(data[,5])

################################################
## Segmentar datos en entrenamiento y prueba
################################################
## Llevaremos a cabo la participación de los datos
## con base a la variable que queremos predecir de
## tal forma que se preserven las mismas proporciones
## en ambas muestras. Utilizaremos el 70% de los datos
## para entrenar y el 30% para calcular nuestro error
## de generalización. Notemos que las clases no están
## balanceadas.
## 1  CYT  463
## 2  ERL    5
## 3  EXC   35
## 4  ME1   44
## 5  ME2   51
## 6  ME3  163
## 7  MIT  244
## 8  NUC  429
## 9  POX   20
## 10 VAC   30
set.seed(117077)
sample <- createDataPartition(data$V10,
                             p = .7,
                             list = FALSE,
                             times = 1)
data.train <- data[sample,]
data.test  <- data[-sample,]

################################################
## Preprocesamiento de variables
################################################

################################################
## Variables con sesgo. 
## Veamos si la variable 6 está sesgada
nearZeroVar(data.train, saveMetrics = TRUE)
## Vemos que tanto la variable 6 como la variable 7
## presentan una distribución sesgada. 
count(data.train[,c(5,9)])

################################################
## Antes que nada, es pertinente escalar los datos
## para asegurar que los resultados no se ven afectados
## por la dimensionalidad de estos.
params    <- preProcess(data.train[,-c(5,9)],
                       method = c("center",
                                  "scale"))

data.train[,-c(5,9)] <- predict(params, data.train[,-c(5,9)])
data.test[,-c(5,9)]  <- predict(params, data.test[,-c(5,9)])

################################################
## Identificar predictores que estén altamente
## correlacionados.
cor_plot <- ggpairs(data=data.train,
                   columns = 1:8,
                   upper = list(continuous = "blank"),
                   lower = list(continuous = "points",
                                combo = "dot"),
                   colour = "V10")
## Las únicas variables que presentan correlación
## lineal son V2 y V3.
cor_mat <- cor(data.train[,-c(5,9)])

## Sin embargo, esta correlación es menor a .6
alta_cor <- findCorrelation(cor_mat, cutoff = .3)

## De hecho son las únicas variables que tienen
## una correlación mayor a abs(.3).


################################################
## Entrenamiento de modelo
################################################
class <- data.train[,9]
pred  <- data.train[,-9]
################################################
## Validación cruzada
control <- trainControl(method  = "repeatedcv",
                       number  = 10,
                       repeats = 10)
################################################
model1  <- train(pred, class,
                method = "nnet",
                trControl = control,
                verbose = FALSE)

################################################
## Error de entrenamiento
################################################
## Predecimos valores de entrenamiento
preds <- predict(model1, newdata = data.train[,-9])

## Error de entrenamiento
err_train <- (sum(preds == data.train[,9])/nrow(data.train))

## Benchmark
data_bench <- count(class)
bench_pred <- sample(data_bench$x, nrow(data.train), replace = TRUE,
                    prob = data_bench$freq/nrow(data.train))

## Error Benchmark
err_bench <- (sum(bench_pred == data.train[,9])/nrow(data.train))

## Vemos una clara mejora
err_train/err_bench


################################################
## Error de prueba
################################################
## Predecimos valores de entrenamiento
preds <- predict(model1, newdata = data.test[,-9])

## Error de entrenamiento
err_train <- (sum(preds == data.test[,9])/nrow(data.test))

## Benchmark
data_bench <- count(class)
bench_pred <- sample(data_bench$x, nrow(data.test), replace = TRUE,
                    prob = data_bench$freq/nrow(data.train))

## Error Benchmark
err_bench <- (sum(bench_pred == data.test[,9])/nrow(data.test))

## Vemos una clara mejora
err_test/err_bench

########################################################################
# TB_assinatura_abordagemA.R
#
# Identificação de assinatura transcriptômica em sangue total para
# diferenciação entre tuberculose ativa e infecção latente: validação
# em múltiplas coortes
#
# Written by   : Gustavo Pires Ramos
# Contexto     : TCC - Bacharelado em Biotecnologia, ICS/UFBA
# Orientação   : Dra. Mariana Araújo-Pereira (IGM - Fiocruz Bahia)
# Coorientação : Isabella Brige Ferreira
# Data         : novembro de 2026 (última revisão: 25/11/2026)
# Ambiente     : R 4.5.2
#
# Dados de entrada (pasta R/meta/), um .xlsx por série do GEO:
#   GSE73408, GSE39940, GSE39939, GSE37250, GSE19439, GSE19442, GSE19444
#   aba de expressão + aba 2 com os metadados da série
#
# Etapas
#   1. Pacotes
#   2. Leitura e filtragem das séries
#   3. Integração das séries e ComBat (ajustado no treino)
#   4. PCA antes e depois do ComBat
#   5. Expressão diferencial (limma/TREAT) e enriquecimento
#   6. Seleção da assinatura (RFE + Random Forest, CV por estudo)
#   7. Avaliação no conjunto de teste (30%)
#   -  Teste complementar: estabilidade do split entre seeds
#   -  Teste complementar: validação por exclusão de estudo (LODO)
#   8. Análise de coexpressão (CEMiTool)
#   9. Figuras do TCC e painel ROC por série
#
# Saídas principais
#   combat_one_ok2.xlsx, metadado_merged_ok2.xlsx
#   validacao_leave_one_dataset_out.xlsx
#   cemitool_resultados.xlsx, cemitool_report/
#   R/meta/figuras_TCC_A/figA_*.png
#
# Notas
#   - seed 42 no pipeline principal. O split 70/30 é refeito
#     com a mesma seed no ComBat, nos DEGs e no RFE, então as três
#     etapas usam exatamente as mesmas amostras de treino.
#   - alguns comentários antigos citam o FCGR1A como assinatura final;
#     eles são de uma versão anterior. A assinatura usada é
#     assinatura_1 = "AIM2" (definida na etapa 6).
########################################################################

# install.packages(c("readxl", "readr", "dplyr"))
# 
# # Instalar pacotes do Bioconductor
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c(
#   "TBSignatureProfiler", 
#   "GEOquery", 
#   "affy", 
#   "CEMiTool"
# ))

########################################################################
# 1. PACOTES
########################################################################

#### Load Packages ####
{
  library(GEOquery)
  library(readxl)
  library(affy)
  library(readr)
  library(dplyr)
  library(openxlsx)
  library(stringr)
  library(h2o)
  library(GA)
  library(caret)
  library(randomForest)
  
}


########################################################################
# 2. LEITURA E FILTRAGEM DAS SÉRIES
#
# Para cada série: lê a expressão e os metadados, remove HIV+,
# outras doenças, TB subclínica, pneumonia e controles saudáveis,
# e mantém só as amostras presentes nas duas tabelas (os outliers
# dos estudos originais já não estão na matriz depositada).
# Abas de expressão: GSE73408 "collapsed_data" (RMA);
# GSE19439/42/44 "log2_data"; demais séries, aba 1 (já em log2).
########################################################################

#### Load the data ####
{
  # GSE37250 Filter
  # Filter HIV positive condition
  GSE37250_All_data <- openxlsx::read.xlsx("R/meta/GSE37250_All_data.xlsx", sheet = 1, rowNames = T) #Dataset already cleaned for OD and HIV positive
  boxplot(GSE37250_All_data, las = 2, outline = F)
  metadado_GSE37250 <- openxlsx::read.xlsx("R/meta/GSE37250_All_data.xlsx", sheet = 2, rowNames = T)
  metadado2_GSE37250 <-openxlsx::read.xlsx("R/meta/GSE37250_All_data.xlsx", sheet = 2, rowNames = T)
  metadado_GSE37250$title <- gsub(" ", ".", metadado_GSE37250$title, fixed = T)
  metadado_GSE37250 <- metadado_GSE37250[!metadado_GSE37250$`hiv.status:ch1` %in% "HIV positive",]
  metadado_GSE37250 <- metadado_GSE37250 %>% 
    filter(!`disease.state:ch1` %in% "other disease")
  #Filter outliers
  common_samples <- intersect(colnames(GSE37250_All_data), metadado_GSE37250$title)
  GSE37250_All_data <- GSE37250_All_data[, colnames(GSE37250_All_data) %in% common_samples]
  metadado_GSE37250 <- metadado_GSE37250[metadado_GSE37250$title %in% common_samples, , drop = FALSE]
  dim(metadado_GSE37250)
  dim(GSE37250_All_data)
  
  # GSE73408 Filter
  # Filter PNA and OD condition
  GSE73408_All_data <- openxlsx::read.xlsx("R/meta/GSE73408_All_data.xlsx", sheet = "collapsed_data", rowNames = T) #Dataset already cleaned for PNA,log aplicado duas vezes na aba 1 do GSE73408
  boxplot(GSE73408_All_data, las = 2, outline = F)
  #summary(as.vector(as.matrix(GSE73408_All_data))) # esperado: ~0,8 a 13,2
  metadado_GSE73408 <- openxlsx::read.xlsx("R/meta/GSE73408_All_data.xlsx", sheet = 2, rowNames = T)
  metadado2_GSE73408 <- openxlsx::read.xlsx("R/meta/GSE73408_All_data.xlsx", sheet = 2, rowNames = T)
  row.names(GSE73408_All_data)
  # Formatar rownames
  # Extrair a primeira parte antes de " // "
  # new_rownames <- as.data.frame(str_split(rownames(GSE73408_All_data), " // ", n = Inf, simplify = TRUE))
  # head(new_rownames)
  # new_rownames <- str_split(rownames(GSE73408_All_data), " // ", n = 3, simplify = TRUE)[, 2]
  
  # Verificar e corrigir duplicatas
  # if (any(duplicated(new_rownames))) {
  #   new_rownames <- make.unique(new_rownames, sep = "_")  # Adiciona _1, _2, etc.
  # }
  # rownames(GSE73408_All_data) <- new_rownames
  
  # Verificar se está correto
  # sum(duplicated(rownames(GSE73408_All_data)))  #verifica se retorna 0, se sim tem genes duplicados
  # colnames(GSE73408_All_data) %in% metadado_GSE73408$title
  
  # Filter outliers
  
  GSE73408_All_data <- GSE73408_All_data[, colnames(GSE73408_All_data) %in% metadado_GSE73408$title]
  metadado_GSE73408 <- metadado_GSE73408[metadado_GSE73408$title %in% colnames(GSE73408_All_data),]
  metadado_GSE73408 <- metadado_GSE73408[order(match(metadado_GSE73408$title, colnames(GSE73408_All_data))),]
  metadado_GSE73408$title == colnames(GSE73408_All_data)
  dim(GSE73408_All_data)
  dim(metadado_GSE73408)
  
  # GSE39940 Filter
  # Filter HIV positive condition and OD condition
  GSE39940_All_data <- openxlsx::read.xlsx("R/meta/GSE39940_All_data.xlsx", sheet = 1, rowNames = T) 
  boxplot(GSE39940_All_data, las = 2, outline = F)
  metadado_GSE39940 <- openxlsx::read.xlsx("R/meta/GSE39940_All_data.xlsx", sheet = 2, rowNames = T) #Filter ´TB status:subclinical´
  metadado2_GSE39940 <- openxlsx::read.xlsx("R/meta/GSE39940_All_data.xlsx", sheet = 2, rowNames = T)
  metadado_GSE39940 <- metadado_GSE39940[
    !(
      metadado_GSE39940$`hiv.status:ch1` %in% "HIV positive" |
        metadado_GSE39940$`disease.status:ch1` %in% "other disease"
    ),
  ]
  #Filter outliers
  common_samples <- intersect(colnames(GSE39940_All_data), metadado_GSE39940$title)
  GSE39940_All_data <- GSE39940_All_data[, colnames(GSE39940_All_data) %in% common_samples]
  metadado_GSE39940 <- metadado_GSE39940[metadado_GSE39940$title %in% common_samples, , drop = FALSE]
  dim(metadado_GSE39940)
  dim(GSE39940_All_data)
  
  dim(metadado_GSE39940)
  dim(GSE39940_All_data)
  
  # GSE39939 Filter
  # Filter other disease condition *criterio para confirmar a cultura TB
  GSE39939_All_data <- openxlsx::read.xlsx("R/meta/GSE39939_All_data.xlsx", sheet = 1, rowNames = T)  
  boxplot(GSE39939_All_data, las = 2, outline = F)
  metadado_GSE39939 <-openxlsx::read.xlsx("R/meta/GSE39939_All_data.xlsx", sheet = 2, rowNames = T) 
  metadado2_GSE39939 <-openxlsx::read.xlsx("R/meta/GSE39939_All_data.xlsx", sheet = 2, rowNames = T)
  metadado_GSE39939 <- metadado_GSE39939[!metadado_GSE39939$`illness:ch1` %in% c("other disease", "other disease (IGRA +)"),]
  metadado_GSE39939 <- metadado_GSE39939[!metadado_GSE39939$`hiv.status:ch1` %in% "HIV positive",]
  ##
  #Filter outliers
  common_samples <- intersect(colnames(GSE39939_All_data), metadado_GSE39939$title)
  GSE39939_All_data <- GSE39939_All_data[, colnames(GSE39939_All_data) %in% common_samples]
  metadado_GSE39939 <- metadado_GSE39939[metadado_GSE39939$title %in% common_samples, , drop = FALSE]
  # Replace the right condition on the metadada and data 
  metadado_GSE39939$title <- mapply(
    function(title, status) {
      gsub("^WB_\\d+(.*)$", paste0(status, "\\1"), title)
    },
    metadado_GSE39939$title,metadado_GSE39939$`illness:ch1`
  )
  
  colnames(GSE39939_All_data) <- mapply(
    function(title, status) {
      gsub("^WB_\\d+(.*)$", paste0(status, "\\1"), title)
    },
    colnames(GSE39939_All_data), metadado_GSE39939$title
  )
  
  colnames(GSE39939_All_data) <- sub("_(KIL_[0-9]+)$", "", colnames(GSE39939_All_data))
  
  dim(metadado_GSE39939)
  dim(GSE39939_All_data)
  
  # GSE19439 Filter
  # Filter healthy control
  GSE19439_All_data <- openxlsx::read.xlsx("R/meta/GSE19439_All_data.xlsx", sheet = "log2_data", rowNames = TRUE) # aba 3 estava em escala linear, não em log2 (o que era mais próximo)
  boxplot(GSE19439_All_data, las = 2, outline = F)
  metadado_GSE19439 <- openxlsx::read.xlsx("R/meta/GSE19439_All_data.xlsx", sheet = 2, rowNames = T) 
  metadado2_GSE19439 <- openxlsx::read.xlsx("R/meta/GSE19439_All_data.xlsx", sheet = 2, rowNames = T)
  metadado_GSE19439 <- metadado_GSE19439[!grepl("^CON_LON_", metadado_GSE19439$title),]
  GSE19439_All_data <- GSE19439_All_data[!grepl("^CON_LON_", colnames(GSE19439_All_data))]
  
  # # Rename expression dataset for condition; location; pacient_id
  
  #Filter Outliers
  common_samples <- intersect(colnames(GSE19439_All_data), metadado_GSE19439$title)
  GSE19439_All_data <- GSE19439_All_data[, colnames(GSE19439_All_data) %in% common_samples]
  metadado_GSE19439 <- metadado_GSE19439[metadado_GSE19439$title %in% common_samples, , drop = FALSE]
  dim(metadado_GSE19439)
  dim(GSE19439_All_data)
  
  # GSE19442 Filter
  # No Filter
  GSE19442_All_data <- openxlsx::read.xlsx("R/meta/GSE19442_All_data.xlsx", sheet = "log2_data", rowNames = T) #aba 1 é a correta com log2, diminui a quantidade de negativos ao invés da aba 3 - collapsed data
  boxplot(GSE19442_All_data, las = 2, outline = F)
  metadado_GSE19442 <- openxlsx::read.xlsx("R/meta/GSE19442_All_data.xlsx", sheet = 2, rowNames = T)
  metadado2_GSE19442 <- openxlsx::read.xlsx("R/meta/GSE19442_All_data.xlsx", sheet = 2, rowNames = T)
  
  common_samples <- intersect(colnames(GSE19442_All_data), metadado_GSE19442$title)
  GSE19442_All_data <- GSE19442_All_data[, colnames(GSE19442_All_data) %in% common_samples]
  metadado_GSE19442 <- metadado_GSE19442[metadado_GSE19442$title %in% common_samples, , drop = FALSE]
  dim(metadado_GSE19442)
  dim(GSE19442_All_data)
  
  
  # GSE19444
  # Filter healthy control
  # All condition included - PTB, LTB, LATENT
  GSE19444_All_data <- openxlsx::read.xlsx("R/meta/GSE19444_All_data.xlsx", sheet = "log2_data", rowNames = T)
  boxplot(GSE19444_All_data, las = 2, outline = F)
  metadado_GSE19444 <- openxlsx::read.xlsx("R/meta/GSE19444_All_data.xlsx", sheet = 2, rowNames = T)
  metadado2_GSE19444 <- openxlsx::read.xlsx("R/meta/GSE19444_All_data.xlsx", sheet = 2, rowNames = T)
  #Excluded HC - Healthy Control
  GSE19444_All_data <- GSE19444_All_data[!grepl("^CON_LON_", colnames(GSE19444_All_data))]
  metadado_GSE19444 <- metadado_GSE19444[!grepl("^CON_LON_", metadado_GSE19444$title),]
  
  # Filter outliers
  common_samples <- intersect(colnames(GSE19444_All_data), metadado_GSE19444$title)
  GSE19444_All_data <- GSE19444_All_data[, colnames(GSE19444_All_data) %in% common_samples]
  metadado_GSE19444 <- metadado_GSE19444[metadado_GSE19444$title %in% common_samples, , drop = FALSE]
  dim(metadado_GSE19444)
  dim(GSE19444_All_data)
}



########################################################################
# 3. INTEGRAÇÃO DAS SÉRIES E CORREÇÃO DE LOTE (ComBat)
#
# - junta as 7 matrizes pelos genes em comum e padroniza os rótulos
#   (Active / Latent tuberculosis); batch = série (1 a 7)
# - split 70/30 estratificado pelo grupo clínico ANTES da correção
# - ComBat paramétrico (Johnson et al., 2007): parâmetros estimados
#   no treino e aplicados ao teste, sem os rótulos clínicos
# Saída: combat_df (genes x amostras), metadado_merged, batch
########################################################################

#### SVA and merge #### 
{
  library(dplyr)
  library(plyr)
  library(sva)
  library(openxlsx)
  
  # Ordenar a lista de acordo com a sequência desejada
  ordered_names <- c("GSE73408","GSE39940","GSE39939","GSE37250","GSE19439","GSE19442","GSE19444")
  
  # Adicionar coluna X com os nomes das linhas para cada dataset
  GSE73408_All_data$X <- row.names(GSE73408_All_data)
  GSE39940_All_data$X <- row.names(GSE39940_All_data)
  GSE39939_All_data$X <- row.names(GSE39939_All_data)
  GSE37250_All_data$X <- row.names(GSE37250_All_data)
  GSE19439_All_data$X <- row.names(GSE19439_All_data)
  GSE19442_All_data$X <- row.names(GSE19442_All_data)
  GSE19444_All_data$X <- row.names(GSE19444_All_data)
  
  expression_matrix <- na.omit(plyr::join_all(list(GSE73408_All_data, GSE39940_All_data, GSE39939_All_data,
                                                   GSE37250_All_data, GSE19439_All_data, GSE19442_All_data,
                                                   GSE19444_All_data), by = "X", type = "full"))
  
  
  # Criar lista ordenada dos datasets
  data_expression_ordered <- list(
    GSE73408_All_data,
    GSE39940_All_data,
    GSE39939_All_data,
    GSE37250_All_data,
    GSE19439_All_data,
    GSE19442_All_data,
    GSE19444_All_data
  )
  
  # Fazer full join de todos os datasets pela coluna "X"
  data_expression_all <- Reduce(function(x, y) dplyr::full_join(x, y, by = "X"), data_expression_ordered)
  
  # Remover coluna "X" para ficar só com dados numéricos
  rownames(data_expression_all) <- data_expression_all$X
  data_expression_all$X <- NULL
  
  #
  metadado_GSE73408$`clinical.group:ch1` <- gsub("^TB$", "Active tuberculosis", metadado_GSE73408$`clinical.group:ch1`)
  metadado_GSE73408$`clinical.group:ch1` <- gsub("^LTBI$", "Latent tuberculosis", metadado_GSE73408$`clinical.group:ch1`)
  metadado_GSE73408$`clinical.group:ch1` <- gsub("^PNA$", "Other", metadado_GSE73408$`clinical.group:ch1`)
  metadado_GSE73408 <- metadado_GSE73408[metadado_GSE73408$`clinical.group:ch1` != "Other", ]
  #
  metadado_GSE39940$`disease.status:ch1` <- gsub("^active tuberculosis", "Active tuberculosis",metadado_GSE39940$`disease.status:ch1`)
  metadado_GSE39940$`disease.status:ch1` <- gsub("^latent TB infection", "Latent tuberculosis",metadado_GSE39940$`disease.status:ch1`)
  #
  metadado_GSE39939$`source_name_ch1` <- gsub ("^Whole Blood from patient with Active TB", "Active tuberculosis", metadado_GSE39939$`source_name_ch1`)
  metadado_GSE39939$`source_name_ch1` <- gsub ("^Whole Blood from patient with LTBI", "Latent tuberculosis", metadado_GSE39939$`source_name_ch1`)
  #
  metadado_GSE37250$`source_name_ch1` <- gsub ("^Whole Blood from patient with Active TB", "Active tuberculosis", metadado_GSE37250$`source_name_ch1`)
  metadado_GSE37250$`source_name_ch1` <- gsub ("^Whole Blood from patient with LTBI", "Latent tuberculosis", metadado_GSE37250$`source_name_ch1`)
  #
  metadado_GSE19439$`source_name_ch1` <- gsub ("^Whole Blood from patient with Active TB", "Active tuberculosis", metadado_GSE19439$`source_name_ch1`)
  metadado_GSE19439$`source_name_ch1` <- gsub ("^Whole Blood from patient with Latent TB", "Latent tuberculosis", metadado_GSE19439$`source_name_ch1`)
  #
  metadado_GSE19442$`source_name_ch1` <- gsub ("^Whole Blood from patient with Active TB", "Active tuberculosis", metadado_GSE19442$`source_name_ch1`)
  metadado_GSE19442$`source_name_ch1` <- gsub ("^Whole Blood from patient with Latent TB", "Latent tuberculosis", metadado_GSE19442$`source_name_ch1`)
  #
  metadado_GSE19444$`source_name_ch1` <- gsub ("^Whole Blood from patient with Active TB", "Active tuberculosis", metadado_GSE19444$`source_name_ch1`)
  metadado_GSE19444$`source_name_ch1` <- gsub ("^Whole Blood from patient with Latent TB", "Latent tuberculosis", metadado_GSE19444$`source_name_ch1`)
  
  # Criar metadados para cada dataset com batch e grupos
  metadado2_GSE73408 <- data.frame(row.names = metadado_GSE73408$title, 
                                   sample = 1:nrow(metadado_GSE73408),
                                   batch = "1",
                                   Groups = metadado_GSE73408$`clinical.group:ch1`)
  
  
  metadado2_GSE39940 <- data.frame(row.names = metadado_GSE39940$title, 
                                   sample = 1:nrow(metadado_GSE39940),
                                   batch = "2",
                                   Groups = metadado_GSE39940$`disease.status:ch1`)
  
  metadado2_GSE39939 <- data.frame(row.names = metadado_GSE39939$title,
                                   sample = 1:nrow(metadado_GSE39939),
                                   batch = "3",
                                   Groups = metadado_GSE39939$`source_name_ch1`)
  
  metadado2_GSE37250 <- data.frame(row.names = metadado_GSE37250$title,
                                   sample = 1:nrow(metadado_GSE37250),
                                   batch = "4",
                                   Groups = metadado_GSE37250$`source_name_ch1`)
  
  metadado2_GSE19439 <- data.frame(row.names = metadado_GSE19439$title, 
                                   sample = 1:nrow(metadado_GSE19439),
                                   batch = "5",
                                   Groups = metadado_GSE19439$`source_name_ch1`)
  
  metadado2_GSE19442 <- data.frame(row.names = metadado_GSE19442$title,
                                   sample = 1:nrow(metadado_GSE19442),
                                   batch = "6",
                                   Groups = metadado_GSE19442$`source_name_ch1`)
  
  metadado2_GSE19444 <- data.frame(row.names = metadado_GSE19444$title,
                                   sample = 1:nrow(metadado_GSE19444),
                                   batch = "7",
                                   Groups = metadado_GSE19444$`source_name_ch1`)
  
  # Combinar todos os metadados
  metadado_merged <- rbind(metadado2_GSE73408, metadado2_GSE39940, metadado2_GSE39939, 
                           metadado2_GSE37250, metadado2_GSE19439, metadado2_GSE19442, metadado2_GSE19444)
  
  
  # Ajustar a ordem do metadado para coincidir com as colunas da matriz de expressão
  #metadado_merged <- metadado_merged[colnames(data_expression_all), ]
  
  common_expression <- intersect(colnames(data_expression_all), rownames(metadado_merged))
  
  # Filtrar expressão e metadado para manter só amostras em comum
  expression_matrix <- data_expression_all[, common_expression]
  metadado_merged   <- metadado_merged[common_expression ,]
  
  print(sum(is.na(metadado_merged$batch)))
  print(sum(is.na(metadado_merged$Groups)))
  
  # Extrair vetor batch e fator Groups
  batch <- as.factor(metadado_merged$batch)
  groups <- factor(metadado_merged$Groups)
  
  # Converter matriz de expressão para numérica (se necessário)
  # expression_matrix <- as.matrix(data_expression_all)
  # mode(expression_matrix) <- "numeric"
  
  metadado_merged <- na.omit(metadado_merged)
  expression_matrix <- as.data.frame(na.omit(expression_matrix))
  
  ### CORREÇÃO DE VAZAMENTO ###
  # alinhamento expressão × metadado depois dos na.omit()
  common_pos_combat <- intersect(colnames(expression_matrix), rownames(metadado_merged))
  expression_matrix <- expression_matrix[, common_pos_combat, drop = FALSE]
  metadado_merged   <- metadado_merged[common_pos_combat, , drop = FALSE]
  stopifnot(identical(colnames(expression_matrix), rownames(metadado_merged)))
  batch  <- as.factor(metadado_merged$batch)
  groups <- factor(metadado_merged$Groups)
  ### FIM DA CORREÇÃO ###
  
  dim(expression_matrix)  # genes x amostras
  length(batch)           # deve ser = ncol(expression_matrix)
  length(groups)          # deve ser = ncol(expression_matrix)
  
  # conferência de escala: todos os estudos em log2 (~0 a 16)
  print(t(sapply(split(colnames(expression_matrix), batch), function(s)
    summary(as.vector(as.matrix(expression_matrix[, s]))))))
  
  head(expression_matrix)
  
  ### CORREÇÃO DE VAZAMENTO ###
  # split treino/teste ANTES da correção de lote; ComBat ajustado só no treino
  # e aplicado ao teste, sem mod (não usa os rótulos)
  
  y_combat_split <- factor(metadado_merged$Groups)
  levels(y_combat_split) <- make.names(levels(y_combat_split))
  
  set.seed(42)
  idx_treino_combat <- createDataPartition(y_combat_split, p = 0.7, list = FALSE)
  
  expression_treino_bruto <- expression_matrix[, idx_treino_combat, drop = FALSE]
  expression_teste_bruto  <- expression_matrix[, -idx_treino_combat, drop = FALSE]
  batch_treino_combat <- droplevels(batch[idx_treino_combat])
  batch_teste_combat  <- droplevels(batch[-idx_treino_combat])
  
  # --- ComBat paramétrico (Johnson et al., 2007), treino → teste ---
  combat_treinar <- function(dat, batch, conv = 1e-4, max_iter = 1000) {
    dat <- as.matrix(dat); batch <- as.character(batch)
    lotes <- sort(unique(batch))
    n_b <- as.numeric(table(batch)[lotes]); n <- ncol(dat)
    
    media_b  <- sapply(lotes, function(b) rowMeans(dat[, batch == b, drop = FALSE]))
    grand    <- as.vector(media_b %*% (n_b / n))
    resid    <- dat - media_b[, match(batch, lotes)]
    var_pool <- rowSums(resid^2) / n
    
    s <- (dat - grand) / sqrt(var_pool)
    
    gamma_hat <- sapply(lotes, function(b) rowMeans(s[, batch == b, drop = FALSE]))
    delta_hat <- sapply(lotes, function(b) apply(s[, batch == b, drop = FALSE], 1, var))
    gamma_star <- gamma_hat; delta_star <- delta_hat
    
    for (b in lotes) {
      sb <- s[, batch == b, drop = FALSE]; nb <- ncol(sb)
      g_bar <- mean(gamma_hat[, b]); t2 <- var(gamma_hat[, b])
      m <- mean(delta_hat[, b]);     s2 <- var(delta_hat[, b])
      a  <- (2 * s2 + m^2) / s2
      bb <- (m * s2 + m^3) / s2
      g_old <- gamma_hat[, b]; d_old <- delta_hat[, b]
      for (it in seq_len(max_iter)) {
        g_new <- (t2 * nb * gamma_hat[, b] + d_old * g_bar) / (t2 * nb + d_old)
        d_new <- (0.5 * rowSums((sb - g_new)^2) + bb) / (nb / 2 + a - 1)
        ch <- max(abs(g_new - g_old) / (abs(g_old) + 1e-12), abs(d_new - d_old) / d_old)
        g_old <- g_new; d_old <- d_new
        if (ch < conv) break
      }
      gamma_star[, b] <- g_new; delta_star[, b] <- d_new
    }
    
    list(grand = grand, var_pool = var_pool, gamma = gamma_star,
         delta = delta_star, genes = rownames(dat))
  }
  
  # genes sem variância em algum estudo do treino são removidos (o ComBat falha com eles)
  sd_ok <- sapply(split(seq_along(batch_treino_combat), batch_treino_combat), function(i)
    apply(as.matrix(expression_treino_bruto[, i]), 1, sd) > 0)
  genes_ok <- rownames(expression_treino_bruto)[rowSums(!sd_ok) == 0]
  
  # parâmetros aprendidos no treino (mesmo nome do objeto original)
  parametros_por_estudo <- combat_treinar(expression_treino_bruto[genes_ok, ], batch_treino_combat)
  
  # mesma assinatura da função original: aplicar_correcao_batch(X, batch_vec)
  aplicar_correcao_batch <- function(X, batch_vec) {
    fit <- parametros_por_estudo
    X <- as.matrix(X)[fit$genes, , drop = FALSE]; batch_vec <- as.character(batch_vec)
    if (!all(batch_vec %in% colnames(fit$gamma)))
      stop("Há lotes sem parâmetros aprendidos no treino.")
    s <- (X - fit$grand) / sqrt(fit$var_pool)
    (s - fit$gamma[, batch_vec, drop = FALSE]) / sqrt(fit$delta[, batch_vec, drop = FALSE]) *
      sqrt(fit$var_pool) + fit$grand
  }
  
  combat_df_treino <- aplicar_correcao_batch(expression_treino_bruto, batch_treino_combat)
  combat_df_teste  <- aplicar_correcao_batch(expression_teste_bruto,  batch_teste_combat)
  
  combat_df <- as.data.frame(cbind(combat_df_treino, combat_df_teste))
  combat_df <- combat_df[, colnames(expression_matrix), drop = FALSE]
  
  cat("OK: ComBat ajustado sem vazamento (treino:", ncol(combat_df_treino),
      "| teste:", ncol(combat_df_teste), "amostras |", length(genes_ok), "genes)\n")
  
  # conferência opcional: diferença para o sva::ComBat no treino deve ser ~0
  # ref_sva <- sva::ComBat(as.matrix(expression_treino_bruto[genes_ok, ]),
  #                        batch = batch_treino_combat, mod = NULL)
  # max(abs(ref_sva - combat_df_treino))
  ### FIM DA CORREÇÃO ###
  
  # Salvar resultado
  readr::write_csv2(combat_df, "combat_one_ok2.xlsx")
  readr::write_csv2(metadado_merged, "metadado_merged_ok2.xlsx")
  
}


########################################################################
# 4. PCA ANTES E DEPOIS DO ComBat
#
# Amostras coloridas por série -> figuras_TCC_A/figA_pca_painel.png
########################################################################

#### PCA ######

#### PCA - painel antes e depois do ComBat ######
{
  library(ggplot2)
  library(patchwork)
  
  nomes_series <- c("1" = "GSE73408", "2" = "GSE39940", "3" = "GSE39939", "4" = "GSE37250",
                    "5" = "GSE19439", "6" = "GSE19442", "7" = "GSE19444")
  serie_pca <- factor(nomes_series[as.character(batch)], levels = nomes_series)
  
  # ---------- PCA antes e depois (mesmos genes nas duas matrizes) ----------
  pca_before <- prcomp(t(expression_matrix[genes_ok, colnames(combat_df)]), scale. = TRUE)
  pca_after  <- prcomp(t(combat_df), scale. = TRUE)
  
  var_explained_before <- round(pca_before$sdev^2 / sum(pca_before$sdev^2) * 100, 2)
  var_explained_after  <- round(pca_after$sdev^2  / sum(pca_after$sdev^2)  * 100, 2)
  
  cat("Antes : PC1 =", var_explained_before[1], "% | PC2 =", var_explained_before[2],
      "% | total =", sum(var_explained_before[1:2]), "%\n")
  cat("Depois: PC1 =", var_explained_after[1], "% | PC2 =", var_explained_after[2],
      "% | total =", sum(var_explained_after[1:2]), "%\n")
  
  # ---------- função de plot (cores e formas definidas para as 7 séries) ----------
  cores_series  <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#666666")
  formas_series <- c(16, 17, 15, 18, 3, 4, 8)
  
  plot_pca <- function(pca_obj, var_exp, titulo) {
    df <- data.frame(PC1 = pca_obj$x[, 1], PC2 = pca_obj$x[, 2], Serie = serie_pca)
    ggplot(df, aes(PC1, PC2, colour = Serie, shape = Serie)) +
      geom_point(size = 1.4, alpha = 0.75) +
      stat_ellipse(level = 0.80, linewidth = 0.5) +
      scale_colour_manual(values = cores_series, name = "Série") +
      scale_shape_manual(values = formas_series, name = "Série") +
      labs(title = titulo,
           x = paste0("PC1 (", format(var_exp[1], decimal.mark = ","), "%)"),
           y = paste0("PC2 (", format(var_exp[2], decimal.mark = ","), "%)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  
  pca_plot_before <- plot_pca(pca_before, var_explained_before, "A) Antes da correção (ComBat)")
  pca_plot_after  <- plot_pca(pca_after,  var_explained_after,  "B) Após a correção (ComBat)")
  
  # ---------- painel lado a lado, legenda única ----------
  painel_pca <- (pca_plot_before | pca_plot_after) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(painel_pca)
  ggsave("R/meta/figuras_TCC_A/figA_pca_painel.png", painel_pca,
         width = 11, height = 5, dpi = 300)
}


########################################################################
# 5. EXPRESSÃO DIFERENCIAL E ENRIQUECIMENTO FUNCIONAL
#
# - limma só com as amostras de treino (mesmo split, seed 42)
# - TREAT: |log2FC| > 1 e FDR < 0,05 -> DEG_Active_vs_Latent
# - Reactome (ReactomePA) e GO BP (clusterProfiler) para os genes
#   com expressão aumentada e diminuída na TB ativa
########################################################################

############ ---------- DIFFRENTIAL EXPRESSION ANALYSYS----------########################
{
  #### Differential Expression analysis ####
  library(limma)
  expr_data <- combat_df
  meta <- metadado_merged
 
  ### CORREÇÃO DE VAZAMENTO ###
  # PROBLEMA ORIGINAL: os DEGs (DEG_Active_vs_Latent) eram calculados usando
  # expr_data/meta com as 542 amostras (treino + teste) - ou seja, a seleção
  # de quais genes entram como candidatos para o RFE já "via" as amostras
  # de teste antes mesmo do split acontecer.
  #
  # SOLUÇÃO: restringir expr_data/meta às amostras de TREINO antes de rodar
  # limma/eBayes. Usa a mesma lógica de set.seed(42) + createDataPartition
  # do bloco "Universal Signature" mais abaixo, então o split aqui é
  # idêntico ao que será usado depois no restante do pipeline.
  
  y_deg_split <- factor(meta$Groups)
  levels(y_deg_split) <- make.names(levels(y_deg_split))
  
  set.seed(42)
  idx_treino_deg <- createDataPartition(y_deg_split, p = 0.7, list = FALSE)
  
  expr_data <- expr_data[, idx_treino_deg, drop = FALSE]
  meta      <- meta[idx_treino_deg, , drop = FALSE]
  
  cat("OK: DEGs serão calculados só com treino (", ncol(expr_data), "amostras )\n")
  ### FIM DA CORREÇÃO ###
  
  group <- factor(meta$Groups)
  
  dge <- limma::getEAWP(as.matrix(expr_data))
  
  # Design matrix (sem intercepto)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  colnames(design) <- make.names(colnames(design))
  
  # Contrast (Active vs Latent)
  
  contr.matrix <- makeContrasts(
    Active_vs_Latent = Active.tuberculosis - Latent.tuberculosis,
    levels = design
  )
  
  # Fit model
  library(edgeR)
  
  fit <- lmFit(expr_data, design)
  fit2 <- contrasts.fit(fit, contr.matrix)
  efit <- eBayes(fit2, trend = TRUE)
  plotSA(efit, main = "Final model: Mean-variance trend")
  
  DEG_Active_vs_Latent <- topTreat(
    efit,
    coef = "Active_vs_Latent",
    adjust.method = "fdr",
    p.value = 0.05,
    number = Inf,
    lfc = 1
  )
  
  # plot the graph of model and raw
  plotSA(efit,main="Final model:Mean-variance trend")
  
  summary(DEG_Active_vs_Latent)
  summary(decideTests(efit))
  
  # > summary(DEG_Active_vs_Latent)
  # logFC           AveExpr             t              P.Value            adj.P.Val               B          
  # Min.   :-1.640   Min.   : 4.267   Min.   :-16.010   Min.   :0.000e+00   Min.   :0.000e+00   Min.   :  4.349  
  # 1st Qu.: 1.006   1st Qu.: 6.192   1st Qu.:  6.746   1st Qu.:0.000e+00   1st Qu.:0.000e+00   1st Qu.: 40.057  
  # Median : 1.107   Median : 8.027   Median : 11.224   Median :0.000e+00   Median :0.000e+00   Median : 52.323  
  # Mean   : 0.740   Mean   : 8.028   Mean   :  7.102   Mean   :1.501e-08   Mean   :3.864e-08   Mean   : 59.987  
  # 3rd Qu.: 1.302   3rd Qu.: 9.652   3rd Qu.: 14.564   3rd Qu.:0.000e+00   3rd Qu.:0.000e+00   3rd Qu.: 79.509  
  # Max.   : 2.595   Max.   :13.257   Max.   : 19.670   Max.   :1.366e-06   Max.   :3.515e-06   Max.   :125.995  
  # >   summary(decideTests(efit))
  # Active_vs_Latent
  # Down               3732
  # NotSig             3690
  # Up                 2397
  # 
  openxlsx::write.xlsx(
    list(
      DEG_Active_vs_Latent = DEG_Active_vs_Latent,
      Input_data = expr_data
    ),
    file = "Active_vs_Latent_TB_DEGs_ok2.xlsx",
    rowNames = TRUE
  )
  
  ############ ---------- CLUSTERING ANALYSYS----------########################
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(openxlsx)
  library(dplyr)
  library(ReactomePA)
  
  # usa os objetos já na memória (sem reler a planilha)
  gene_list <- data.frame(SYMBOL = rownames(DEG_Active_vs_Latent), DEG_Active_vs_Latent)
  gene_list <- gene_list[order(gene_list$logFC), ]
  
  gene_log  <- setNames(gene_list$logFC, gene_list$SYMBOL)
  fc_vector <- gene_log
  
  sig_genes  <- DEG_Active_vs_Latent            # já filtrado: FDR < 0.05 e |logFC| > 1
  up_genes   <- rownames(sig_genes)[sig_genes$logFC > 0]
  down_genes <- rownames(sig_genes)[sig_genes$logFC < 0]
  length(up_genes); length(down_genes)
  
  # conversão para ENTREZ
  up_entrez   <- bitr(up_genes,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  # universo = TODOS os genes testados (não só os DEGs)
  universo_entrez <- bitr(rownames(expr_data), fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  gene_clusters <- list(Upregulated = up_entrez, Downregulated = down_entrez)
  
  #inspecionar objeto, e ver min e max
  ck <- compareCluster(
    geneCluster = gene_clusters,
    fun = "enrichPathway",
    organism = "human",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  path_data <- enrichplot::dotplot(
    ck,
    showCategory = 20,
    title = "Enriquecimento de vias - Reactome",
    font.size = 15
  )
  
  p_reactome_pt <- path_data +
    scale_y_discrete(labels = traduz_via) +
    scale_x_discrete(labels = traduz_grupo) +
    labs(x = "Grupo de DEGs", fill = "p ajustado", colour = "p ajustado", size = "Razão gênica")
  
  print(p_reactome_pt)
  
  # tradução dos nomes das vias (Reactome)
  traducao_vias <- c(
    "Interferon Signaling"            = "Sinalização por interferon",
    "Neutrophil degranulation"        = "Desgranulação de neutrófilos",
    "Interferon alpha/beta signaling" = "Sinalização por interferon alfa/beta",
    "Antigen activates B Cell Receptor (BCR) leading to generation of second messengers" =
      "Ativação do receptor de células B (BCR) por antígeno com geração de segundos mensageiros"
  )
  traduz_via <- function(x) {
    y <- ifelse(x %in% names(traducao_vias), traducao_vias[x], x)
    stringr::str_wrap(y, width = 30)          # mesma quebra de linha do dotplot
  }
  traduz_grupo <- function(x) {               # mantém o "(50)" e o "(11)"
    x <- sub("Upregulated", "Superexpressos", x)
    sub("Downregulated", "Subexpressos", x)
  }
  
  p_reactome_pt <- path_data +
    scale_y_discrete(labels = traduz_via) +
    scale_x_discrete(labels = traduz_grupo) +
    labs(x = "Grupo de DEGs", fill = "p ajustado", colour = "p ajustado", size = "Razão gênica")
  
  print(p_reactome_pt)
  
  #
  
  fc_vector <- gene_list$logFC
  names(fc_vector) <- gene_list$SYMBOL 
  
  # cnetplot(ck, categorySize="pvalue", foldChange=fc_vector)
  
  go_enrich <- compareCluster(
    geneCluster   = gene_clusters,        # list de ENTREZID, formato certo para compareCluster
    fun           = "enrichGO",
    universe      = gene_list$SYMBOL |> 
      (\(s) clusterProfiler::bitr(s, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)(),
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  )
  
  #
  
  go_down <- enrichGO(
    gene          = down_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",     # <- corrigido, era "SYMBOL"
    readable      = TRUE,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10
  )
  
  go_up <- enrichGO(
    gene          = up_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",     # <- confirme que está assim, não "SYMBOL"
    readable      = TRUE,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10
  )
  
  upsetplot(go_up)
  
}

########################################################################
# 6. SELEÇÃO DA ASSINATURA (RFE + RANDOM FOREST)
#
# - preditores: DEGs; mesmo split 70/30 (seed 42)
# - CV agrupada por estudo: cada fold deixa uma série de fora
# - RFE com 2, 3, 5, 10, 15, 20, 30 e 50 genes
# - candidatos comparados por regressão logística (AUC média,
#   AUC mínima e desvio-padrão entre os estudos)
# - modelo final RF com assinatura_1 + teste de permutação (100x)
########################################################################

#### Universal Signature - Feature selection #####
library(pROC)
library(caret)
library(randomForest)
library(pheatmap)
library(dplyr)

####---- FORMATAÇÃO DATASET ----####
t_combat_df <- as.data.frame(t(combat_df))

stopifnot(all(rownames(t_combat_df) == rownames(metadado_merged)))
cat("OK: amostras de t_combat_df e metadado_merged alinhadas\n")

y <- factor(metadado_merged$Groups)
levels(y) <- make.names(levels(y))
table(y)

genes_deg <- rownames(DEG_Active_vs_Latent)
genes_deg <- intersect(genes_deg, colnames(t_combat_df))

cat("Genes diferencialmente expressos disponíveis:", length(genes_deg), "\n")

combat_df_full <- t_combat_df[, genes_deg]

set.seed(42)
idx_treino <- createDataPartition(y, p = 0.7, list = FALSE)

X_treino <- combat_df_full[idx_treino, ]
y_treino <- y[idx_treino]
X_teste  <- combat_df_full[-idx_treino, ]
y_teste  <- y[-idx_treino]

table(y_treino)
table(y_teste)

batch_treino <- metadado_merged$batch[idx_treino]

#### ---- CRIAÇÃO DE FOLDS ---- ####
set.seed(42)
estudos_unicos <- sort(unique(batch_treino))

folds_agrupados <- lapply(estudos_unicos, function(estudo_fora) {
  which(batch_treino != estudo_fora)
})
names(folds_agrupados) <- paste0("Fold_sem_", estudos_unicos)

n_estudos_esperado <- length(unique(batch_treino))
if (length(folds_agrupados) != n_estudos_esperado) {
  stop("folds_agrupados tem ", length(folds_agrupados), " folds, esperado ",
       n_estudos_esperado, ". Reinicie a sessão e rode o script do zero.")
}
cat("OK: folds_agrupados tem", length(folds_agrupados), "folds\n")

controle_rfe <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  index = folds_agrupados,
  verbose = FALSE
)

#### ---- RFE ---- ####
tamanhos_testados <- c(2, 3, 5, 10, 15, 20, 30, 50)
tamanhos_testados <- tamanhos_testados[tamanhos_testados < ncol(X_treino)]

resultado_rfe <- rfe(
  x = X_treino,
  y = y_treino,
  sizes = tamanhos_testados,
  rfeControl = controle_rfe
)

print(resultado_rfe)
plot(resultado_rfe, type = c("g", "o"),
     main = "RFE - Acurácia por tamanho de assinatura (CV agrupada por estudo)")

tamanho_parcimonioso <- pickSizeTolerance(
  resultado_rfe$results,
  metric = "Accuracy",
  maximize = TRUE,
  tol = 15
)
cat("Tamanho de assinatura mais parcimonioso:", tamanho_parcimonioso, "\n")

melhor <- max(resultado_rfe$results$Accuracy)
pior   <- min(resultado_rfe$results$Accuracy)
resultado_rfe$results$tolerancia_pct <- 100 * (melhor - resultado_rfe$results$Accuracy) / (melhor - pior)
print(resultado_rfe$results[, c("Variables", "Accuracy", "AccuracySD", "tolerancia_pct")])


# sig_2 <- extrair_assinatura(resultado_rfe, 2)
# sig_3 <- extrair_assinatura(resultado_rfe, 3)
# 
# # gene que entra no tamanho 3
# setdiff(sig_3, sig_2)
# 
# A curva de acurácia não é monotônica em RFE — é esperado que 
# tamanhos intermediários oscilem. O que seria preocupante é se a acurácia crescesse monotonicamente com mais genes,
# pois indicaria overfitting.

# ------------------------------------------------------------------
# Função genérica: extrair assinatura de um tamanho arbitrário
# ------------------------------------------------------------------
extrair_assinatura <- function(resultado_rfe, tamanho) {
  resultado_rfe$variables %>%
    filter(Variables == tamanho) %>%
    group_by(var) %>%
    summarise(freq = n()) %>%
    arrange(desc(freq)) %>%
    slice_head(n = tamanho) %>%
    pull(var)
}

# Extrai assinatura parcimoniosa e força para 1 gene (FCGR1A)
# FCGR1B excluído por ser pseudogene; FCGR1A apresenta AUC superior
assinatura_rfe <- extrair_assinatura(resultado_rfe, tamanho_parcimonioso)
cat("Assinatura RFE original:", paste(assinatura_rfe, collapse = ", "), "\n")
#Assinatura RFE original: GBP6, AIM2 
assinatura_1 <-  "AIM2"
cat("Assinatura final (1 gene):", assinatura_1, "\n")

cat("\nEstudos presentes no treino (batch_treino)\n")
print(table(batch_treino))
cat("Número de estudos únicos no treino:", length(unique(batch_treino)), "\n")

cat("\n=== Mapeamento fold -> estudo(s) ===\n")
for (i in seq_along(folds_agrupados)) {
  idx_no_fold <- folds_agrupados[[i]]
  idx_fora    <- setdiff(seq_along(batch_treino), idx_no_fold)
  estudos_teste <- unique(batch_treino[idx_fora])
  cat("Fold", i, "- testado no(s) estudo(s):", paste(estudos_teste, collapse = ", "),
      "- N amostras teste:", length(idx_fora), "\n")
}

cat("\n=== Distribuição de classes por estudo (base completa) ===\n")
print(table(metadado_merged$batch, y))

# print(table(metadado_merged$batch, y))
# y
# Active.tuberculosis Latent.tuberculosis
# 1                  34                  35
# 2                  65                  53
# 3                  49                  14
# 4                  95                  79
# 5                  13                  17
# 6                  20                  28
# 7                  19                  21

resultado_rfe$results$SE <- resultado_rfe$results$AccuracySD / sqrt(7)
resultado_rfe$results$IC_inferior <- resultado_rfe$results$Accuracy - 1.96 * resultado_rfe$results$SE
resultado_rfe$results$IC_superior <- resultado_rfe$results$Accuracy + 1.96 * resultado_rfe$results$SE
print(resultado_rfe$results[, c("Variables", "Accuracy", "IC_inferior", "IC_superior")])

#### ---- TREINAR O MODELO FINAL (RF) COM ASSINATURA ---- ####

controle_treino <- trainControl(
  method = "cv",
  index = folds_agrupados,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

if (length(controle_treino$index) != n_estudos_esperado) {
  stop("controle_treino$index tem ", length(controle_treino$index),
       " folds, esperado ", n_estudos_esperado, ". Objeto desatualizado.")
}
cat("OK: controle_treino usa", length(controle_treino$index), "folds agrupados\n")

train <- caret::train

modelo_rf_1gene <- caret::train(
  x = X_treino[, assinatura_1, drop = FALSE],
  y = y_treino,
  method = "rf",
  metric = "ROC",
  trControl = controle_treino,
  tuneLength = 5,
  importance = TRUE
)

imp <- varImp(modelo_rf_1gene, scale = FALSE)$importance
print(imp)

# COMPARAÇÃO CANDIDATOS PARA ASSINATURA
X_cand <- t_combat_df[idx_treino, ]          # todos os genes, só amostras de treino

candidatos <- list(
  FCGR1A           = "FCGR1A",
  GBP6             = "GBP6",
  AIM2             = "AIM2",
  GBP6_AIM2        = c("GBP6", "AIM2"),
  FCGR1A_GBP6_AIM2 = c("FCGR1A", "GBP6", "AIM2")
)

comparacao <- t(sapply(candidatos, function(g) {
  set.seed(42)
  m <- caret::train(x = X_cand[, g, drop = FALSE], y = y_treino,
                    method = "glm", metric = "ROC", trControl = controle_treino)
  c(ROC_medio = mean(m$resample$ROC),     # média entre os 7 estudos deixados de fora
    ROC_min   = min(m$resample$ROC),      # pior estudo
    ROC_dp    = sd(m$resample$ROC))
}))
round(comparacao, 3)

# round(comparacao, 3)
# ROC_medio ROC_min ROC_dp
# FCGR1A               0.915   0.821  0.074
# GBP6                 0.905   0.709  0.092
# AIM2                 0.926   0.842  0.048
# GBP6_AIM2            0.933   0.806  0.063
# FCGR1A_GBP6_AIM2     0.929   0.803  0.065

library(ggplot2)

imp <- varImp(modelo_rf_1gene, scale = FALSE)$importance

# com 1 gene, o caret devolve classes nas linhas: transpõe
if (!all(rownames(imp) %in% assinatura_1)) {
  imp <- as.data.frame(t(imp))
  rownames(imp) <- assinatura_1
}
print(imp)

imp_df <- data.frame(
  gene       = rownames(imp),
  importance = rowMeans(imp)
)

ggplot(imp_df, aes(x = importance, y = reorder(gene, importance))) +
  geom_point(size = 3, color = "steelblue") +
  geom_segment(aes(x = 0, xend = importance, yend = gene), color = "steelblue") +
  labs(x = "Mean Decrease Accuracy", y = NULL,
       title = paste0("Importância dos genes na assinatura (", length(assinatura_1), " gene)")) +
  theme_classic()

# ------------------------------------------------------------------
# Teste de permutação — afirma que FCGR1A é relevante estatisticamente
# ------------------------------------------------------------------
set.seed(42)
n_permutacoes <- 100
roc_nulo <- numeric(n_permutacoes)

for (i in 1:n_permutacoes) {
  y_embaralhado <- sample(y_treino)
  
  modelo_nulo <- caret::train(
    x = X_treino[, assinatura_1, drop = FALSE],
    y = y_embaralhado,
    method = "rf",
    metric = "ROC",
    trControl = controle_treino
  )
  roc_nulo[i] <- max(modelo_nulo$results$ROC)
}

roc_observado     <- max(modelo_rf_1gene$results$ROC)
p_valor_permutacao <- mean(roc_nulo >= roc_observado)

cat("ROC observado:", round(roc_observado, 3), "\n")
cat("ROC nulo (média):", round(mean(roc_nulo), 3), "\n")
cat("p-valor (permutação):", p_valor_permutacao, "\n")

hist(roc_nulo, main = "Distribuição nula (rótulos embaralhados)",
     xlab = "ROC", xlim = c(0.3, 1))
abline(v = roc_observado, col = "red", lwd = 2)

########################################################################
# 7. AVALIAÇÃO NO CONJUNTO DE TESTE (30%)
#
# Probabilidade de TB ativa com limiar 0,5. AUC, acurácia (IC 95%),
# sensibilidade e especificidade no teste completo e por série,
# heatmap da assinatura e curvas ROC.
########################################################################

# ------------------------------------------------------------------
# Função de avaliação
# ------------------------------------------------------------------
avaliar_em_teste <- function(modelo, assinatura, X_alvo, y_alvo, titulo = "") {
  
  X_alvo_sub  <- X_alvo[, assinatura, drop = FALSE]
  pred_classe <- predict(modelo, newdata = X_alvo_sub)
  pred_prob   <- predict(modelo, newdata = X_alvo_sub, type = "prob")
  
  nivel_positivo <- levels(y_alvo)[2]
  
  cm <- confusionMatrix(pred_classe, y_alvo, positive = nivel_positivo)
  
  roc_obj <- tryCatch(
    pROC::roc(response = y_alvo, predictor = pred_prob[[nivel_positivo]],
              levels = levels(y_alvo), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  
  list(
    titulo = titulo,
    cm     = cm,
    roc    = roc_obj,
    auc    = if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA,
    n      = length(y_alvo)
  )
}

# ------------------------------------------------------------------
# Avaliação no grupo TESTE COMPLETO (holdout de 30%)
# ------------------------------------------------------------------
res_teste_1 <- avaliar_em_teste(modelo_rf_1gene, assinatura_1, X_teste, y_teste,
                                titulo = "Assinatura 1 gene - Teste completo")

cat("\n=== Assinatura 1 gene — Grupo Teste ===\n")
print(res_teste_1$cm)
cat("AUC:", round(res_teste_1$auc, 3), "\n")

# AUC: 0.821 FCRG1A
# AUC: 0.854 AIM2

# Confusion Matrix and Statistics FCR1A
# 
# Reference
# Prediction            Active.tuberculosis Latent.tuberculosis
# Active.tuberculosis                  63                  17
# Latent.tuberculosis                  25                  57
# 
# Accuracy : 0.7407             
# 95% CI : (0.6661, 0.8063)   
# No Information Rate : 0.5432             
# P-Value [Acc > NIR] : 1.761e-07          
# 
# Kappa : 0.482              
# 
# Mcnemar's Test P-Value : 0.2801             
#                                              
#             Sensitivity : 0.7703             
#             Specificity : 0.7159             
#          Pos Pred Value : 0.6951             
#          Neg Pred Value : 0.7875             
#              Prevalence : 0.4568             
#          Detection Rate : 0.3519             
#    Detection Prevalence : 0.5062             
#       Balanced Accuracy : 0.7431             
#                                              
#        'Positive' Class : Latent.tuberculosis

# ------------------------------------------------------------------
# Avaliação SEPARADA por batch/estudo
# ------------------------------------------------------------------
batch_teste <- metadado_merged$batch[-idx_treino]

avaliar_por_batch <- function(modelo, assinatura, X_teste, y_teste, batch_teste) {
  
  estudos <- sort(unique(batch_teste))
  
  resultados <- lapply(estudos, function(est) {
    idx   <- which(batch_teste == est)
    y_sub <- y_teste[idx]
    
    if (length(unique(y_sub)) < 2) {
      return(data.frame(
        batch = est, n = length(idx),
        n_classe1 = sum(y_sub == levels(y_teste)[1]),
        n_classe2 = sum(y_sub == levels(y_teste)[2]),
        acuracia = NA, AUC = NA,
        obs = "classe única no estudo - ROC não calculável"
      ))
    }
    
    r <- avaliar_em_teste(modelo, assinatura, X_teste[idx, , drop = FALSE], y_sub,
                          titulo = paste("Batch:", est))
    
    data.frame(
      batch = est, n = length(idx),
      n_classe1 = sum(y_sub == levels(y_teste)[1]),
      n_classe2 = sum(y_sub == levels(y_teste)[2]),
      acuracia = unname(r$cm$overall["Accuracy"]),
      AUC = r$auc,
      obs = ""
    )
  })
  
  bind_rows(resultados)
}

tabela_por_batch_1 <- avaliar_por_batch(modelo_rf_1gene, assinatura_1,
                                        X_teste, y_teste, batch_teste)

cat("\n=== Desempenho por batch — assinatura 1 gene ===\n")
print(tabela_por_batch_1)

# checagem batch 7
idx_b7 <- which(batch_teste == 7)

boxplot(
  X_teste[idx_b7, "FCGR1A"] ~ y_teste[idx_b7],
  main = "FCGR1A — Batch 7",
  xlab = "Classe",
  ylab = "Expressão"
)

# checagem bach 2
idx_b2 <- which(batch_teste == 2)

boxplot(
  X_teste[idx_b2, "FCGR1A"] ~ y_teste[idx_b2],
  main = "FCGR1A — Batch 2",
  xlab = "Classe",
  ylab = "Expressão"
)


# ------------------------------------------------------------------
# HEATMAP da assinatura no grupo teste
# ------------------------------------------------------------------

library(pheatmap)
library(viridis)

# z-score por gene e limite em ±2.5 (evita as faixas pretas)
mat_z <- t(scale(t(mat)))
mat_z[mat_z >  2.5] <-  2.5
mat_z[mat_z < -2.5] <- -2.5

breaks <- seq(-2.5, 2.5, length.out = 101)

# o nome da lista TEM que ser igual ao nome da coluna em anotacao_colunas
annotation_colors <- list(
  Grupo = c("Active.tuberculosis" = "#E74C3C",
            "Latent.tuberculosis" = "#2ECC71")
)

titulo_heat <- paste0("Heatmap ", paste(assinatura_1, collapse = ", "), " - Split Test")

p_heat <- pheatmap(
  mat_z,
  color             = magma(100),
  scale             = "none",          # já escalado acima
  breaks            = breaks,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  annotation_col    = anotacao_colunas,
  annotation_colors = annotation_colors,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  cellheight        = 25,              # altura fixa da linha do gene
  cellwidth         = 3,               # largura de cada amostra
  fontsize          = 9,
  border_color      = NA,
  main              = titulo_heat,
  silent            = TRUE
)

# exibir na tela
grid::grid.newpage(); grid::grid.draw(p_heat$gtable)

# salvar (o png() vem ANTES do desenho e o dev.off() depois)
png(paste0("heatmap_", paste(assinatura_1, collapse = "_"), "_teste.png"),
    width = 3000, height = 700, res = 300)
grid::grid.draw(p_heat$gtable)
dev.off()

png("heatmap_AIM2_heat.png", width = 3000, height = 600, res = 300)

dev.off()

library(ggplot2)
library(patchwork)

df_box <- data.frame(
  Expressao = as.vector(mat),
  Grupo     = y_teste[ordem]
)

p_box <- ggplot(df_box, aes(x = Grupo, y = Expressao, fill = Grupo)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("Active.tuberculosis"  = "#E74C3C",
                               "Latent.tuberculosis" = "#2ECC71")) +
  labs(y = "Expression FCGR1A", x = NULL) +
  theme_classic() +
  theme(legend.position = "none")

anotacao_colunas <- data.frame(
  Batch = factor(batch_teste),
  Group  = y_teste,
  row.names = colnames(mat)
)

while (!is.null(dev.list())) dev.off()
graphics.off()


# ------------------------------------------------------------------
# CURVA ROC — grupo teste completo
# ------------------------------------------------------------------

plot(res_teste_1$roc, col = "darkgreen", lwd = 2,
     main = paste0("Curva ROC - Assinatura 1 gene / AIM2 (Teste)\nAUC = ",
                   round(res_teste_1$auc, 3)))


#

cm_por_batch <- function(modelo, assinatura, X_teste, y_teste, batch_teste) {
  
  estudos <- sort(unique(batch_teste))
  
  lapply(estudos, function(est) {
    idx   <- which(batch_teste == est)
    y_sub <- y_teste[idx]
    
    if (length(unique(y_sub)) < 2) {
      cat("\n=== Batch", est, "— classe única, CM não calculável ===\n")
      return(invisible(NULL))
    }
    
    X_sub       <- X_teste[idx, assinatura, drop = FALSE]
    pred_classe <- predict(modelo, newdata = X_sub)
    
    cm <- confusionMatrix(pred_classe, y_sub,
                          positive = levels(y_sub)[2])
    print(cm)
  })
  
  invisible(NULL)
}

cm_por_batch(modelo_rf_1gene, assinatura_1, X_teste, y_teste, batch_teste)

#

roc_por_batch_expressao <- function(assinatura, X_completo, y_completo,
                                    batch_completo) {
  
  # vetor com os nomes dos estudos ordenados (ex: "1","2",...,"7")
  estudos   <- sort(unique(batch_completo))
  
  # número total de estudos
  n_estudos <- length(estudos)
  
  # número de colunas do grid de plots
  ncols     <- 3
  
  # número de linhas necessárias para acomodar todos os estudos
  nrows     <- ceiling(n_estudos / ncols)
  
  # configura o dispositivo gráfico em grid (nrows x ncols)
  # mar define as margens de cada plot (bottom, left, top, right)
  par(mfrow = c(nrows, ncols), mar = c(4, 4, 3, 1))
  
  # itera sobre cada estudo
  lapply(estudos, function(est) {
    
    # índices das amostras pertencentes ao estudo atual
    idx   <- which(batch_completo == est)
    
    # rótulos reais das amostras do estudo atual
    y_sub <- y_completo[idx]
    
    # se o estudo tem apenas uma classe, ROC não é calculável
    # plota um painel vazio com aviso e pula para o próximo estudo
    if (length(unique(y_sub)) < 2) {
      plot.new()
      title(main = paste("Batch", est, "— classe única"))
      return(invisible(NULL))
    }
    
    # extrai a expressão bruta do FCGR1A para as amostras do estudo
    # sem passar pelo modelo — [[é o sinal biológico puro do gene
    expressao_sub <- X_completo[idx, assinatura]
    
    # calcula a curva ROC usando expressão do gene como preditor direto
    # direction = "<" significa que valores maiores predizem a classe positiva
    # tryCatch evita que um erro quebre toda a função
    roc_obj <- tryCatch(
      pROC::roc(response  = y_sub,
                predictor = expressao_sub,
                levels    = levels(y_sub),
                direction = ">",   # AIM2 maior na TB ativa (1º nível = Active)
                quiet     = TRUE),
      error = function(e) NULL  # retorna NULL se houver erro
    )
    
    # só plota se a ROC foi calculada com sucesso
    if (!is.null(roc_obj)) {
      
      # extrai o valor numérico do AUC
      auc_val <- as.numeric(pROC::auc(roc_obj))
      
      # plota a curva ROC do estudo
      plot(roc_obj,
           col  = "darkgreen",  # cor da curva
           lwd  = 2,            # espessura da linha
           main = paste0("Batch ", est, " (n=", length(idx), ")\nAUC = ",
                         round(auc_val, 3)))
      
      # adiciona a diagonal de referência (classificador aleatório AUC=0.5)
      abline(a = 1, b = -1, lty = 2, col = "gray60")
    }
  })
  
  # restaura o dispositivo gráfico para 1 plot por janela
  par(mfrow = c(1, 1))
  
  # retorna NULL invisivelmente — a função age pelo efeito colateral dos plots
  invisible(NULL)
}

# N total de cada batch (treino + teste)
X_completo     <- combat_df_full
y_completo     <- y
batch_completo <- metadado_merged$batch

roc_por_batch_expressao(assinatura_1, X_completo, y_completo, batch_completo)

png(paste0("ROC_por_batch_", paste(assinatura_1, collapse = "_"), ".png"),
    width = 3600, height = 3600, res = 300)
roc_por_batch_expressao(assinatura_1, X_completo, y_completo, batch_completo)
dev.off()

# ------------------------------------------------------------------
#   GRID FINAL
# ------------------------------------------------------------------

# curva ROC com IC
roc_df <- data.frame(
  especificidade = 1 - res_teste_1$roc$specificities,
  sensibilidade  = res_teste_1$roc$sensitivities
)

# calcula IC
roc_ci <- pROC::ci.se(res_teste_1$roc,
                      specificities = seq(0, 1, 0.01),
                      quiet = TRUE)

ci_df <- data.frame(
  especificidade = 1 - as.numeric(rownames(roc_ci)),
  inferior       = roc_ci[, 1],
  superior       = roc_ci[, 3]
)

p_roc <- ggplot(roc_df, aes(x = especificidade, y = sensibilidade)) +
  geom_ribbon(data = ci_df,
              aes(x = especificidade, ymin = inferior, ymax = superior),
              fill = "#0d9138", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(color = "#0d9138", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  annotate("text", x = 0.75, y = 0.10,
           label = paste0("AUC = ", round(res_teste_1$auc, 3),
                          " (", round(pROC::ci(res_teste_1$roc)[1], 3),
                          "-",  round(pROC::ci(res_teste_1$roc)[3], 3), ")"),
           size = 3.5, color = "#0d9138", fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "ROC Curve — Test Set") +
  theme_classic(base_size = 10) +
  theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1)

# heatmap
p_heat <- pheatmap(
  mat,
  color             = magma(100),
  scale             = "row",
  breaks            = breaks,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  annotation_col    = anotacao_colunas,
  annotation_colors = annotation_colors,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 9,
  fontsize          = 9,
  border_color      = NA,
  main              = "Heatmap FCGR1A - Split Test",
  silent            = TRUE
)

# grid final
png("painel_AIM2_2.png", width = 3500, height = 2000, res = 300)

gridExtra::grid.arrange(
  p_heat$gtable,
  p_box,
  p_roc,
  layout_matrix = rbind(c(1, 1, 1),
                        c(2, 3, NA)),
  heights = c(1.5, 2)
)

dev.off()



########################################################################
# TESTE COMPLEMENTAR: ESTABILIDADE DO SPLIT 70/30 ENTRE SEEDS DIFERENTES
# (assinatura final fixa: FCGR1A, 1 gene)
#
# Objetivo: checar se a estratificação de pacientes (qual amostra cai em
# treino vs. teste) influencia muito o desempenho do modelo final, usando
# a assinatura JÁ ESCOLHIDA (FCGR1A) - sem refazer ComBat, DEGs ou RFE.
# Só varia o split e retreina o modelo final (1 gene) para cada seed.
#
# Script standalone (sem source()) - complementar ao pipeline principal,
# não sobrescreve nenhum objeto dele (usa sufixo "_ts" em tudo).
#
# Pré-requisitos (já devem existir na sessão, vindos do script principal):
#   combat_df_full   - matriz amostras x genes (já corrigida, restrita aos DEGs)
#   y                - fator de classe (Active/Latent), mesma ordem
#   metadado_merged  - metadados, com a coluna $batch
#   assinatura_1     - "FCGR1A" (a assinatura final já escolhida)
########################################################################

library(caret)
library(pROC)
library(dplyr)

train <- caret::train  # evita conflito de nomes (ex.: com h2o)

# ------------------------------------------------------------------
# Checagem de pré-requisitos
# ------------------------------------------------------------------

objetos_necessarios <- c("combat_df_full", "y", "metadado_merged", "assinatura_1")
faltando <- objetos_necessarios[!sapply(objetos_necessarios, exists)]

if (length(faltando) > 0) {
  stop("Objetos ausentes na sessão: ", paste(faltando, collapse = ", "),
       ". Rode o pipeline principal antes deste teste complementar.")
}

cat("OK: todos os objetos necessários encontrados na sessão.\n")
cat("Assinatura fixa usada neste teste:", assinatura_1, "\n\n")


# ------------------------------------------------------------------
# Seeds a testar (edite esta linha para mudar quais seeds avaliar)
# ------------------------------------------------------------------

seeds_teste_seed <- c(1, 42, 123, 2024, 777)


# ------------------------------------------------------------------
# Loop: para cada seed, refaz só o split + modelo final (1 gene fixo)
# ------------------------------------------------------------------

resultados_teste_seed <- list()

for (s in seeds_teste_seed) {
  
  cat("Rodando seed", s, "...\n")
  
  set.seed(s)
  idx_treino_ts <- createDataPartition(y, p = 0.7, list = FALSE)
  
  X_treino_ts <- combat_df_full[idx_treino_ts, assinatura_1, drop = FALSE]
  y_treino_ts <- y[idx_treino_ts]
  X_teste_ts  <- combat_df_full[-idx_treino_ts, assinatura_1, drop = FALSE]
  y_teste_ts  <- y[-idx_treino_ts]
  
  batch_treino_ts <- metadado_merged$batch[idx_treino_ts]
  
  # folds agrupados por estudo (leave-one-study-out), só para o treino do modelo
  estudos_ts <- sort(unique(batch_treino_ts))
  folds_ts <- lapply(estudos_ts, function(e) which(batch_treino_ts != e))
  names(folds_ts) <- paste0("Fold_sem_", estudos_ts)
  
  controle_treino_ts <- trainControl(method = "cv", index = folds_ts, classProbs = TRUE,
                                     summaryFunction = twoClassSummary, savePredictions = "final")
  
  modelo_ts <- train(x = X_treino_ts, y = y_treino_ts, method = "rf", metric = "ROC",
                     trControl = controle_treino_ts, tuneLength = 5, importance = TRUE)
  
  pred_classe_ts <- predict(modelo_ts, newdata = X_teste_ts)
  pred_prob_ts   <- predict(modelo_ts, newdata = X_teste_ts, type = "prob")
  nivel_positivo_ts <- levels(y_teste_ts)[2]
  
  cm_ts <- confusionMatrix(pred_classe_ts, y_teste_ts, positive = nivel_positivo_ts)
  
  roc_ts <- pROC::roc(response = y_teste_ts, predictor = pred_prob_ts[[nivel_positivo_ts]],
                      levels = levels(y_teste_ts), direction = "auto", quiet = TRUE)
  
  resultados_teste_seed[[as.character(s)]] <- data.frame(
    seed = s,
    n_treino = length(y_treino_ts),
    n_teste = length(y_teste_ts),
    AUC = round(as.numeric(pROC::auc(roc_ts)), 3),
    acuracia = round(unname(cm_ts$overall["Accuracy"]), 3),
    sensibilidade = round(unname(cm_ts$byClass["Sensitivity"]), 3),
    especificidade = round(unname(cm_ts$byClass["Specificity"]), 3)
  )
}


# ------------------------------------------------------------------
# Tabela resumo final
# ------------------------------------------------------------------

tabela_estabilidade_split <- bind_rows(resultados_teste_seed)

cat("\n=== Estabilidade do split 70/30 (assinatura fixa:", assinatura_1, ") ===\n")
print(tabela_estabilidade_split)

cat("\nAUC — média:", round(mean(tabela_estabilidade_split$AUC), 3),
    "| desvio-padrão:", round(sd(tabela_estabilidade_split$AUC), 3),
    "| min:", round(min(tabela_estabilidade_split$AUC), 3),
    "| max:", round(max(tabela_estabilidade_split$AUC), 3), "\n")

cat("Acurácia — média:", round(mean(tabela_estabilidade_split$acuracia), 3),
    "| desvio-padrão:", round(sd(tabela_estabilidade_split$acuracia), 3), "\n")


# ------------------------------------------------------------------
# Estratificação de N amostras por batch (treino vs. teste, por classe)
# ------------------------------------------------------------------

tabela_estratificacao <- data.frame(
  batch = metadado_merged$batch,
  grupo = y,
  conjunto = ifelse(seq_along(y) %in% idx_treino, "Treino", "Teste")
)

# Tabela cruzada: batch x conjunto (treino/teste), com contagem total
tabela_batch_conjunto <- table(tabela_estratificacao$batch, tabela_estratificacao$conjunto)
cat("=== N de amostras por batch: Treino vs. Teste ===\n")
print(tabela_batch_conjunto)

# Tabela cruzada: batch x conjunto x classe (mais detalhada)
tabela_completa <- tabela_estratificacao %>%
  dplyr::count(batch, conjunto, grupo) %>%
  tidyr::pivot_wider(names_from = c(conjunto, grupo), values_from = n, values_fill = 0)

cat("\n=== N de amostras por batch, conjunto e classe ===\n")
print(tabela_completa)



##############################################################################################

########################################################################
# TESTE COMPLEMENTAR: VALIDAÇÃO POR DATASET (LEAVE-ONE-DATASET-OUT, LODO)
#
# Objetivo: em vez do split 70/30 por PACIENTE (misturando os 7 estudos
# no treino e no teste), testar se a assinatura generaliza para um ESTUDO
# INTEIRO nunca visto. A cada rodada um dos 7 estudos fica de fora como
# teste, e a correção de batch, os DEGs e o RFE são refeitos usando só os
# outros 6 - nested completo, sem vazamento entre datasets.
#
# Correção de batch: usa sva::ComBat DE VERDADE (com mod = grupo) para
# ajustar os 6 estudos de treino - ganha a proteção do sinal biológico e
# o shrinkage empírico Bayesiano. O estudo excluído, porém, é um batch
# que NUNCA existiu no ajuste (não dá pra "congelar e projetar" como no
# script principal, onde o mesmo batch aparece nos dois lados) - por
# isso ele é projetado com a mesma reescala location-scale mais simples
# de antes, mas agora mirando no alvo certo: a média/DP da matriz JÁ
# corrigida pelo ComBat real do treino, em vez de uma média crua ingênua.
#

########################################################################

library(dplyr)
library(sva)
library(limma)
library(edgeR)
library(caret)
library(randomForest)
library(pROC)
library(openxlsx)
library(ggplot2)

train <- caret::train  # evita conflito de nomes (ex.: com h2o)

# ------------------------------------------------------------------
# Checagem de pré-requisitos
# ------------------------------------------------------------------
objetos_necessarios_lodo <- c("expression_matrix", "metadado_merged")
faltando_lodo <- objetos_necessarios_lodo[!sapply(objetos_necessarios_lodo, exists)]
if (length(faltando_lodo) > 0) {
  stop("Objetos ausentes na sessão: ", paste(faltando_lodo, collapse = ", "),
       ". Rode o script principal até a secao 'DeSeq2 e merge' (ANTES do ",
       "ComBat) antes deste script de validacao por dataset.")
}

stopifnot(all(colnames(expression_matrix) == rownames(metadado_merged)))
cat("OK: expression_matrix e metadado_merged alinhados por amostra\n")

y_lodo     <- factor(metadado_merged$Groups)
levels(y_lodo) <- make.names(levels(y_lodo))
batch_lodo <- factor(as.character(metadado_merged$batch))
estudos_lodo <- sort(unique(as.character(batch_lodo)))

cat("Estudos disponiveis para leave-one-dataset-out:",
    paste(estudos_lodo, collapse = ", "), "\n")
cat("\n=== Distribuicao de classes por estudo ===\n")
print(table(batch_lodo, y_lodo))

########################################################################
#### FUNÇÕES AUXILIARES ####
########################################################################

# ------------------------------------------------------------------
# Correção de batch: ComBat DE VERDADE (sva) nos 6 estudos de treino,
# com mod = grupo (protege sinal biológico + shrinkage empírico
# Bayesiano). O estudo excluído nunca existiu nesse ajuste - não dá pra
# "congelar e projetar" parâmetros que não existem para um batch novo -
# então ele é projetado com uma reescala location-scale simples, usando
# só a média/DP DELE MESMO (nenhum rótulo Active/Latente é usado), mas
# mirando no alvo certo: a média/DP da matriz JÁ corrigida pelo ComBat
# real do treino.
# ------------------------------------------------------------------
ajustar_correcao_batch_lodo <- function(X_treino, batch_treino, grupo_treino, X_novo_estudo) {

  # 1) ComBat de verdade no treino, com mod = grupo
  mod_treino <- model.matrix(~ grupo_treino)
  X_treino_corrigido <- sva::ComBat(
    dat = as.matrix(X_treino), batch = batch_treino, mod = mod_treino,
    par.prior = TRUE, prior.plots = FALSE
  )

  #  projeta o estudo novo (nunca visto no ajuste acima): sem rotulo,
  #  sem parametro de batch aprendido para ele - so reescala media/DP
  #  dele mesmo para o alvo real do ComBat (nao mais uma media crua)
  media_alvo <- rowMeans(X_treino_corrigido)
  dp_alvo    <- apply(X_treino_corrigido, 1, sd)

  media_novo <- rowMeans(X_novo_estudo)
  dp_novo    <- apply(X_novo_estudo, 1, sd)
  dp_novo[dp_novo == 0] <- 1e-6
  X_novo_corrigido <- ((X_novo_estudo - media_novo) / dp_novo) * dp_alvo + media_alvo

  list(treino = X_treino_corrigido, teste = X_novo_corrigido,
       media_alvo = media_alvo, dp_alvo = dp_alvo)
}

# ------------------------------------------------------------------
# DEGs (limma), calculados SÓ com as amostras de treino do fold
# ------------------------------------------------------------------
calcular_degs_lodo <- function(expr_treino, grupo_treino) {

  grupo_treino <- factor(grupo_treino)
  levels(grupo_treino) <- make.names(levels(grupo_treino))

  design <- model.matrix(~0 + grupo_treino)
  colnames(design) <- make.names(levels(grupo_treino))

  niveis <- levels(grupo_treino)
  if (length(niveis) != 2) {
    stop("calcular_degs_lodo espera exatamente 2 classes, encontrou: ",
         paste(niveis, collapse = ", "))
  }

  contr.matrix <- makeContrasts(
    contrasts = paste0(niveis[1], "-", niveis[2]),
    levels = design
  )

  fit  <- lmFit(expr_treino, design)
  fit2 <- contrasts.fit(fit, contr.matrix)
  efit <- eBayes(fit2)

  topTreat(efit, coef = 1, adjust.method = "fdr", p.value = 0.05,
           number = Inf, lfc = 1)
}

# ------------------------------------------------------------------
# RFE nested: CV agrupada pelos estudos de TREINO do fold (leave-one-
# estudo-restante-out), igual ao script principal, só que restrito aos
# 6 estudos daquele fold.
# ------------------------------------------------------------------
rfe_por_fold <- function(X_treino, y_treino, batch_treino_fold,
                         tamanhos = c(2, 3, 5, 10, 15, 20, 30, 50)) {

  estudos_fold <- sort(unique(as.character(batch_treino_fold)))
  folds_fold <- lapply(estudos_fold, function(e) {
    which(as.character(batch_treino_fold) != e)
  })
  names(folds_fold) <- paste0("Fold_sem_", estudos_fold)

  controle_rfe_fold <- rfeControl(
    functions = rfFuncs, method = "cv", index = folds_fold, verbose = FALSE
  )

  tamanhos <- tamanhos[tamanhos < ncol(X_treino)]
  if (length(tamanhos) == 0) tamanhos <- max(2, ncol(X_treino) - 1)

  resultado <- rfe(x = X_treino, y = y_treino, sizes = tamanhos,
                   rfeControl = controle_rfe_fold)

  tamanho_parcimonioso <- tryCatch(
    pickSizeTolerance(resultado$results, metric = "Accuracy",
                      maximize = TRUE, tol = 15),
    error = function(e) resultado$bestSubset
  )

  assinatura <- resultado$variables %>%
    dplyr::filter(Variables == tamanho_parcimonioso) %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(freq)) %>%
    dplyr::slice_head(n = tamanho_parcimonioso) %>%
    dplyr::pull(var)

  list(resultado_rfe = resultado, assinatura = assinatura,
       tamanho = tamanho_parcimonioso)
}

# ------------------------------------------------------------------
# Treina o modelo final do fold (CV agrupada nos 6 estudos de treino)
# e avalia no estudo excluído (teste externo daquele fold). Retorna a
# linha-resumo (para a tabela), o objeto ROC completo (para plotar) e a
# matriz de confusão completa (para inspecionar por dataset).
# ------------------------------------------------------------------
treinar_avaliar_fold <- function(X_treino, y_treino, batch_treino_fold,
                                 assinatura, X_teste, y_teste, estudo_out) {

  assinatura <- intersect(assinatura, colnames(X_treino))
  if (length(assinatura) == 0) {
    warning("Fold ", estudo_out, ": nenhum gene da assinatura disponivel - ",
            "fold ignorado.")
    return(NULL)
  }

  estudos_fold <- sort(unique(as.character(batch_treino_fold)))
  folds_fold <- lapply(estudos_fold, function(e) {
    which(as.character(batch_treino_fold) != e)
  })
  names(folds_fold) <- paste0("Fold_sem_", estudos_fold)

  controle_treino_fold <- trainControl(
    method = "cv", index = folds_fold, classProbs = TRUE,
    summaryFunction = twoClassSummary, savePredictions = "final"
  )

  modelo <- caret::train(
    x = X_treino[, assinatura, drop = FALSE], y = y_treino,
    method = "rf", metric = "ROC", trControl = controle_treino_fold,
    tuneLength = 5, importance = TRUE
  )

  X_teste_sub <- X_teste[, assinatura, drop = FALSE]
  pred_classe <- predict(modelo, newdata = X_teste_sub)
  pred_prob   <- predict(modelo, newdata = X_teste_sub, type = "prob")
  nivel_positivo <- levels(y_teste)[2]

  linha_base <- data.frame(
    estudo_out = estudo_out,
    n_treino = length(y_treino),
    n_teste = length(y_teste),
    n_genes_assinatura = length(assinatura),
    genes_assinatura = paste(assinatura, collapse = ";")
  )

  if (length(unique(y_teste)) < 2) {
    linha <- cbind(linha_base, data.frame(
      AUC = NA, acuracia = NA, sensibilidade = NA, especificidade = NA,
      obs = "classe unica no estudo de teste - AUC/CM nao calculavel"
    ))
    return(list(linha = linha, roc = NULL, cm = NULL))
  }

  cm <- confusionMatrix(pred_classe, y_teste, positive = nivel_positivo)
  roc_obj <- tryCatch(
    pROC::roc(response = y_teste, predictor = pred_prob[[nivel_positivo]],
              levels = levels(y_teste), direction = "auto", quiet = TRUE),
    error = function(e) NULL
  )

  linha <- cbind(linha_base, data.frame(
    AUC = if (!is.null(roc_obj)) round(as.numeric(pROC::auc(roc_obj)), 3) else NA,
    acuracia = round(unname(cm$overall["Accuracy"]), 3),
    sensibilidade = round(unname(cm$byClass["Sensitivity"]), 3),
    especificidade = round(unname(cm$byClass["Specificity"]), 3),
    obs = ""
  ))

  list(linha = linha, roc = roc_obj, cm = cm)
}

# ------------------------------------------------------------------
# Plota as curvas ROC de todos os folds sobrepostas - uma curva por
# estudo excluído, pra comparar visualmente a generalização entre
# datasets (mesma assinatura em todas as curvas de um mesmo gráfico).
# ------------------------------------------------------------------
plotar_roc_por_dataset <- function(rocs_por_fold, titulo) {
  rocs_validos <- rocs_por_fold[!sapply(rocs_por_fold, is.null)]
  if (length(rocs_validos) == 0) {
    warning("Nenhuma curva ROC valida para plotar.")
    return(invisible(NULL))
  }

  cores <- grDevices::rainbow(length(rocs_validos))
  nomes <- names(rocs_validos)

  plot(rocs_validos[[1]], col = cores[1], lwd = 2, main = titulo)
  if (length(rocs_validos) > 1) {
    for (i in 2:length(rocs_validos)) {
      plot(rocs_validos[[i]], col = cores[i], lwd = 2, add = TRUE)
    }
  }

  legenda <- sprintf("Estudo %s (AUC=%.3f)", nomes,
                     sapply(rocs_validos, function(r) as.numeric(pROC::auc(r))))
  legend("bottomright", legend = legenda, col = cores, lwd = 2, cex = 0.8)
  invisible(NULL)
}

# ------------------------------------------------------------------
# Imprime a matriz de confusão completa de cada estudo excluído
# (mesma ideia do cm_por_batch do script principal)
# ------------------------------------------------------------------
cm_por_dataset <- function(cms_por_fold) {
  for (estudo_out in names(cms_por_fold)) {
    cm <- cms_por_fold[[estudo_out]]
    cat("\n=== Matriz de confusão - estudo excluído:", estudo_out, "===\n")
    if (is.null(cm)) {
      cat("(classe única no teste - matriz de confusão não calculável)\n")
    } else {
      print(cm)
    }
  }
  invisible(NULL)
}

########################################################################
#### LOOP 1 - nested: descobre a assinatura em cada fold e mede o    ####
#### desempenho de cada fold usando a assinatura DAQUELE fold        ####
########################################################################

resultados_rfe_por_fold   <- list()
rocs_rfe_por_fold         <- list()
cms_rfe_por_fold          <- list()
assinaturas_por_fold      <- list()
dados_corrigidos_por_fold <- list()   # reaproveitado no LOOP 2
genes_deg_por_fold        <- list()

set.seed(42)

for (estudo_out in estudos_lodo) {

  cat("\n=====================================================\n")
  cat("Fold: estudo excluido =", estudo_out, "\n")
  cat("=====================================================\n")

  idx_treino <- which(as.character(batch_lodo) != estudo_out)
  idx_teste  <- which(as.character(batch_lodo) == estudo_out)

  X_treino_bruto <- expression_matrix[, idx_treino, drop = FALSE]
  X_teste_bruto  <- expression_matrix[, idx_teste,  drop = FALSE]
  y_treino_fold  <- droplevels(y_lodo[idx_treino])
  y_teste_fold   <- droplevels(y_lodo[idx_teste])
  batch_treino_fold <- droplevels(batch_lodo[idx_treino])

  # 1) correcao de batch (ComBat de verdade no treino, com mod = grupo),
  #    projetando o teste sem usar o rotulo dele
  correcao <- tryCatch(
    ajustar_correcao_batch_lodo(X_treino_bruto, batch_treino_fold,
                                y_treino_fold, X_teste_bruto),
    error = function(e) {
      warning("Fold ", estudo_out, ": erro no ComBat - ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(correcao)) next
  dados_corrigidos_por_fold[[estudo_out]] <- correcao

  # 2) DEGs, calculados so com o treino (rotulo do estudo excluido nunca
  #    entra aqui)
  degs_fold <- tryCatch(
    calcular_degs_lodo(correcao$treino, y_treino_fold),
    error = function(e) {
      warning("Fold ", estudo_out, ": erro ao calcular DEGs - ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(degs_fold) || nrow(degs_fold) < 2) {
    warning("Fold ", estudo_out, ": menos de 2 DEGs encontrados - fold pulado.")
    next
  }

  genes_deg_por_fold[[estudo_out]] <- rownames(degs_fold)
  genes_candidatos <- intersect(rownames(degs_fold), rownames(correcao$treino))
  cat("DEGs encontrados no treino (sem o estudo", estudo_out, "):",
      length(genes_candidatos), "\n")

  X_treino_deg <- as.data.frame(t(correcao$treino[genes_candidatos, , drop = FALSE]))
  X_teste_deg  <- as.data.frame(t(correcao$teste[genes_candidatos, , drop = FALSE]))

  # RFE nested (CV agrupada pelos 6 estudos de treino do fold)
  rfe_fold <- tryCatch(
    rfe_por_fold(X_treino_deg, y_treino_fold, batch_treino_fold),
    error = function(e) {
      warning("Fold ", estudo_out, ": erro no RFE - ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(rfe_fold)) next

  assinaturas_por_fold[[estudo_out]] <- rfe_fold$assinatura
  cat("Assinatura escolhida neste fold:",
      paste(rfe_fold$assinatura, collapse = ", "), "\n")

  # treina o modelo final do fold com a assinatura DO PROPRIO fold e
  # avalia no estudo excluido
  resultado_fold <- treinar_avaliar_fold(
    X_treino = X_treino_deg, y_treino = y_treino_fold,
    batch_treino_fold = batch_treino_fold, assinatura = rfe_fold$assinatura,
    X_teste = X_teste_deg, y_teste = y_teste_fold, estudo_out = estudo_out
  )
  if (is.null(resultado_fold)) next

  resultados_rfe_por_fold[[estudo_out]] <- resultado_fold$linha
  rocs_rfe_por_fold[[estudo_out]]       <- resultado_fold$roc
  cms_rfe_por_fold[[estudo_out]]        <- resultado_fold$cm
}

tabela_rfe_por_fold <- dplyr::bind_rows(resultados_rfe_por_fold)
cat("\n=== Desempenho por dataset (assinatura re-selecionada em cada fold) ===\n")
print(tabela_rfe_por_fold)

########################################################################
#### CONSENSO - qual(is) gene(s) sao mais consistentes entre os     ####
#### folds
########################################################################

freq_genes <- sort(table(unlist(assinaturas_por_fold)), decreasing = TRUE)
cat("\n=== Frequencia dos genes selecionados entre os",
    length(assinaturas_por_fold), "folds validos ===\n")
print(freq_genes)

assinatura_final_lodo <- names(freq_genes)[freq_genes == max(freq_genes)]
cat("\nGene(s) mais consistente(s) entre os folds (assinatura final fixa):",
    paste(assinatura_final_lodo, collapse = ", "),
    "- selecionado(s) em", max(freq_genes), "de",
    length(assinaturas_por_fold), "folds\n")

########################################################################
#### LOOP 2 - avaliacao "honesta" por dataset, usando a MESMA       ####
#### assinatura final (fixa) em todos os 7 folds - sem re-selecionar ####
#### genes. Reaproveita a correcao de batch ja calculada no LOOP 1  ####
#### (cada fold ja tinha o estudo de teste projetado sem vazamento). ####
########################################################################

resultados_fixo_por_fold <- list()
rocs_fixo_por_fold       <- list()
cms_fixo_por_fold        <- list()

for (estudo_out in estudos_lodo) {

  correcao <- dados_corrigidos_por_fold[[estudo_out]]
  if (is.null(correcao)) next

  idx_treino <- which(as.character(batch_lodo) != estudo_out)
  idx_teste  <- which(as.character(batch_lodo) == estudo_out)
  y_treino_fold <- droplevels(y_lodo[idx_treino])
  y_teste_fold  <- droplevels(y_lodo[idx_teste])
  batch_treino_fold <- droplevels(batch_lodo[idx_treino])

  genes_disponiveis <- intersect(assinatura_final_lodo, rownames(correcao$treino))
  if (length(genes_disponiveis) == 0) {
    warning("Fold ", estudo_out, ": assinatura final nao disponivel neste fold.")
    next
  }

  X_treino_fixo <- as.data.frame(t(correcao$treino[genes_disponiveis, , drop = FALSE]))
  X_teste_fixo  <- as.data.frame(t(correcao$teste[genes_disponiveis, , drop = FALSE]))

  resultado_fold <- treinar_avaliar_fold(
    X_treino = X_treino_fixo, y_treino = y_treino_fold,
    batch_treino_fold = batch_treino_fold, assinatura = genes_disponiveis,
    X_teste = X_teste_fixo, y_teste = y_teste_fold, estudo_out = estudo_out
  )
  if (is.null(resultado_fold)) next

  resultados_fixo_por_fold[[estudo_out]] <- resultado_fold$linha
  rocs_fixo_por_fold[[estudo_out]]       <- resultado_fold$roc
  cms_fixo_por_fold[[estudo_out]]        <- resultado_fold$cm
}

tabela_fixo_por_fold <- dplyr::bind_rows(resultados_fixo_por_fold)
cat("\n=== Desempenho por dataset - ASSINATURA FIXA (",
    paste(assinatura_final_lodo, collapse = ", "), ") ===\n")
print(tabela_fixo_por_fold)

if (nrow(tabela_fixo_por_fold) > 0 && any(!is.na(tabela_fixo_por_fold$AUC))) {
  cat("\nAUC entre datasets - media:",
      round(mean(tabela_fixo_por_fold$AUC, na.rm = TRUE), 3),
      "| desvio-padrao:", round(sd(tabela_fixo_por_fold$AUC, na.rm = TRUE), 3),
      "| min:", round(min(tabela_fixo_por_fold$AUC, na.rm = TRUE), 3),
      "| max:", round(max(tabela_fixo_por_fold$AUC, na.rm = TRUE), 3), "\n")
}

########################################################################
#### AVALIAÇÕES POR BATCH/DATASET: matriz de confusão completa de    ####
#### cada estudo excluído + curvas ROC sobrepostas para comparar     ####
#### visualmente a generalização entre datasets                     ####
########################################################################

cat("\n\n################ MATRIZES DE CONFUSÃO POR DATASET ################\n")
cat("--- Assinatura re-selecionada em cada fold ---\n")
cm_por_dataset(cms_rfe_por_fold)
cat("\n--- Assinatura final fixa (", paste(assinatura_final_lodo, collapse = ", "), ") ---\n")
cm_por_dataset(cms_fixo_por_fold)

png("lodo_roc_por_dataset.png", width = 2400, height = 2000, res = 300)
plotar_roc_por_dataset(
  rocs_fixo_por_fold,
  paste0("Curvas ROC por estudo excluído - assinatura fixa: ",
         paste(assinatura_final_lodo, collapse = ", "))
)
dev.off()

##### >>> FIGURA TCC - LODO: espaço da assinatura + curvas ROC (abordagem A) #####
library(ggplot2)
library(dplyr)
library(patchwork)

pasta_fig <- "R/meta/figuras_TCC_A"
dir.create(pasta_fig, showWarnings = FALSE, recursive = TRUE)

nomes_gse <- c("1" = "GSE73408", "2" = "GSE39940", "3" = "GSE39939", "4" = "GSE37250",
               "5" = "GSE19439", "6" = "GSE19442", "7" = "GSE19444")
rotulo_estudo <- function(e) ifelse(e %in% names(nomes_gse), nomes_gse[e], e)

genes_xy <- assinatura_final_lodo[1:2]            # AIM2 e LHFPL2

# ---- 1) dados: amostras do estudo excluído em cada rodada ----
amostras <- list(); curvas <- list(); pontos <- list(); aucs <- c()
for (estudo_out in names(rocs_fixo_por_fold)) {
  roc_fold <- rocs_fixo_por_fold[[estudo_out]]
  correcao <- dados_corrigidos_por_fold[[estudo_out]]
  if (is.null(roc_fold) || is.null(correcao)) next
  
  idx_teste <- which(as.character(batch_lodo) == estudo_out)
  grupo     <- ifelse(grepl("^Active", as.character(y_lodo[idx_teste])), "TB ativa", "ILTB")
  prob_tb   <- 1 - roc_fold$predictor              # o fold guardou P(Latent)
  if (length(prob_tb) != length(idx_teste)) stop("Tamanhos diferentes no fold ", estudo_out)
  
  nome <- rotulo_estudo(estudo_out)
  amostras[[nome]] <- data.frame(estudo = nome, grupo = grupo,
                                 x = as.numeric(correcao$teste[genes_xy[1], ]),
                                 y = as.numeric(correcao$teste[genes_xy[2], ]))
  roc_tb <- pROC::roc(response = factor(grupo, levels = c("ILTB", "TB ativa")),
                      predictor = prob_tb, levels = c("ILTB", "TB ativa"),
                      direction = "<", quiet = TRUE)
  aucs[nome] <- as.numeric(pROC::auc(roc_tb))
  curvas[[nome]] <- data.frame(estudo = nome, fpr = 1 - roc_tb$specificities,
                               tpr = roc_tb$sensitivities) %>% arrange(fpr, tpr)
  pontos[[nome]] <- data.frame(estudo = nome,                       # limiar 0,5
                               fpr = mean(prob_tb[grupo == "ILTB"] >= 0.5),
                               tpr = mean(prob_tb[grupo == "TB ativa"] >= 0.5))
}
df_amostras <- bind_rows(amostras); df_curvas <- bind_rows(curvas); df_pontos <- bind_rows(pontos)

ordem  <- names(aucs)
cores  <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#D55E00",
                     "#0072B2", "#CC79A7", "#7F7F7F")[seq_along(ordem)], ordem)
rotulo_leg <- setNames(sprintf("%s (AUC = %.3f)", ordem, aucs[ordem]), ordem)
df_amostras$estudo <- factor(df_amostras$estudo, levels = ordem)
df_curvas$estudo   <- factor(df_curvas$estudo,   levels = ordem)
df_pontos$estudo   <- factor(df_pontos$estudo,   levels = ordem)

print(round(aucs, 3))   # deve bater com a Figura 8 atual

# ---- 2) tema minimalista ----
tema_min <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.6),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(), legend.key = element_blank(),
        legend.background = element_blank(), legend.text = element_text(size = 9))

# ---- 3) painel A: amostras (cheio = TB ativa; vazado = ILTB) ----
p_esp <- ggplot(df_amostras, aes(x, y, colour = estudo, shape = grupo)) +
  geom_point(size = 1.9, alpha = 0.55, stroke = 0.7) +
  scale_colour_manual(values = cores, labels = rotulo_leg) +
  scale_shape_manual(values = c("TB ativa" = 16, "ILTB" = 1)) +
  guides(colour = "none",
         shape  = guide_legend(override.aes = list(size = 2.5, colour = "grey30", alpha = 1))) +
  labs(title = "Amostras por série excluída", x = genes_xy[1], y = genes_xy[2]) +
  tema_min + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

# ---- 4) painel B: curvas ROC + ponto no limiar 0,5 ----
p_roc <- ggplot(df_curvas, aes(fpr, tpr, colour = estudo)) +
  geom_abline(linetype = "dotted", colour = "grey60") +
  geom_step(linewidth = 0.8, direction = "vh") +
  geom_point(data = df_pontos, size = 2.4) +
  scale_colour_manual(values = cores, labels = rotulo_leg) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0.02, 0.02)) +
  coord_equal() +
  labs(title = "Curvas ROC", x = "1 − Especificidade", y = "Sensibilidade") +
  tema_min + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
                   legend.key.height = unit(0.4, "cm"))

p_lodo_fig <- p_esp | p_roc
print(p_lodo_fig)
ggsave(file.path(pasta_fig, "figA_lodo_espaco_roc.png"), p_lodo_fig,
       width = 11, height = 5.2, dpi = 300)

pasta_fig <- "C:/Users/Admin/Documents/R/meta/figuras_TCC_A"
dir.create(pasta_fig, showWarnings = TRUE, recursive = TRUE)
dir.exists(pasta_fig)   # tem que dar TRUE

ggsave(file.path(pasta_fig, "figA_lodo_espaco_roc_min.png"), p_lodo_fig,
       width = 11, height = 5.2, dpi = 300)

########################################################################
#### SAÍDAS: planilha + gráfico comparando desempenho entre datasets ####
########################################################################

openxlsx::write.xlsx(
  list(
    Desempenho_assinatura_por_fold = tabela_rfe_por_fold,
    Frequencia_genes               = as.data.frame(freq_genes),
    Desempenho_assinatura_fixa     = tabela_fixo_por_fold
  ),
  file = "validacao_leave_one_dataset_out.xlsx"
)

if (nrow(tabela_fixo_por_fold) > 0) {
  p_lodo <- ggplot(tabela_fixo_por_fold, aes(x = estudo_out, y = AUC)) +
    geom_col(fill = "#2E86AB") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste0("Validacao Leave-One-Dataset-Out - assinatura fixa: ",
                     paste(assinatura_final_lodo, collapse = ", ")),
      subtitle = "Cada barra = 1 estudo inteiro nunca visto no treino daquela rodada",
      x = "Estudo excluido (teste)", y = "AUC"
    ) +
    theme_classic()
  print(p_lodo)
  ggsave("lodo_auc_por_dataset.png", p_lodo, width = 8, height = 5, dpi = 300)
}

cat("\nConcluido. Compare 'tabela_fixo_por_fold' (validacao por dataset, sem",
    "vazamento) com o resultado de teste do script principal (split 70/30 por",
    "paciente) para ver se a diferenca de estrategia de validacao muda a",
    "conclusao sobre generalizacao da assinatura.\n")

########################################################################
# 8. ANÁLISE DE COEXPRESSÃO (CEMiTool)
#
# Matriz corrigida com as 542 amostras; rede signed, Spearman.
# GSEA dos módulos (TB ativa x ILTB), módulo do FCGR1A, correlação
# do FCGR1A com os demais genes e GO BP dos correlacionados.
# Saídas: cemitool_resultados.xlsx e cemitool_report/
########################################################################

library(CEMiTool)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------------
# 1. PREPARAÇÃO DOS DADOS PARA O CEMITOOL
# ------------------------------------------------------------------

# amostras na mesma ordem nas duas tabelas?
all(rownames(metadado_merged) == colnames(combat_df))   # deve dar TRUE
# se der FALSE:
# metadado_merged <- metadado_merged[colnames(combat_df), ]

# lotes e grupos
table(metadado_merged$batch, metadado_merged$Groups)

# escala de cada lote (todos devem estar em log2, ~0 a 16)
estudos <- split(colnames(expression_matrix), metadado_merged$batch)
t(sapply(estudos, function(s) summary(as.vector(as.matrix(expression_matrix[, s])))))

# CEMiTool espera genes em linhas e amostras em colunas
expr_cemitool <- as.data.frame(combat_df)

meta_cemitool <- data.frame(
  SampleName = colnames(expr_cemitool),
  Class      = metadado_merged$Groups
)

cat("Dimensões da matriz para CEMiTool:", dim(expr_cemitool), "\n")
cat("Distribuição de classes:\n")
print(table(meta_cemitool$Class))

# ------------------------------------------------------------------
# 2. RODAR O CEMITOOL
# ------------------------------------------------------------------

rodar_cem <- function(...) {
  cemitool(
    expr               = expr_cemitool,
    annot              = meta_cemitool,
    sample_name_column = "SampleName",
    class_column       = "Class",
    filter             = TRUE,
    filter_pval        = 0.1,
    apply_vst          = FALSE,      # dados já em log2
    network_type       = "signed",
    cor_method         = "spearman",
    verbose            = TRUE,
    ...
  )
}

cem <- rodar_cem()
if (nmodules(cem) == 0) cem <- rodar_cem(force_beta = TRUE)   # plano B

cat("\nMódulos encontrados:", nmodules(cem), "\n")
print(show(cem))

# GSEA dos módulos (ativação em TB ativa x latente)
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")

# ------------------------------------------------------------------
# EM QUAL MÓDULO ESTÁ O FCGR1A
# ------------------------------------------------------------------

modulos_df <- module_genes(cem)

modulo_fcgr1a <- modulos_df %>%
  filter(genes == "FCGR1A") %>%
  pull(modules)

if (length(modulo_fcgr1a) == 0) {
  cat("\nFCGR1A não foi alocado em nenhum módulo (pode ter sido filtrado)\n")
} else {
  cat("\nFCGR1A pertence ao módulo:", modulo_fcgr1a, "\n")
  
  genes_modulo_fcgr1a <- modulos_df %>%
    filter(modules == modulo_fcgr1a) %>%
    pull(genes)
  
  cat("Número de genes no módulo", modulo_fcgr1a, ":", length(genes_modulo_fcgr1a), "\n")
  cat("Primeiros genes do módulo:", paste(head(genes_modulo_fcgr1a, 10), collapse = ", "), "\n")
}

# ------------------------------------------------------------------
# 4. CORRELAÇÃO DO FCGR1A COM TODOS OS GENES
# ------------------------------------------------------------------

cor_fcgr1a <- cor(t(as.matrix(expr_cemitool)),
                  as.numeric(expr_cemitool["FCGR1A", ]),
                  method = "spearman")[, 1]

cor_df <- data.frame(gene = names(cor_fcgr1a), rho = cor_fcgr1a) %>%
  filter(gene != "FCGR1A") %>%
  arrange(desc(rho))

cat("\n=== Top 20 genes positivamente correlacionados com FCGR1A ===\n")
print(head(cor_df, 20))
cat("\n=== Top 20 genes negativamente correlacionados com FCGR1A ===\n")
print(tail(cor_df, 20))

# ------------------------------------------------------------------
# 5. ENRIQUECIMENTO — genes correlacionados up e down
# ------------------------------------------------------------------

genes_up_cor   <- cor_df %>% filter(rho >  0.5) %>% pull(gene)
genes_down_cor <- cor_df %>% filter(rho < -0.5) %>% pull(gene)

cat("\nGenes correlacionados positivamente (rho > 0.5):", length(genes_up_cor), "\n")
cat("Genes correlacionados negativamente (rho < -0.5):", length(genes_down_cor), "\n")

up_entrez_cor <- clusterProfiler::bitr(
  genes_up_cor, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
)$ENTREZID

down_entrez_cor <- clusterProfiler::bitr(
  genes_down_cor, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
)$ENTREZID

go_up_cor <- enrichGO(
  gene = up_entrez_cor, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE
)

go_down_cor <- enrichGO(
  gene = down_entrez_cor, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE
)

if (!is.null(go_up_cor) && nrow(go_up_cor) > 0) {
  print(dotplot(go_up_cor, showCategory = 20,
                title = "GO BP — Genes positivamente correlacionados com FCGR1A"))
}
if (!is.null(go_down_cor) && nrow(go_down_cor) > 0) {
  print(dotplot(go_down_cor, showCategory = 20,
                title = "GO BP — Genes negativamente correlacionados com FCGR1A"))
}

# ------------------------------------------------------------------
# 6. ATIVIDADE DO MÓDULO DO FCGR1A
# ------------------------------------------------------------------

# NES de cada módulo em cada classe (resultado da GSEA)
nes_modulos <- cem@enrichment$nes
print(nes_modulos)

# escore do módulo por amostra (média dos z-scores dos genes do módulo)
if (length(modulo_fcgr1a) > 0 && modulo_fcgr1a != "Not.Correlated") {
  
  z_mod <- t(scale(t(as.matrix(expr_cemitool[genes_modulo_fcgr1a, ]))))
  
  activity_df <- data.frame(
    sample = colnames(expr_cemitool),
    score  = colMeans(z_mod),
    grupo  = metadado_merged$Groups
  )
  
  p_activity <- ggplot(activity_df, aes(x = grupo, y = score, fill = grupo)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.1, size = 0.5, alpha = 0.4) +
    scale_fill_manual(values = c("Active tuberculosis" = "#E74C3C",
                                 "Latent tuberculosis" = "#2ECC71")) +
    labs(title = paste0("Escore do módulo ", modulo_fcgr1a, " (contém FCGR1A)"),
         x = NULL, y = "Escore do módulo (z-score médio)") +
    theme_classic() +
    theme(legend.position = "none")
  
  print(p_activity)
}

# ------------------------------------------------------------------
# 7. SALVA RESULTADOS
# ------------------------------------------------------------------

openxlsx::write.xlsx(
  list(
    Modulos           = modulos_df,
    NES_modulos       = data.frame(modulo = rownames(nes_modulos), nes_modulos,
                                   row.names = NULL),
    Pvalor_GSEA       = data.frame(modulo = rownames(cem@enrichment$pval),
                                   cem@enrichment$pval, row.names = NULL),
    Correlacao_FCGR1A = cor_df
  ),
  file = "cemitool_resultados.xlsx"
)

# relatório HTML do CEMiTool
generate_report(cem, directory = "cemitool_report", force = TRUE)
cat("\nResultados salvos em cemitool_resultados.xlsx e cemitool_report/\n")


########################################################################
# 9. FIGURAS DO TCC - ABORDAGEM A
#
# Rodar depois do pipeline completo, na mesma sessão.
# Objetos usados: lista logo abaixo.
########################################################################

#### FIGURAS DO TCC - ABORDAGEM A ####
#   resultado_rfe, modelo_rf_1gene, assinatura_1, X_treino, y_treino, controle_treino,
#   X_teste, y_teste, batch_teste, tabela_por_batch_1, res_teste_1,
#   cem, modulo_fcgr1a
# Saída: R/meta/figuras_TCC_A/ (mesmos nomes das figuras do Word)

library(ggplot2)
library(dplyr)
library(caret)
library(patchwork)
library(viridis)

pasta_fig <- "R/meta/figuras_TCC_A"
dir.create(pasta_fig, showWarnings = FALSE, recursive = TRUE)

nomes_series <- c("1" = "GSE73408", "2" = "GSE39940", "3" = "GSE39939", "4" = "GSE37250",
                  "5" = "GSE19439", "6" = "GSE19442", "7" = "GSE19444")
rotulo_grupo <- function(y) factor(ifelse(y == "Active.tuberculosis" | y == "Active tuberculosis",
                                          "TB ativa", "ILTB"), levels = c("TB ativa", "ILTB"))
cores_grupo  <- c("TB ativa" = "#E74C3C", "ILTB" = "#2ECC71")
virgula      <- function(x) format(x, decimal.mark = ",", nsmall = 3)

# ------------------------------------------------------------------
# FIGURA X - Acurácia da RFE por número de genes (IC 95%)
# ------------------------------------------------------------------
n_folds <- length(controle_treino$index)
rfe_df <- resultado_rfe$results %>%
  mutate(SE = AccuracySD / sqrt(n_folds),
         IC_inf = Accuracy - 1.96 * SE,
         IC_sup = Accuracy + 1.96 * SE,
         pos = seq_along(Variables))

tam_sel <- tamanho_parcimonioso
genes_sel <- paste(extrair_assinatura(resultado_rfe, tam_sel), collapse = " + ")
ponto_sel <- rfe_df %>% filter(Variables == tam_sel)

p_rfe <- ggplot(rfe_df, aes(pos, Accuracy)) +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = 0.15, colour = "#2166AC") +
  geom_line(colour = "#2166AC", linewidth = 0.8) +
  geom_point(colour = "#2166AC", size = 2.5) +
  geom_point(data = ponto_sel, shape = 21, size = 5, stroke = 1.2,
             colour = "#E74C3C", fill = NA) +
  annotate("segment", x = ponto_sel$pos + 2.2, y = min(rfe_df$IC_inf) - 0.004,
           xend = ponto_sel$pos + 0.12, yend = ponto_sel$Accuracy - 0.004,
           colour = "#E74C3C", arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text", x = ponto_sel$pos + 0.3, y = min(rfe_df$IC_inf) - 0.008, hjust = 0,
           size = 3, colour = "#E74C3C",
           label = paste0("Tamanho selecionado: ", tam_sel, " genes (", genes_sel, ")")) +
  scale_x_continuous(breaks = rfe_df$pos, labels = rfe_df$Variables) +
  scale_y_continuous(labels = function(x) format(x, decimal.mark = ",", nsmall = 3)) +
  labs(x = "Número de genes", y = "Acurácia (validação cruzada por estudo)") +
  theme_classic()
ggsave(file.path(pasta_fig, "figA_rfe_acuracia.png"), p_rfe, width = 6.5, height = 4, dpi = 300)

# ------------------------------------------------------------------
# FIGURA Y - AUC do AIM2 por série no conjunto de teste
# ------------------------------------------------------------------
auc_serie <- tabela_por_batch_1 %>%
  mutate(serie = factor(paste0(nomes_series[as.character(batch)], " (n = ", n, ")"),
                        levels = rev(paste0(nomes_series[as.character(batch)], " (n = ", n, ")"))))
auc_total <- res_teste_1$auc

p_auc <- ggplot(auc_serie, aes(AUC, serie)) +
  geom_segment(aes(x = 0.5, xend = AUC, yend = serie), colour = "grey75") +
  geom_point(aes(size = n), colour = "#1B7837") +
  geom_text(aes(label = virgula(round(AUC, 3))), nudge_x = 0.03, size = 3) +
  geom_vline(xintercept = auc_total, linetype = "dashed", colour = "#E74C3C") +
  annotate("text", x = auc_total, y = nrow(auc_serie) + 0.6, hjust = 0.5, size = 3,
           colour = "#E74C3C", label = paste0("AUC no teste completo = ", virgula(round(auc_total, 3)))) +
  scale_x_continuous(limits = c(0.5, 1.02), labels = function(x) format(x, decimal.mark = ",")) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = paste0("AUC do ", assinatura_1), y = NULL) +
  theme_classic()
ggsave(file.path(pasta_fig, "figA_auc_por_serie.png"), p_auc, width = 6.5, height = 4, dpi = 300)

# ------------------------------------------------------------------
# FIGURA Z - Matrizes de confusão do AIM2 e do FCGR1A (conjunto de teste)
# ------------------------------------------------------------------
# modelo de 1 gene com o FCGR1A, treinado da mesma forma que o do AIM2
set.seed(42)
modelo_rf_fcgr1a <- caret::train(
  x = X_treino[, "FCGR1A", drop = FALSE], y = y_treino,
  method = "rf", metric = "ROC", trControl = controle_treino,
  tuneLength = 5, importance = TRUE
)
res_teste_fcgr1a <- avaliar_em_teste(modelo_rf_fcgr1a, "FCGR1A", X_teste, y_teste,
                                     titulo = "FCGR1A - Teste completo")
cat("AUC FCGR1A no teste:", round(res_teste_fcgr1a$auc, 3), "\n")
print(res_teste_fcgr1a$cm)

tabela_cm <- function(res, nome) {
  as.data.frame(res$cm$table) %>%
    mutate(Predito = rotulo_grupo(Prediction), Real = rotulo_grupo(Reference),
           modelo = nome)
}
cm_df <- bind_rows(tabela_cm(res_teste_1, assinatura_1),
                   tabela_cm(res_teste_fcgr1a, "FCGR1A")) %>%
  mutate(modelo = factor(modelo, levels = c(assinatura_1, "FCGR1A")))

p_cm <- ggplot(cm_df, aes(Real, Predito, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Freq, colour = Freq > 45), size = 6, fontface = "bold") +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", limits = c(0, 90), guide = "none") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~ modelo) +
  labs(x = "Classe real", y = "Classe predita") +
  theme_minimal() + theme(strip.text = element_text(face = "bold", size = 11))
ggsave(file.path(pasta_fig, "figA_matriz_confusao.png"), p_cm, width = 8, height = 3.6, dpi = 300)

# ------------------------------------------------------------------

# ------------------------------------------------------------------
# FIGURA W - Heatmap e boxplot do AIM2 no conjunto de teste
# (só por grupo; heatmap com ~1/3 da altura do boxplot)
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)

pasta_fig <- "R/meta/figuras_TCC_A"
dir.create(pasta_fig, showWarnings = FALSE, recursive = TRUE)

set.seed(1)
genes_fig <- c("AIM2", "FCGR1A")                      # ordem: de cima para baixo
faltando  <- setdiff(genes_fig, colnames(X_teste))
if (length(faltando) > 0) stop("Genes ausentes em X_teste: ", paste(faltando, collapse = ", "))

rotulo_grupo <- function(y) factor(ifelse(y %in% c("Active.tuberculosis", "Active tuberculosis"),
                                          "TB ativa", "ILTB"), levels = c("TB ativa", "ILTB"))
cores_grupo  <- c("TB ativa" = "#E74C3C", "ILTB" = "#2ECC71")

# tabela larga: uma coluna por gene + z-score de cada gene (no conjunto de teste)
df_teste <- data.frame(amostra = rownames(X_teste), grupo = rotulo_grupo(y_teste),
                       X_teste[, genes_fig, drop = FALSE], check.names = FALSE)
for (g in genes_fig) df_teste[[paste0("z_", g)]] <- as.numeric(scale(df_teste[[g]]))

# mesma ordem de amostras nas duas faixas: grupo e, dentro dele, média dos z-scores
df_teste <- df_teste %>%
  mutate(z_medio = rowMeans(across(starts_with("z_")))) %>%
  arrange(grupo, desc(z_medio)) %>%
  mutate(pos = row_number())
n_por_grupo <- table(df_teste$grupo)
n_tb <- n_por_grupo[["TB ativa"]]

# formato longo para heatmap e boxplot
df_long <- df_teste %>%
  select(amostra, grupo, pos, all_of(genes_fig)) %>%
  pivot_longer(all_of(genes_fig), names_to = "gene", values_to = "expressao") %>%
  group_by(gene) %>%
  mutate(z = pmax(pmin(as.numeric(scale(expressao)), 2.5), -2.5)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, levels = genes_fig))

p_grupo <- ggplot(df_teste, aes(pos, 1, fill = grupo)) +
  geom_raster() +
  scale_fill_manual(values = cores_grupo, name = NULL,
                    labels = paste0(names(n_por_grupo), " (n = ", n_por_grupo, ")")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Grupo") + theme_void() +
  theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 9,
                                    margin = ggplot2::margin(r = 5)),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"))

# uma faixa por gene, uma embaixo da outra (o gene vira o rótulo do eixo y)
p_heat <- ggplot(df_long, aes(pos, gene, fill = z)) +
  geom_raster() +
  geom_hline(yintercept = 1.5, colour = "white", linewidth = 1.2) +
  geom_vline(xintercept = n_tb + 0.5, colour = "white", linewidth = 0.8) +
  scale_fill_viridis_c(option = "magma", limits = c(-2.5, 2.5), name = "Expressão (z-score)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(genes_fig), expand = c(0, 0)) +
  labs(x = paste0("Amostras do conjunto de teste (n = ", nrow(df_teste), ")"), y = NULL) +
  theme_void() +
  theme(axis.text.y  = element_text(size = 9, hjust = 1, face = "bold",
                                    margin = ggplot2::margin(r = 5)),
        axis.title.x = element_text(size = 9, margin = ggplot2::margin(t = 4)),
        legend.position = "bottom", legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"), legend.title = element_text(size = 8, vjust = 0.8))

# boxplots lado a lado, um painel por gene (cada um com o próprio eixo y)
p_box <- ggplot(df_long, aes(grupo, expressao, fill = grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.5) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = cores_grupo) +
  labs(x = NULL, y = "Expressão (log2)") +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))

# espaçadores acima e abaixo centralizam o heatmap na altura dos boxplots
lado_esq <- plot_spacer() / p_grupo / p_heat / plot_spacer() +
  plot_layout(heights = c(0.45, 0.13, 1.6, 0.75))
p_w <- (lado_esq | p_box) + plot_layout(widths = c(4, 2))

print(p_w)
ggsave(file.path(pasta_fig, "figA_heatmap_boxplot_AIM2_FCGR1A.png"), p_w,
       width = 12, height = 4.8, dpi = 300)


# ------------------------------------------------------------------
# FIGURA V - NES dos módulos do CEMiTool (GSEA)
# ------------------------------------------------------------------
descricao_modulos <- c(M1 = "neutrófilos / inflamação", M2 = "linfócitos T e B",
                       M6 = "células citotóxicas")
descricao_modulos[modulo_fcgr1a] <- "interferon (FCGR1A)"

nes_df <- cem@enrichment$nes %>%
  dplyr::select(pathway, NES = `Active tuberculosis`) %>%
  dplyr::left_join(cem@enrichment$padj %>%
                     dplyr::select(pathway, FDR = `Active tuberculosis`),
                   by = "pathway") %>%
  dplyr::mutate(rotulo = ifelse(pathway %in% names(descricao_modulos),
                                paste0(pathway, " – ", descricao_modulos[pathway]), pathway),
                rotulo = reorder(rotulo, NES),
                direcao = ifelse(NES > 0, "Aumentado na TB ativa", "Suprimido na TB ativa"),
                significativo = FDR < 0.05,
                texto_fdr = ifelse(significativo,
                                   paste0("FDR = ", format(FDR, digits = 2, scientific = TRUE)),
                                   "n.s."))

p_nes <- ggplot(nes_df, aes(NES, rotulo, fill = direcao, alpha = significativo)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0) +
  geom_text(aes(label = texto_fdr, hjust = ifelse(NES > 0, -0.1, 1.1)),
            size = 2.6, colour = "grey30", alpha = 1) +
  scale_fill_manual(values = c("Aumentado na TB ativa" = "#E74C3C",
                               "Suprimido na TB ativa" = "#2166AC"), name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.35)) +
  labs(x = "NES na TB ativa (GSEA)", y = NULL) +
  theme_classic() + theme(legend.position = "bottom")
ggsave(file.path(pasta_fig, "figA_cemitool_nes.png"), p_nes, width = 7, height = 4.5, dpi = 300)

cat("\nFiguras salvas em", normalizePath(pasta_fig), ":\n")
print(list.files(pasta_fig, pattern = "^figA_"))


########################################################################
# 9b. PAINEL ROC POR SÉRIE NO CONJUNTO DE TESTE
#
# Curvas ROC do AIM2 em cada série do teste; sombra = IC 95% da
# sensibilidade (bootstrap). No final, confere as AUCs do painel
# com a Tabela 3 (tabela_por_batch_1).
########################################################################

########### PAINEL ROCS POR BATCH NO CONJUNTO DE TESTE ###############

#### FIGURA - Curvas ROC do AIM2 por série no conjunto de teste (com IC 95%) ####
# Rode DEPOIS do pipeline da Abordagem A, na MESMA sessão do R.
# Usa: modelo_rf_1gene, assinatura_1, X_teste, y_teste, batch_teste, tabela_por_batch_1

library(ggplot2)
library(dplyr)
library(pROC)

# "modelo"    -> probabilidade do Random Forest (mesmos valores da Tabela 3)
# "expressao" -> expressão bruta do AIM2, sem modelo (análise complementar)
preditor_tipo <- "modelo"
n_boot <- 2000            # reamostragens bootstrap para a sombra do IC

pasta_fig <- "R/meta/figuras_TCC_A"
dir.create(pasta_fig, showWarnings = FALSE, recursive = TRUE)

nomes_series <- c("1" = "GSE73408", "2" = "GSE39940", "3" = "GSE39939", "4" = "GSE37250",
                  "5" = "GSE19439", "6" = "GSE19442", "7" = "GSE19444")
virgula <- function(x, n = 3) formatC(x, format = "f", digits = n, decimal.mark = ",")

# TB ativa como classe positiva
nivel_tb   <- grep("Active", levels(y_teste), value = TRUE)
nivel_iltb <- setdiff(levels(y_teste), nivel_tb)

# preditor de cada amostra do teste
if (preditor_tipo == "modelo") {
  pred_teste <- predict(modelo_rf_1gene,
                        newdata = X_teste[, assinatura_1, drop = FALSE],
                        type = "prob")[[nivel_tb]]
} else {
  pred_teste <- X_teste[, assinatura_1[1]]
}

faz_roc <- function(y, p) {
  pROC::roc(response = y, predictor = p,
            levels = c(nivel_iltb, nivel_tb), direction = "<", quiet = TRUE)
}
pontos_roc <- function(r, rotulo) {
  data.frame(painel = rotulo,
             fpr = rev(1 - r$specificities),
             tpr = rev(r$sensitivities))
}
# faixa do IC 95% da sensibilidade em cada especificidade (bootstrap)
faixa_ic <- function(r, rotulo) {
  esp <- seq(0, 1, by = 0.02)
  ci  <- suppressWarnings(pROC::ci.se(r, specificities = esp,
                                      boot.n = n_boot, progress = "none"))
  data.frame(painel = rotulo,
             fpr = 1 - esp,
             inf = ci[, 1],
             sup = ci[, 3])
}

set.seed(42)   # bootstrap reprodutível
curvas <- list(); faixas <- list(); rotulos <- list(); conferencia <- list()

# uma curva por série
for (b in sort(unique(batch_teste))) {
  i <- which(batch_teste == b)
  r <- faz_roc(y_teste[i], pred_teste[i])
  nome <- sprintf("%s\nn = %d (%d TB / %d ILTB)",
                  nomes_series[as.character(b)], length(i),
                  sum(y_teste[i] == nivel_tb), sum(y_teste[i] == nivel_iltb))
  curvas[[nome]]  <- pontos_roc(r, nome)
  faixas[[nome]]  <- faixa_ic(r, nome)
  rotulos[[nome]] <- data.frame(painel = nome,
                                lab = paste0("AUC = ", virgula(as.numeric(pROC::auc(r)))))
  conferencia[[nome]] <- data.frame(batch = b, AUC_painel = as.numeric(pROC::auc(r)))
}

# teste completo, com IC 95% da AUC (DeLong)
r_total  <- faz_roc(y_teste, pred_teste)
ic_total <- pROC::ci.auc(r_total, method = "delong")
nome_tot <- sprintf("Teste completo\nn = %d (%d TB / %d ILTB)",
                    length(y_teste), sum(y_teste == nivel_tb), sum(y_teste == nivel_iltb))
curvas[[nome_tot]]  <- pontos_roc(r_total, nome_tot)
faixas[[nome_tot]]  <- faixa_ic(r_total, nome_tot)
rotulos[[nome_tot]] <- data.frame(painel = nome_tot,
                                  lab = sprintf("AUC = %s\n(IC 95%%: %s–%s)",
                                                virgula(ic_total[2]), virgula(ic_total[1]), virgula(ic_total[3])))

ordem <- names(curvas)
curvas_df  <- dplyr::bind_rows(curvas);  curvas_df$painel  <- factor(curvas_df$painel,  levels = ordem)
faixas_df  <- dplyr::bind_rows(faixas);  faixas_df$painel  <- factor(faixas_df$painel,  levels = ordem)
rotulos_df <- dplyr::bind_rows(rotulos); rotulos_df$painel <- factor(rotulos_df$painel, levels = ordem)
curvas_df$tipo <- ifelse(curvas_df$painel == nome_tot, "total", "serie")
faixas_df$tipo <- ifelse(faixas_df$painel == nome_tot, "total", "serie")

cores <- c(serie = "#1B7837", total = "#2166AC")

p_roc_serie <- ggplot() +
  geom_ribbon(data = faixas_df, aes(x = fpr, ymin = inf, ymax = sup, fill = tipo),
              alpha = 0.18, show.legend = FALSE) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_path(data = curvas_df, aes(fpr, tpr, colour = tipo),
            linewidth = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = cores) +
  scale_fill_manual(values = cores) +
  geom_text(data = rotulos_df, aes(x = 0.97, y = 0.05, label = lab),
            hjust = 1, vjust = 0, size = 3.2, fontface = "bold", lineheight = 0.9) +
  facet_wrap(~ painel, ncol = 4) +
  coord_equal() +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = function(x) virgula(x, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = function(x) virgula(x, 1)) +
  labs(x = "1 – Especificidade", y = "Sensibilidade") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = 8.5, face = "bold"),
        panel.grid.minor = element_blank())

print(p_roc_serie)

arquivo <- ifelse(preditor_tipo == "modelo",
                  "figA_roc_por_serie_teste.png",
                  "figA_roc_por_serie_teste_expressao.png")
ggsave(file.path(pasta_fig, arquivo), p_roc_serie, width = 9, height = 5.6, dpi = 300)

#### CONFERÊNCIA: AUCs do painel x Tabela 3 (tabela_por_batch_1) ####
conf <- dplyr::bind_rows(conferencia) %>%
  dplyr::left_join(tabela_por_batch_1 %>% dplyr::select(batch, AUC_tabela = AUC), by = "batch") %>%
  dplyr::mutate(serie = nomes_series[as.character(batch)],
                diferenca = round(AUC_painel - AUC_tabela, 4))
print(conf)
cat("AUC teste completo:", round(as.numeric(pROC::auc(r_total)), 3),
    "| IC 95%:", round(ic_total[1], 3), "-", round(ic_total[3], 3), "\n")
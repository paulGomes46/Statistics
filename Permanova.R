# =========================================
# 0. Packages
# =========================================
library(vegan)
library(readxl)

# =========================================
# 1. Import Excel
# =========================================

file_path <- "~/PhD/ML/Beeswarm/Significance.xlsx"

data <- read_excel(file_path, sheet = "Limo")

# Vérification rapide
str(data)
head(data)

# =========================================
# 2. Séparer classe et matrice
# =========================================

class <- as.factor(data$Class)

# Retirer colonne Class
X <- data[, colnames(data) != "Class"]

# Vérifier dimensions
dim(X)
length(class)

# =========================================
# 3. (PAS de log2 car déjà log2)
# =========================================

# Scaling recommandé
X_scaled <- scale(X)

# =========================================
# 4. Matrice de distance
# =========================================

dist_matrix <- dist(X_scaled, method = "euclidean")

# =========================================
# 5. PERMANOVA globale
# =========================================

set.seed(123)

permanova <- adonis2(dist_matrix ~ class,
                     permutations = 9999)

print(permanova)

# =========================================
# 6. Test dispersion (OBLIGATOIRE)
# =========================================

disp <- betadisper(dist_matrix, class)

anova(disp)
permutest(disp, permutations = 9999)

# =========================================
# 7. Pairwise PERMANOVA
# =========================================

pairwise_permanova <- function(dist_matrix, group, perm = 9999) {
  
  combn_levels <- combn(levels(group), 2)
  
  results <- data.frame()
  
  for (i in 1:ncol(combn_levels)) {
    
    g1 <- combn_levels[1, i]
    g2 <- combn_levels[2, i]
    
    idx <- group %in% c(g1, g2)
    
    sub_dist <- as.dist(as.matrix(dist_matrix)[idx, idx])
    sub_group <- droplevels(group[idx])
    
    res <- adonis2(sub_dist ~ sub_group, permutations = perm)
    
    results <- rbind(results,
                     data.frame(Comparison = paste(g1, "vs", g2),
                                F = res$F[1],
                                R2 = res$R2[1],
                                p = res$`Pr(>F)`[1]))
  }
  
  return(results)
}

pairwise_results <- pairwise_permanova(dist_matrix, class)

print(pairwise_results)
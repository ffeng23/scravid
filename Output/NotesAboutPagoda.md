notes about scde

pagoda.pathway.wPCA returns

$xv: local normalized $score; of each component relative to sampled PC1 sd
        avar <- pmax(0, (xp$sd^2-mean(z[, 1]^2))/sd(z[, 1]^2))
        xv <- t(xp$scores)
        xv <- xv/apply(xv, 1, sd)*sqrt(avar)
    also there is only values for PC1, but not other PCs, since the randomized PC only for PC1.

$xp: wPCA results, $weight, $rotation, $scores, $var, $sd and $total_variance. $var and $sd contains use-defined number of PCs.
        

$z  : this is not z. this n randomized run of pca based on sampled gene from the matrix with same number of genes as $ngene
        this is sd. Not z values.
$n : number of genes in the pathway

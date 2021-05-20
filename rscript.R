# load libs
library(vroom)
library(dplyr)
library(tidyr)
library(ggplot2)

# path to file
ifile <- "results/body.tmp"

# que caracter se metera en los missing GT
miss_char <- ".|."

# Cuanto missing data en general quieres?
missing_fraction <- 0.01        # 0.01 es el 1 % de las celdas van a contener missing data

# Cuantas muestras quieres que tengan pesimo missing data
# este N de muestras tendran un BAD_SAMPLE % de missing data
bad_samples <- 5
bad_sample_miss_fraction <- 0.5      # 0.5 es el 50% de missing data

# Cuantas variantes quieres que tengan pesimo missing data
# en fraccion
bad_variants <- 0.05                    # 0.05 el 5% de las variantes tendran mal missing data
bad_variants_miss_fraction <- 0.5       # 0.5 es el 50% de missing data

# read data
vcfbody.df <- vroom( file = ifile, col_names = FALSE)

# keep main cols
maincols.df <- vcfbody.df %>%
  select( X1:X9)

# Keep sample data cols
samplecols.df <- vcfbody.df %>%
  select( -(X1:X9) ) 

matrix.tmp <- as.matrix(samplecols.df)

# Calc the number of cells to randomly turn into missing  
n_missing_cells <- round( length(matrix.tmp) * missing_fraction )

# create a vector of N random numbers between 1 and length(matrix.tmp)
missing_coordinates.v <- 1:length(matrix.tmp) %>% sample( n_missing_cells )

# replace missing in the matrix
matrix.tmp[ missing_coordinates.v ] <- miss_char

### Alter a sample(s) ----
# chose N sample cols at random
# create a vector of N random numbers between 1 and ncol(matrix.tmp)
missing_cols.v <- 1:ncol(matrix.tmp) %>% sample( bad_samples )

# a bad loop...
for ( n in 1:length( missing_cols.v ) ) {
  
  coln <- missing_cols.v[n]
  
  #get the data
  coltmp <- matrix.tmp [ , coln]
  
  # Calc the number of cells to randomly turn into missing  
  n_miss <- round( length(coltmp) * bad_sample_miss_fraction )
  
  # create a vector of N random numbers between 1 and length(matrix.tmp)
  miss_coord <- 1:length(coltmp) %>% sample( n_miss )
  
  # replace missing in the matrix
  matrix.tmp[ miss_coord, coln ] <- miss_char
  
}

### Alter a position(s) ----
# chose N variant rows at random

# Calc the number of rows to randomly turn bad
n_bad_rows <- round( nrow(matrix.tmp) * bad_variants )

# create a vector of N random numbers between 1 and ncol(matrix.tmp)
missing_rows.v <- 1:nrow(matrix.tmp) %>% sample( n_bad_rows )

# a bad loop...
for ( n in 1:length( missing_rows.v ) ) {
  
  rown <- missing_rows.v[n]
  
  #get the data
  rowtmp <- matrix.tmp [ rown, ]
  
  # Calc the number of cells to randomly turn into missing  
  n_miss2 <- round( length(rowtmp) * bad_variants_miss_fraction )
  
  # create a vector of N random numbers between 1 and length(matrix.tmp)
  miss_coord2 <- 1:length(rowtmp) %>% sample( n_miss2 )
  
  # replace missing in the matrix
  matrix.tmp[ rown, miss_coord2 ] <- miss_char
  
}

# redo whole vcf
samplecols.df <- as.data.frame( matrix.tmp )

# rebind
final_vcfbody.df <- cbind(maincols.df, samplecols.df)

# remote annotation field
final_vcfbody.df$X8 <- "."


## Save modified data
write.table( x = final_vcfbody.df, file = "results/body_with_missingGT.tmp", append = F, quote = F, sep = "\t", row.names = F, col.names = F)

# create an idcol ----
plotable.df <- final_vcfbody.df %>%
  mutate( varid = paste(X1, X2, X3, X4, X5, sep = "_") ) %>%
  select( varid, 10:ncol(vcfbody.df))

# a Heatmap maybe...
long_data.df <- plotable.df %>%
  pivot_longer( cols = -varid, names_to = "Sample", values_to = "GT") %>%
  mutate( tile = ifelse( test = GT == ".|.", yes = "missing", no = "gt") )

# a heatmap
heato.p <- ggplot( data = long_data.df, mapping = aes( x = Sample,
                                                       y = varid,
                                                       fill = tile)) +
  geom_tile() +
  labs( title = "missing data structure",
        subtitle = paste("in file: ", ifile),
        x = "Samples",
        y = "Variants",
        fill = "") +
  scale_fill_manual( values = c("missing" = "red4", "gt" = "white") ) +
  theme_linedraw() +
  theme( axis.text = element_blank() )



# save
ggsave(filename = "elheat.png", plot = heato.p, device = "png",
       width = 7, height = 10, units = "in", dpi = 300)

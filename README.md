# eccStructONT
The repository contains scripts for analysis of Nanopore sequencing of eccDNA.

It has the following useful scripts:

# 1. Jupyter notebook 'eccDraw_HS1_vs_mergedZA.ipynb'
It can be used for calculation of enrichment of eccDNA reads in genome windows of two samples (e.g. heat-stress and ZA sample). 
The variables requred to be assigned before the analysis are in the section "VARIABLES (CHANGE THEM!)". 

# 2. Script 'eccDNA_struct.py'
This script is used to get matrix that can be passed to the heatmap visualization by ComplexHeatmap R package. Example code is below:

```
from eccDNApac.eccDNA_struct import *

fastq = './merged_rep2_eccDNA_stress_ddm1.fastq' 
genome = './GCF_000001735.4_TAIR10.1_genomic.fna'
prefix = 'ddm1_HS2'
TE_coordinate =   'NC_003076.8:4,208,084..4,213,085' #ONSEN3 #EVD: 'NC_003076.8:5,629,978..5,635,310'
output_folder = './structure_output_files'
collected_reads = f'{output_folder}/collected_TE_reads' 
getMatrixAndDrawHist(genome, fastq, TE_coordinate, prefix, output_folder, runTH=True, 
                     window = 50, min_occurence_for_path_to_include_inHM = 1)

```

To draw heatmap with eccDNA structures the follwing code example can be used (in R):

```
ddm1_evd_dt = fread('./structure_output_files/MATRIX_NC_003076.8:4208082..4213083_WT_HS2_ONSEN3_matrix.tab', header = T)

ddm1_evd_dt = ddm1_evd_dt[apply(ddm1_evd_dt[,-c(1)], 1, max, na.rm=TRUE) > 3]
rn = ddm1_evd_dt$V1
ddm1_evd_dt = ddm1_evd_dt[,-c(1)]

max <- apply(ddm1_evd_dt, 1, max, na.rm=TRUE)
mx_count = as.matrix(ddm1_evd_dt)
rownames(mx_count) = rn


png("heatmap_ONSEN3_WT_HS.tiff",width=10,height=12,units="in",res=1200)

white.line <- function(j, i, x, y, w, h, fill) {  }

Heatmap(mx_count, 
        row_title_gp = gpar(fontsize =  10),
        row_names_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 12),
        #column_names_gp = gpar(fontsize = 5),
        cluster_columns = F,
        show_row_names=T,
        show_column_names=F,
        row_names_side="left",
        right_annotation = rowAnnotation(count = anno_barplot(max)),
        heatmap_legend_param=list(color_bar="continuous", 
                                  legend_direction="horizontal", legend_width=unit(3,"cm"),
                                  title_position="topcenter", title_gp=gpar(fontsize=6),
                                  labels_gp = gpar(fontsize = 9)
        ),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (mx_count[i, j] == 0){
        grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "white", fill = "white"))
                # grid.text('=', x, y, gp = gpar(col='grey'))
                }
            grid.lines(x = c(x - width/2, x + width / 2), y = c(y + height / 2, y + height / 2), gp = gpar(col = 'grey', lwd = 2))
            grid.lines(x = c(x - width/2, x + width / 2), y = c(y - height / 2, y - height / 2), gp = gpar(col = 'grey', lwd = 2))

            },#cluster_rows = row_dend, 
        show_row_dend = T, col = wes_palette("Zissou1", 100, type = "continuous"), km=6)
dev.off()

```

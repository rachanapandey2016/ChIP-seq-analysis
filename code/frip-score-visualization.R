
# Load required package
library(ggplot2)
library(readxl)

setwd("C:/PhD_Courses/Fourth Semester/GCD 8141 Computational genomics/Project-2/visualization")

frip_data <- read_excel("C:/PhD_Courses/Fourth Semester/GCD 8141 Computational genomics/Project-2/visualization/Frip_scores.xlsx", col_names=FALSE)

# Assign column names
colnames(frip_data) <- c("Treatment_Conditions", "FRiP_Score")

# Check the first few rows
head(frip_data)

# Create bar plot
pdf(file= "Frip_score_visualization.pdf", height=7, width=9)
ggplot(frip_data, aes(x=Treatment_Conditions, y=FRiP_Score, fill=Treatment_Conditions)) +
   geom_bar(stat="identity", width=0.6) +
   theme_minimal() +
   labs(title="FRiP Score Comparison across Samples", 
        x="Sample Name", 
        y="FRiP Score") +
   theme(
      axis.text.x = element_text(angle=45, hjust=1, size=12, face="bold"),
      axis.text.y = element_text(size=11, face="bold"),
      axis.title.x = element_text(size=12, face="bold"),
      axis.title.y = element_text(size=12, face="bold"),
      plot.title = element_text(size=12, face="bold", hjust=0.5)
   ) +
   scale_fill_brewer(palette="Set2") +  # Optional color scheme
   ylim(0, 0.5)  # Set y-axis limit
dev.off()

setwd("C:/PhD_Courses/Fourth Semester/GCD 8141 Computational genomics/Project-2/Analysis")
library(ggVennDiagram)

files <- list(
   "Control" = read.table("Control_shRNA_summits.bed"),
   "GATA3" = read.table("GATA3_shRNA_summits.bed"),
   "GATA3_JUN" = read.table("GATA3_shRNA_JUN_OE_summits.bed"),
   "JUN" = read.table("JUN_OE_DOX_summits.bed")
)


ggVennDiagram(files) + theme_minimal()


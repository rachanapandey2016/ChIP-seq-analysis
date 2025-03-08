
# Load required package
library(ggplot2)
library(readxl)

setwd("path_to_your_directory_where_FRiP_score_excel_file_is_located")

frip_data <- read_excel("Frip_scores.xlsx", col_names=FALSE)

# Assign column names (In case your excel file doesnot have column names or else this step is not needed)
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

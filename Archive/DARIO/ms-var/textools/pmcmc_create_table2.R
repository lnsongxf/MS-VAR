rm(list=ls())
source("/if/appl/R/Functions/IFfunctions.r") # loads necessary packages and IF helper functions
library(tidyverse)
library(dplyr)
setwd(getwd())


args <- commandArgs(trailingOnly = TRUE)
print(paste0('Target directory: ' , args[1]))
print(paste0('Model name: ', args[2]))

targetDir <- paste0(args[1])
texFile <- paste0(targetDir, args[2],".tex")
csvInput <- paste0(targetDir,args[2],".csv")
tableTitle <-  gsub("/", " ", args[3], fixed = TRUE)
tableTitle <-  gsub("_", ", ", tableTitle, fixed = TRUE)

pmode <- read.csv(csvInput)


texText=cat(paste0("\\documentclass{article}\n",
                   
        "% Packages\n",
        "\\usepackage{amsmath}\n",
        "\\usepackage{amssymb}\n",
        "\\usepackage{caption}\n",
        "\\usepackage[paperwidth=5in,paperheight=5in,centering,margin=0.1in,nofoot]{geometry}\n",
        "\\usepackage{booktabs}\n",
        "\\newcommand{\\ra}[1]{\\renewcommand{\\arraystretch}{#1}}\n",
        "\\begin{document}\n",
        "\\pagenumbering{gobble}\n",
        "%<*tag>\n",
        "\\begin{table}[ht!]\\centering\\ra{1.3}\n",
        "\\caption{", tableTitle, "}\n",
        "\\vspace*{0.2cm}\n",
        "\\begin{tabular}{c l c c}\n",
        "\\hline \\hline \n",
        "\\textbf{Coefficient} & \\textbf{Description} & \\textbf{Bad Regime} & \\textbf{Good Regime} \\\\\n",
        "\\hline\n",
        paste0("$\\alpha_y(s_t)$ & Constant & ",pmode$c_3_1_sync_2[1]," & ",pmode$c_3_1_sync_1[1]," \\\\\n"),
        paste0("& & [",pmode$c_3_1_sync_2[2],",",pmode$c_3_1_sync_2[3],"] & [",pmode$c_3_1_sync_1[2],",",pmode$c_3_1_sync_1[3],"] \\\\\n"),
        paste0("$\\beta_y(s_t)$ & Financial Factor & ",-pmode$a0_3_1_sync_2[1]," & ",-pmode$a0_3_1_sync_1[1]," \\\\\n"),
        paste0("& & [",-pmode$a0_3_1_sync_2[3],",",-pmode$a0_3_1_sync_2[2],"] & [",-pmode$a0_3_1_sync_1[3],",",-pmode$a0_3_1_sync_1[2],"] \\\\\n"),
        paste0("$\\gamma_y(s_t)$ & Macro Factor & ",-pmode$a0_3_2_sync_2[1]," & ",-pmode$a0_3_2_sync_1[1]," \\\\\n"),
        paste0("& & [",-pmode$a0_3_2_sync_2[3],",",-pmode$a0_3_2_sync_2[2],"] & [",-pmode$a0_3_2_sync_1[3],",",-pmode$a0_3_2_sync_1[2],"] \\\\\n"),
        paste0("$\\sigma_y(s_t)$ & Shock to GDP & ",pmode$s_3_3_sync_2[1]," & ",pmode$s_3_3_sync_1[1]," \\\\\n"),
        paste0("& & [",pmode$s_3_3_sync_2[2],",",pmode$s_3_3_sync_2[3],"] & [",pmode$s_3_3_sync_1[2],",",pmode$s_3_3_sync_1[3],"] \\\\\n"),
        # Lag GDP
        paste0("$ y_{t-1}$ & Lag of GDP & ",pmode$a1_3_3_sync_2[1]," & ",pmode$a1_3_3[1]," \\\\\n"),
        paste0("& & [",pmode$a1_3_3_sync_2[2],",",pmode$a1_3_3_sync_2[3],"] & [",pmode$a1_3_3_sync_1[2],",",pmode$a1_3_3_sync_1[3],"] \\\\\n"),
        # Lag FF
        paste0("$ f_{t-1}$ & Lag of FF & ",pmode$a1_3_1_sync_2[1]," & ",pmode$a1_3_1_sync_1[1]," \\\\\n"),
        paste0("& & [",pmode$a1_3_1_sync_2[2],",",pmode$a1_3_1_sync_2[3],"] & [",pmode$a1_3_1_sync_1[2],",",pmode$a1_3_1_sync_1[3],"] \\\\\n"),    
        # Lag MF
        paste0("$ m_{t-1}$ & Lag of MF & ",pmode$a1_3_2_sync_2[1]," & ",pmode$a1_3_2_sync_1[1]," \\\\\n"),
        paste0("& & [",pmode$a1_3_2_sync_2[2],",",pmode$a1_3_2_sync_2[3],"] & [",pmode$a1_3_2_sync_1[2],",",pmode$a1_3_2_sync_1[3],"] \\\\\n"),    
        "\\hline\n",
        "\\end{tabular}\n",                   
        "\\end{table}\n",
        "%</tag>\n",
        "\\end{document}"),file=texFile)

systemCallText = paste("pdflatex -output-directory",targetDir, texFile)
system(systemCallText)

# Author ---------------------------
# Patricia Tran
# August 6, 2019

# Purpose ---------------------------
# This script create a N,C, S and other cycle. At the genome level, presence (red) and absence (black) arrows are used to show the presence of pathways. The pathway cutoff is 0.75, but can be changed when running METABOLIC.
# The community level script creates a similar figure but summarizes the coverage data as a percentage of the community, all arrows are black, and rather, we should look at the percentages instead to determine the presence/potential for that pathway.
# We recommend manually looking through "METABOLIC_result.xlsx" to double check the genes and annotations.

# How to use from the command line ---------------------------
# Open R and type:
# Rscript draw_biogeochemical_cycles.R R_input [name of output folder] [TRUE or FALSE]
# Rscript draw_biogeochemical_cycles.R R_input Output TRUE
# or 
# Rscript draw_biogeochemical_cycles.R R_input Output FALSE
# If no summary file was created in earlier step
# Note that no slash is used at end of input and output folders arguments


# In a folder you have a sample of 99 files: 99 are genome names, and 1 is a summary file named Total.R_input.txt
# Loop through the file names in the folder and use as input in the drawNcycle, drawScycle and drawCcycle below.

# Receive arguments from command line ---------------------------
userprefs <- commandArgs(trailingOnly = TRUE)
R_input <- userprefs[1] # Path to folder with all the summary files R_input.txt
plots.folder.path <- userprefs[2] # Name of new directory to make to store things
summary <- userprefs[3] # Value: TRUE or FALSE, whether you have the Total.R_input.txt file or not in the R_input folder.

if (length(userprefs) > 3){
  mirror.location <- userprefs[4]
}else{
  mirror.location <- "https://cran.mtu.edu"
}

# Create path to folder when plots will be stored ---------------------------
biogeochemcycles.plots.folder <- paste(plots.folder.path, "draw_biogeochem_cycles", sep = "/")

library.path <- .libPaths() 

# Create directory to hold plots ---------------------------
make.plot.directory <- function(FolderPath){
  if (!dir.exists(FolderPath)){
    dir.create(FolderPath)
    cat("made folder: ", FolderPath, "\n")
  }
}


# Nitrogen Cycle Figure ---------------------------
## Nitrogen cycle - Individual ---------------------------
drawNcycle.single <- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  ## Open file connection:
  plot.name <- paste(plot.folder,"/", as.character(name.of.genome),".draw_nitrogen_cycle_single.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5,onefile=FALSE)
  
  library(diagram)
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = paste("Nitrogen Cycle:",name.of.genome)) # Add a title
  elpos <- coordinates (c(1, 2, 2, 2, 1)) # Put the coordinates
  
  straightarrow(from = elpos[1, ], to = elpos[3, ], lty = 1, lcol = input[10,2]) #N-S-01:Nitrogen fixation
  straightarrow(from = elpos[3, ], to = elpos[7, ], lty = 1, lcol = input[11,2]) #N-S-02:Ammonia oxidation
  straightarrow(from = elpos[7, ], to = elpos[8, ], lty = 1, lcol = input[12,2]) #N-S-03:Nitrite oxidation
  straightarrow(from = elpos[8, ], to = elpos[6, ], lty = 1, lcol = input[13,2]) #N-S-04:Nitrate reduction
  straightarrow(from = elpos[6, ], to = elpos[4, ], lty = 1, lcol = input[14,2]) #N-S-05:Nitrite reduction
  straightarrow(from = elpos[4, ], to = elpos[2, ], lty = 1, lcol = input[15,2]) #N-S-06:Nitric oxide reduction
  straightarrow(from = elpos[2, ], to = elpos[1, ], lty = 1, lcol = input[16,2]) #N-S-07:Nitrous oxide reduction
  splitarrow(from = elpos[c(3,6), ], to = elpos[1, ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = input[18,2]) #N-S-09:Anammox
  straightarrow(from = elpos[8, ], to = elpos[3, ], lty = 1, lcol = input[17,2])#N-S-08:Nitrite ammonification 
  
  textrect (elpos[1, ], 0.05, 0.05, lab = expression(paste(N['2'],' ',(0))), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression(paste(N['2'],O,' ',(+1))), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression(paste(NH['4'])^'+'* ' '(-3)), cex = 1.5)
  textrect (elpos[4, ], 0.07, 0.05, lab = expression(paste(NO,' ',(+2))), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* ' '(+3)), cex = 1.5)
  textrect (elpos[7, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* ' '(+3)), cex = 1.5)
  textrect (elpos[8, ], 0.07, 0.05, lab = expression(paste(NO['3'])^'-'* ' '(+5)), cex = 1.5)

  # To know where to write the names (mid = c(X,Y), I ploted elpost (plot(elpos) to get approximate coordinates))
  textplain(mid = c(0.7, 0.85), 
            lab = c("Step1: Nitrogen fixation")) 
  textplain(mid =  c(0.88, 0.5), 
            lab = c("Step2: Ammonia oxidation")) 
  textplain(mid = c(0.73, 0.15), 
            lab = c("Step3: Nitrite oxidation")) 
  textplain(mid = c(0.30, 0.15), 
            lab = c("Step4: Nitrate reduction")) 
  textplain(mid = c(0.12, 0.4), 
            lab = c("Step5: Nitrite reduction")) 
  textplain(mid = c(0.12, 0.59), 
            lab = c("Step6: Nitric oxide reduction")) 
  textplain(mid = c(0.25, 0.85), 
            lab = c("Step7: Nitrous oxide reduction")) 
  textplain(mid = c(0.55, 0.65), 
            lab = c("Step9: Anammox")) 
  textplain(mid = c(0.47, 0.35), 
            lab = c("Step8: Nitrite ammonification"))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}

## Nitrogen Cycle Summary Figure: no coloured arrows, but have the Nb. of Genome annd the coverage next to the arrows: ---------------------------
drawNcycle.total <- function(R_input, OutputFolder){
  input.total <- R_input
  plot.folder <- OutputFolder
  
  # Open file connection
  plot.name <- paste(plot.folder,"/draw_nitrogen_cycle_total.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5,onefile=FALSE)
  
  library(diagram)
  openplotmat()
  #par(mar = c(2, 2, 2, 2))
  openplotmat(main = "Nitrogen Cycle: Summary Figure") # Add a tiitl
  elpos <- coordinates (c(1, 2, 2, 2, 1)) # Put the coordinat
  
  straightarrow(from = elpos[1, ], to = elpos[3, ], lty = 1, lcol = 1) #N-S-01:Nitrogen fixation
  straightarrow(from = elpos[3, ], to = elpos[7, ], lty = 1, lcol = 1) #N-S-02:Ammonia oxidation
  straightarrow(from = elpos[7, ], to = elpos[8, ], lty = 1, lcol = 1) #N-S-03:Nitrite oxidation
  straightarrow(from = elpos[8, ], to = elpos[6, ], lty = 1, lcol = 1) #N-S-04:Nitrate reduction
  straightarrow(from = elpos[6, ], to = elpos[4, ], lty = 1, lcol = 1) #N-S-05:Nitrite reduction
  straightarrow(from = elpos[4, ], to = elpos[2, ], lty = 1, lcol = 1) #N-S-06:Nitric oxide reduction
  straightarrow(from = elpos[2, ], to = elpos[1, ], lty = 1, lcol = 1) #N-S-07:Nitrous oxide reduction
  splitarrow(from = elpos[c(3,6), ], to = elpos[1, ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = 1) #N-S-08:Nitrite ammonification
  straightarrow(from = elpos[8, ], to = elpos[3, ], lty = 1, lcol = 1) #N-S-09:Anammox
  
  textrect (elpos[1, ], 0.05, 0.05, lab = expression(paste(N['2'],' ',(0))), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression(paste(N['2'],O,' ',(+1))), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression(paste(NH['4'])^'+'* ' '(-3)), cex = 1.5)
  textrect (elpos[4, ], 0.07, 0.05, lab = expression(paste(NO,' ',(+2))), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* ' '(+3)), cex = 1.5)
  textrect (elpos[7, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* ' '(+3)), cex = 1.5)
  textrect (elpos[8, ], 0.07, 0.05, lab = expression(paste(NO['3'])^'-'* ' '(+5)), cex = 1.5)
  
  # To know where to write the names (mid = c(X,Y), I ploted elpost (plot(elpos) to get approximate coordinates))
  # Make spacing between lines smaller: 
  par(lheight=0.05)
  
  textplain(mid = c(0.7, 0.85), 
            lab = c("Step1: Nitrogen fixation",
                    paste("Genomes:",input.total$Nb.Genome[10]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[10],"%"))) 
  textplain(mid = c(0.88, 0.5), 
            lab = c("Step2: Ammonia oxidation",
                    paste("Genomes:",input.total$Nb.Genome[11]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[11],"%"))) 
  textplain(mid = c(0.73, 0.15), 
            lab = c("Step3: Nitrite oxidation",
                    paste("Genomes:",input.total$Nb.Genome[12]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[12],"%"))) 
  textplain(mid = c(0.30, 0.15), 
            lab = c("Step4: Nitrate reduction",
                    paste("Genomes:",input.total$Nb.Genome[13]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[13],"%"))) 
  textplain(mid = c(0.12, 0.40), 
            lab = c("Step5: Nitrite reduction",
                    paste("Genomes:",input.total$Nb.Genome[14]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[14],"%"))) 
  textplain(mid = c(0.12, 0.59), 
            lab = c("Step6: Nitric oxide reduction",
                    paste("Genomes:",input.total$Nb.Genome[15]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[15],"%"))) 
  textplain(mid = c(0.28, 0.85), 
            lab = c("Step7: Nitrous oxide reduction",
                    paste("Genomes:",input.total$Nb.Genome[16]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[16],"%"))) 
  textplain(mid = c(0.55, 0.65), 
            lab = c("Step9: Anammox",
                    paste("Genomes:",input.total$Nb.Genome[18]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[18],"%"))) 
  textplain(mid = c(0.47, 0.35), 
            lab = c("Step8: Nitrite ammonification",
                    paste("Genomes:",input.total$Nb.Genome[17]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[17],"%")))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


## Sulfur Cycle - Individual ---------------------------
drawScycle.single <- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/", as.character(name.of.genome),".draw_sulfur_cycle_single.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
  
  library(diagram)
  openplotmat()
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = paste("Sulfur Cycle:", name.of.genome)) # Add a title
  
  elpos <- coordinates (c(1, 3, 3, 3, 1)) # Put the coordinate
  elpos
  straightarrow(from = elpos[1, ], to = elpos[4, ], lty = 1, lcol = input[24,2]) #S-S-01:Sulfide oxidation
  straightarrow(from = elpos[4, ], to = elpos[1, ], lty = 1, lcol = input[25,2]) #S-S-02:Sulfur reduction
  straightarrow(from = elpos[4, ], to = elpos[10, ], lty = 1, lcol = input[26,2]) #S-S-03:Sulfur oxidation
  straightarrow(from = elpos[10, ], to = elpos[11, ], lty = 1, lcol = input[27,2]) #S-S-04:Sulfite oxidation
  straightarrow(from = elpos[11, ], to = elpos[5, ], lty = 1, lcol = input[28,2]) #S-S-05:Sulfate reduction
  straightarrow(from = elpos[5, ], to = elpos[1, ], lty = 1, lcol = input[29,2]) #S-S-06:Sulfite reduction
  straightarrow(from = elpos[6, ], to = elpos[11, ], lty = 1, lcol = input[30,2]) #S-S-07:Thiosulfate oxidation
  splitarrow(from = elpos[6, ], to = elpos[c(1,10), ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = input[31,2]) #S-S-08:Thiosulfate disproportionation
  splitarrow(from = elpos[6, ], to = elpos[c(4,11), ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = input[32,2]) #S-S-09:Thiosulfate disproportionation 2
  
  #https://stackoverflow.com/questions/17083362/colorize-parts-of-the-title-in-a-plot
  
  textrect (elpos[1, ], 0.07, 0.05, lab = expression(paste(H['2'],S,' ',(-2))), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(paste(S,' ',(0))), cex = 1.5)
  textrect (elpos[5, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'*' '(+4)), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(S['2']*O['3'])^'2-'* ' '(+2)), cex = 1.5)
  textrect (elpos[10, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'*' '(+4)), cex = 1.5)
  textrect (elpos[11, ], 0.07, 0.05, lab = expression(paste(SO['4'])^'2-'*' '(+6)), cex = 1.5)
  
  par(lheight=0.08)
  
  textplain(mid = c(0.8, 0.85), 
            lab = c("Step1: Sulfide oxidation"))
  textplain(mid = c(0.68, 0.69), 
            lab = c("Step2: Sulfur reduction"))
  textplain(mid = c(0.93, 0.5), 
            lab = c("Step3: Sulfur oxidation"))
  textplain(mid = c(0.75, 0.15), 
            lab = c("Step4: Sulfite oxidation"))
  textplain(mid = c(0.25, 0.23), 
            lab = c("Step5: Sulfate reduction"))
  textplain(mid = c(0.23, 0.75), 
            lab = c("Step6: Sulfite reduction"))
  textplain(mid = c(0.605, 0.32), 
            lab = c("Step7: Thiosulfate oxidation"))
  textplain(mid = c(0.45, 0.61), 
            lab = c("Step8: Thiosulfate disproportionation 1"))
  textplain(mid = c(0.68, 0.40), 
            lab = c("Step9: Thiosulfate disproportionation 2"))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


## Sulfur Cycle - Total ---------------------------
drawScycle.total <- function(R_input, OutputFolder){
  input.total <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection
  plot.name <- paste(plot.folder, "/draw_sulfur_cycle_total.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
  
  library(diagram)
  openplotmat()
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = "Sulfur Cycle : Summary Figure") # Add a tiitle
  
  elpos <- coordinates (c(1, 3, 3, 3, 1)) # Put the coordinate
  
  straightarrow(from = elpos[1, ], to = elpos[4, ], lty = 1, lcol = 1) #S-S-01:Sulfide oxidation
  straightarrow(from = elpos[4, ], to = elpos[1, ], lty = 1, lcol = 1) #S-S-02:Sulfur reduction
  straightarrow(from = elpos[4, ], to = elpos[10, ], lty = 1, lcol = 1) #S-S-03:Sulfur oxidation
  straightarrow(from = elpos[10, ], to = elpos[11, ], lty = 1, lcol = 1) #S-S-04:Sulfite oxidation
  straightarrow(from = elpos[11, ], to = elpos[5, ], lty = 1, lcol = 1) #S-S-05:Sulfate reduction
  straightarrow(from = elpos[5, ], to = elpos[1, ], lty = 1, lcol = 1) #S-S-06:Sulfite reduction
  straightarrow(from = elpos[6, ], to = elpos[11, ], lty = 1, lcol = 1) #S-S-07:Thiosulfate oxidation
  splitarrow(from = elpos[6, ], to = elpos[c(1,10), ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = 1) #S-S-08:Thiosulfate disproportionation 1
  splitarrow(from = elpos[6, ], to = elpos[c(4,11), ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = 1) #S-S-09:Thiosulfate disproportionation 2
  
  #https://stackoverflow.com/questions/17083362/colorize-parts-of-the-title-in-a-plot
  
  textrect (elpos[1, ], 0.07, 0.05, lab = expression(paste(H['2'],S,' ',(-2))), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(paste(S,' ',(0))), cex = 1.5)
  textrect (elpos[5, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'*' '(+4)), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(S['2']*O['3'])^'2-'* ' '(+2)), cex = 1.5)
  textrect (elpos[10, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'*' '(+4)), cex = 1.5)
  textrect (elpos[11, ], 0.07, 0.05, lab = expression(paste(SO['4'])^'2-'*' '(+6)), cex = 1.5)
  

  par(lheight=0.01)
  
  textplain(mid = c(0.8, 0.85), 
            lab = c("Step1: Sulfide oxidation",
                    paste("Genomes:",input.total$Nb.Genome[24]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[24],"%"))) 
  textplain(mid = c(0.68, 0.68), 
            lab = c("Step2: Sulfur reduction",
                    paste("Genomes:",input.total$Nb.Genome[25]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[25],"%"))) 
  textplain(mid = c(0.93, 0.5), 
            lab = c("Step3: Sulfur",
                    "oxidation",
                    paste("Genomes:",input.total$Nb.Genome[26]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[26],"%")))
  textplain(mid = c(0.75, 0.12), 
            lab = c("Step4: Sulfite oxidation",
                    paste("Genomes:",input.total$Nb.Genome[27]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[27],"%")))
  textplain(mid = c(0.25, 0.20), 
            lab = c("Step5: Sulfate reduction",
                    paste("Genomes:",input.total$Nb.Genome[28]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[28],"%"))) 
  textplain(mid = c(0.23, 0.75), 
            lab = c("Step6: Sulfite reduction",
                    paste("Genomes:",input.total$Nb.Genome[29]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[29],"%")))
  textplain(mid = c(0.41, 0.35), 
            lab = c("Step7: Thiosulfate oxidation",
                    paste("Genomes:",input.total$Nb.Genome[30]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[30],"%")))
  textplain(mid = c(0.45, 0.62), 
            lab = c("Step8: Thiosulfate disproportionation 1",
                    paste("Genomes:",input.total$Nb.Genome[31]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[31],"%")))
  textplain(mid = c(0.68, 0.40), 
            lab = c("Step9: Thiosulfate disproportionation 2",
                    paste("Genomes:",input.total$Nb.Genome[32]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[32],"%")))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
  
}


## Carbon Cycle - Individual ---------------------------
drawCcycle.single <- function(R_input, OutputFolder){
  
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/", as.character(name.of.genome),".draw_carbon_cycle_single.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
  
  library(diagram)
  openplotmat()
  par(mar = c(1, 1, 1, 1))
  openplotmat(main = paste("Carbon Cycle:", name.of.genome)) # Add a title
  elpos <- coordinates (c(1, 2, 1, 2, 1, 1)) # Put the coordinate
  elpos
  curvedarrow(from = elpos[1, ], to = elpos[8, ], curve = -0.5, lty = 1, lcol = input[1,2]) #C-S-01:Organic carbon oxidation
  curvedarrow(from = elpos[8, ], to = elpos[1, ], curve = -0.5, lty = 1, lcol = input[2,2]) #C-S-02:Carbon fixation
  curvedarrow(from = elpos[3, ], to = elpos[8, ], curve = -0.2, lty = 1, lcol = input[3,2]) #C-S-03:Ethanol oxidation
  curvedarrow(from = elpos[2, ], to = elpos[8, ], curve = 0.2, lty = 1, lcol = input[4,2]) #C-S-04:Acetate oxidation
  straightarrow(from = elpos[1, ], to = elpos[4, ], lty = 1, lcol = input[5,2]) #C-S-05:Hydrogen generation
  straightarrow(from = elpos[2, ], to = elpos[3, ], lty = 1, lcol = input[6,2]) #C-S-06:Fermentation
  straightarrow(from = elpos[3, ], to = elpos[2, ], lty = 1, lcol = input[6,2]) #C-S-06:Fermentation
  straightarrow(from = elpos[4, ], to = elpos[7, ], lty = 1, lcol = input[7,2]) #C-S-07:Methanogenesis
  curvedarrow(from = elpos[2, ], to = elpos[7, ], curve = 0.2, lty = 1, lcol = input[7,2]) #C-S-07:Methanogenesis
  curvedarrow(from = elpos[7, ], to = elpos[8, ], curve = 0.1, lty = 1, lcol = input[8,2]) #C-S-08:Methanotrophy
  curvedarrow(from = elpos[8, ], to = elpos[7, ], curve = 0.1, lty = 1, lcol = input[7,2]) #C-S-07:Methanogenesis
  straightarrow(from = elpos[4, ], to = elpos[6, ], lty = 1, lcol = input[9,2]) #C-S-09:Hydrogen oxidation
  
  textrect (elpos[1, ], 0.12, 0.05, lab = expression("Organic carbon"), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression("Acetate"), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression("Ethanol"), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(paste(H['2'])), cex = 1.5)
  textrect (elpos[6, ], 0.05, 0.05, lab = expression(paste(H['2'])*"O"), cex = 1.5)
  textrect (elpos[7, ], 0.05, 0.05, lab = expression(paste(CH['4'])), cex = 1.5)
  textrect (elpos[8, ], 0.05, 0.05, lab = expression(paste(CO['2'])), cex = 1.5)
  
  textplain(mid = c(0.1, 0.7), lab = "Step2: Carbon fixation")
  textplain(mid = c(0.85, 0.9), lab = "Step1: Organic carbon oxidation")
  textplain(mid = c(0.9, 0.55), lab ="Step3: Ethanol oxidation")
  textplain(mid = c(0.2, 0.4), lab = "Step4: Acetate oxidation")
  textplain(mid = c(0.61, 0.68), lab = "Step5: Hydrogen generation")
  textplain(mid = c(0.42, 0.78), lab = "Step6: Fermentation")
  textplain(mid = c(0.4, 0.4), lab = "Step7: Methanogenesis")
  textplain(mid = c(0.39, 0.17), lab = "Step8: Methanotrophy")
  textplain(mid = c(0.68, 0.55), lab = "Step9: Hydrogen oxidation")
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


## Carbon cycle - Total ---------------------------
drawCcycle.total <- function(R_input, OutputFolder){
  
  input.total <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/draw_carbon_cycle_total.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
  
  library(diagram)
  openplotmat()
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = "Carbon Cycle") # Add a title
  elpos <- coordinates (c(1, 2, 1, 2, 1, 1)) # Put the coordinate  mx = , my = -0.1)
  elpos
  curvedarrow(from = elpos[1, ], to = elpos[8, ], curve = -0.5, lty = 1, lcol = 1) #C-S-01:Organic carbon oxidation
  curvedarrow(from = elpos[8, ], to = elpos[1, ], curve = -0.5, lty = 1, lcol = 1) #C-S-02:Carbon fixation
  curvedarrow(from = elpos[3, ], to = elpos[8, ], curve = -0.2, lty = 1, lcol = 1) #C-S-03:Ethanol oxidation
  curvedarrow(from = elpos[2, ], to = elpos[8, ], curve = 0.2, lty = 1, lcol = 1) #C-S-04:Acetate oxidation
  straightarrow(from = elpos[1, ], to = elpos[4, ], lty = 1, lcol = 1) #C-S-05:Hydrogen generation
  straightarrow(from = elpos[2, ], to = elpos[3, ], lty = 1, lcol = 1) #C-S-06:Fermentation
  straightarrow(from = elpos[3, ], to = elpos[2, ], lty = 1, lcol = 1) #C-S-06:Fermentation
  straightarrow(from = elpos[4, ], to = elpos[7, ], lty = 1, lcol = 1) #C-S-07:Methanogenesis
  curvedarrow(from = elpos[2, ], to = elpos[7, ], curve = 0.2, lty = 1, lcol = 1) #C-S-07:Methanogenesis
  curvedarrow(from = elpos[7, ], to = elpos[8, ], curve = 0.1, lty = 1, lcol = 1) #C-S-08:Methanotrophy
  curvedarrow(from = elpos[8, ], to = elpos[7, ], curve = 0.1, lty = 1, lcol = 1) #C-S-07:Methanogenesis
  straightarrow(from = elpos[4, ], to = elpos[6, ], lty = 1, lcol = 1) #C-S-09:Hydrogen oxidation
  
  textrect (elpos[1, ], 0.16, 0.05, lab = expression("Organic Carbon"), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression("Acetate"), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression("Ethanol"), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(paste(H['2'])), cex = 1.5)
  textrect (elpos[6, ], 0.05, 0.05, lab = expression(paste(H['2'])*"O"), cex = 1.5)
  textrect (elpos[7, ], 0.05, 0.05, lab = expression(paste(CH['4'])), cex = 1.5)
  textrect (elpos[8, ], 0.05, 0.05, lab = expression(paste(CO['2'])), cex = 1.5)
  
  textplain(mid = c(0.1, 0.80), lab = c("Step2: Carbon fixation",
            paste("Genomes:",input.total$Nb.Genome[2]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[2],"%")))
  
  textplain(mid = c(0.85, 0.9), lab = c("Step1: Organic carbon oxidation",
            paste("Genomes:",input.total$Nb.Genome[1]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[1],"%")))
  textplain(mid = c(0.9, 0.55), lab =c("Step3: Ethanol oxidation",
            paste("Genomes:",input.total$Nb.Genome[3]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[3],"%")))
  textplain(mid = c(0.2, 0.4), lab = c("Step4: Acetate oxidation",
            paste("Genomes:",input.total$Nb.Genome[4]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[4],"%")))
  textplain(mid = c(0.61, 0.68), lab = c("Step5: Hydrogen generation",
            paste("Genomes:",input.total$Nb.Genome[5]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[5],"%")))
  textplain(mid = c(0.42, 0.78), lab = c("Step6: Fermentation",
            paste("Genomes:",input.total$Nb.Genome[6]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[6],"%")))
  textplain(mid = c(0.4, 0.4), lab = c("Step7: Methanogenesis",
            paste("Genomes:",input.total$Nb.Genome[7]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[7],"%")))
  textplain(mid = c(0.39, 0.17), lab = c("Step8: Methanotrophy",
            paste("Genomes:",input.total$Nb.Genome[8]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[8],"%")))
  textplain(mid = c(0.68, 0.55), lab = c("Step9: Hydrogen oxidation",
            paste("Genomes:",input.total$Nb.Genome[9]),
            paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[9],"%")))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


## Other cycles - Individual ---------------------------
drawOthercycles.single<- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/", as.character(name.of.genome),".draw_other_cycle_single.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
  
  library(diagram)
  openplotmat()
  par(mar = c(1, 1, 1, 1))
  
  openplotmat(main = paste("Other cycles",name.of.genome)) # Add a title
  elpos <- coordinates (c(5, 5, 2, 2)) # Put the coordinate
  elpos
  curvedarrow(from = elpos[2, ], to = elpos[7, ], curve = 0.1, lty = 1, lcol = input[19,2]) #O-S-01:Iron reduction
  curvedarrow(from = elpos[7, ], to = elpos[2, ], curve = 0.1, lty = 1, lcol = input[20,2]) #O-S-02:Iron oxidation
  curvedarrow(from = elpos[11, ], to = elpos[13, ], curve = 0.1, lty = 1, lcol = input[21,2]) #O-S-03:Arsenate reduction
  curvedarrow(from = elpos[13, ], to = elpos[11, ], curve = 0.1, lty = 1, lcol = input[22,2]) #O-S-04:Arsenite oxidation
  straightarrow(from = elpos[12, ], to = elpos[14, ], lty = 1, lcol = input[23,2]) #O-S-05:Selenate reduction

  textrect (elpos[2, ], 0.07, 0.05, lab = expression(Fe^'3+'), cex = 1.5)
  textrect (elpos[7, ], 0.07, 0.05, lab = expression(Fe^'2+'), cex = 1.5)
  textrect (elpos[11, ], 0.07, 0.05, lab = expression(As^'5+'), cex = 1.5)
  textrect (elpos[13, ], 0.07, 0.05, lab = expression(As^'3+'), cex = 1.5)
  textrect (elpos[12, ], 0.10, 0.05, lab = expression(paste(SeO['4'])^'2-'), cex=1.5)
  textrect (elpos[14, ], 0.07, 0.05, lab = expression(Se^'0'), cex=1.5)
  
  par(lheight=0.01)
  
  textplain(mid = c(0.15, 0.75), lab = "Step1: Iron reduction")
  textplain(mid = c(0.43, 0.75), lab ="Step2: Iron oxidation")
  textplain(mid = c(0.12, 0.25), lab ="Step3: Arsenate reduction")
  textplain(mid = c(0.4, 0.25), lab = "Step4: Arsenite oxidation")
  textplain(mid = c(0.88, 0.25), lab = "Step5: Selenate reduction")
  dev.off()
  cat("made plot: ", plot.name, "\n")

}  


## Other cycles - Total ---------------------------
drawOthercycles.total<- function(R_input, OutputFolder){
  input.total <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/draw_other_cycle_total.pdf", sep="")
  pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
  
  library(diagram)
  openplotmat()
  par(mar = c(1, 1, 1, 1), lheight=0.5)
  
  openplotmat(main = paste("Other cycles: Summary Figure")) # Add a title
  elpos <- coordinates (c(5, 5, 2, 2)) # Put the coordinate
  elpos
  curvedarrow(from = elpos[2, ], to = elpos[7, ], curve = 0.1, lty = 1, lcol = 1) #O-S-01:Metal reduction
  curvedarrow(from = elpos[7, ], to = elpos[2, ], curve = 0.1, lty = 1, lcol = 1) #O-S-02:Iron oxidation
  curvedarrow(from = elpos[11, ], to = elpos[13, ], curve = 0.1, lty = 1, lcol = 1) #O-S-03:Arsenate reduction
  curvedarrow(from = elpos[13, ], to = elpos[11, ], curve = 0.1, lty = 1, lcol = 1) #O-S-04:Arsenite oxidation
  straightarrow(from = elpos[12, ], to = elpos[14, ], lty = 1, lcol = 1) #C-S-05:Selenate reduction
  
  textrect (elpos[2, ], 0.07, 0.05, lab = expression(Fe^'3+'), cex = 1.5)
  textrect (elpos[7, ], 0.07, 0.05, lab = expression(Fe^'2+'), cex = 1.5)
  textrect (elpos[11, ], 0.07, 0.05, lab = expression(As^'5+'), cex = 1.5)
  textrect (elpos[13, ], 0.07, 0.05, lab = expression(As^'3+'), cex = 1.5)
  textrect (elpos[12, ], 0.10, 0.05, lab = expression(paste(SeO['4'])^'2-'), cex=1.5)
  textrect (elpos[14, ], 0.07, 0.05, lab = expression(Se^'0'), cex=1.5)
  
  textplain(mid = c(0.15, 0.75),
            lheight=0.5,
            lab = c("Step1: Iron reduction",
                    paste("Genomes:",input.total$Nb.Genome[19]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[19],"%")))

  textplain(mid = c(0.43, 0.75),
            lheight=0.5,
            lab = c("Step2: Iron oxidation",
                    paste("Genomes:",input.total$Nb.Genome[20]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[20],"%")))
					
  textplain(mid = c(0.1, 0.25), 
            lab = c("Step3: Arsenate reduction",
                    paste("Genomes:",input.total$Nb.Genome[21]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[21],"%")))
  
  textplain(mid = c(0.4, 0.25), 
            lab = c("Step4: Arsenite oxidation",
                    paste("Genomes:",input.total$Nb.Genome[22]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[22],"%")))
  
  textplain(mid = c(0.90, 0.25), 
            lab = c("Step5: Selenate reduction",
                    paste("Genomes:",input.total$Nb.Genome[23]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[23],"%")))
  
  dev.off()
  cat("made plot: ", plot.name, "\n")
  
}  


##### Use Functions

## Depending on whether TRUE or FALSE is selected, choose which funcctions to run ---------------------------

make.plot.directory(FolderPath = plots.folder.path)

print(biogeochemcycles.plots.folder)
dir.create(biogeochemcycles.plots.folder)

files <- list.files(path=R_input, pattern="*.txt", full.names=TRUE, recursive=FALSE)
file.total <- list.files(path=R_input, pattern="Total.R_input.txt", full.names=TRUE)

# Remove the total file from the "files" list of individual genome:
files <- setdiff(files, file.total)

# files are path to files not the actual files!
print(paste("There are:",length(files), "genomes to process", sep=" "))

for (i in 1:length(files)){
  input <- read.table(files[i], sep="\t")
  name.of.genome <- as.character(files[i])
  
  name.of.genome <- tail(unlist(strsplit(name.of.genome, "/")),n=1)
  name.of.genome <- unlist(strsplit(name.of.genome, ".R_input.txt"))

  print(name.of.genome)
  
  colnames(input) <- c("Reaction","PresenceAbsence")
  input$PresenceAbsence <- input$PresenceAbsence + 1
    
  drawNcycle.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  drawScycle.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  drawCcycle.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  drawOthercycles.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  print("Next genome")
}


## Generating the summary figures (total community) ---------------------------

if (summary==TRUE){
  print(paste("There is:",length(file.total), "total summary file", sep=" "))
  
  # Total
  # in the "Total R input file, the 2nd column is the number of genomes and the 3rd column is the Percentage of them.
  input.total <- read.table(file.total, sep="\t")
  colnames(input.total) <- c("Reaction","Nb.Genome","Genome.Coverage.Percentages")
  # change the genome coverage percentages to actual percentages:
  input.total$Genome.Coverage.Percentages <- input.total$Genome.Coverage.Percentages*100
  # Round to 1 digit:
  input.total$Genome.Coverage.Percentages.Round <- round(input.total$Genome.Coverage.Percentages, digits = 1)
  
  print("Making the summary figures:")
  drawNcycle.total(R_input = input.total, OutputFolder = biogeochemcycles.plots.folder)
  drawScycle.total(R_input = input.total, OutputFolder = biogeochemcycles.plots.folder)
  drawCcycle.total(R_input = input.total, OutputFolder = biogeochemcycles.plots.folder)
  drawOthercycles.total(R_input = input.total, OutputFolder = biogeochemcycles.plots.folder)
}else{
  print("Summary=FALSE, no files to process, not generating summary figures")
}



print("Done! The Biogeochemical Cycle Plots have been created :-)")

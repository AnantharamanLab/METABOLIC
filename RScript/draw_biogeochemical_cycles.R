# Patricia Tran
# July 23. 2019

# METABOLIC DIAGRAM GENERATOR:
# usage from Command line
# Rscript draw_biogeochemical_cycles.R R_input/ Output


# In a folder you have a sample of 99 files: 99 are MAG names, and 1 is a summary file named Total.R_input.txt
# Loop through the file names in the folder and use as input for "R_input" in the drawNcycle, drawScycle and drawCcycle below.

# Receive arguments from command line:
userprefs <- commandArgs(trailingOnly = TRUE)
R_input <- userprefs[1] # Path to folder with all the summary files R_input.txt
plots.folder.path <- userprefs[2] # Name of new directory to make to store things

if (length(userprefs) > 2){
  mirror.location <- userprefs[3]
}else{
  mirror.location <- "https://cran.mtu.edu"
}

biogeochemcycles.plots.folder <- paste(plots.folder.path, "draw_biogeochem_cycles", sep = "/")

library.path <- .libPaths() 

# create directories to hold plots
make.plot.directory <- function(FolderPath){
  if (!dir.exists(FolderPath)){
    dir.create(FolderPath)
    cat("made folder: ", FolderPath, "\n")
  }
}


### NITROGEN CYCLE ### 
# Nitrogen
drawNcycle.single <- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  ## Open file connection:
  plot.name <- paste(plot.folder,"/", as.character(name.of.MAG),"_draw_nitrogen_cycle_single.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = paste("Nitrogen Cycle:",name.of.MAG)) # Add a title
  elpos <- coordinates (c(1, 2, 2, 2, 1), mx = 0.1, my = -0.1) # Put the coordinates
  
  straightarrow(from = elpos[1, ], to = elpos[3, ], lty = 1, lcol = input[11,2]) #N-S-01:Nitrogen fixation
  straightarrow(from = elpos[3, ], to = elpos[7, ], lty = 1, lcol = input[12,2]) #N-S-02:Ammonia oxidation
  straightarrow(from = elpos[7, ], to = elpos[8, ], lty = 1, lcol = input[13,2]) #N-S-03:Nitrite oxidation
  straightarrow(from = elpos[8, ], to = elpos[6, ], lty = 1, lcol = input[14,2]) #N-S-04:Nitrate reduction
  straightarrow(from = elpos[6, ], to = elpos[4, ], lty = 1, lcol = input[15,2]) #N-S-05:Nitrite reduction
  straightarrow(from = elpos[4, ], to = elpos[2, ], lty = 1, lcol = input[16,2]) #N-S-06:Nitric oxide reduction
  straightarrow(from = elpos[2, ], to = elpos[1, ], lty = 1, lcol = input[17,2]) #N-S-07:Nitrous oxide reduction
  splitarrow(from = elpos[c(3,6), ], to = elpos[1, ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = input[18,2]) #N-S-08:Nitrite ammonification
  straightarrow(from = elpos[8, ], to = elpos[3, ], lty = 1, lcol = input[19,2]) #N-S-09:Anammox
  
  textrect (elpos[1, ], 0.05, 0.05, lab = expression(N['2'](0)), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression(paste(N['2'],O,(+1))), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression(paste(NH['4'])^'+'* (-3)), cex = 1.5)
  textrect (elpos[4, ], 0.07, 0.05, lab = expression(NO(+2)), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* (+3)), cex = 1.5)
  textrect (elpos[7, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* (+3)), cex = 1.5)
  textrect (elpos[8, ], 0.07, 0.05, lab = expression(paste(NO['3'])^'-'* (+5)), cex = 1.5)

  # To know where to write the names (mid = c(X,Y), I ploted elpost (plot(elpos) to get approximate coordinates))
  textplain(mid = c(0.8, 0.8), 
            lab = c("N-S-01:Nitrogen fixation")) 
  textplain(mid = c(1, 0.4), 
            lab = c("N-S-02:Ammonia oxidation")) 
  textplain(mid = c(0.75, 0.10), 
            lab = c("N-S-03:Nitrite oxidation")) 
  textplain(mid = c(0.45, 0.10), 
            lab = c("N-S-04:Nitrate reduction")) 
  textplain(mid = c(0.3, 0.3), 
            lab = c("N-S-05:Nitrite reduction")) 
  textplain(mid = c(0.3, 0.5), 
            lab = c("N-S-06:Nitric oxide reduction")) 
  textplain(mid = c(0.45, 0.75), 
            lab = c("N-S-07:Nitrous oxide reduction")) 
  textplain(mid = c(0.6, 0.6), 
            lab = c("N-S-09:Anammox")) 
  textplain(mid = c(0.6, 0.2), 
            lab = c("N-S-08:Nitrite ammonification"))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


# Nitrogen Cycle Summary Figure: no coloured arrows, but have the Nb. of Genome annd the Completeness next to the arrows:
drawNcycle.total <- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  # Open file connection
  plot.name <- paste(plot.folder,"/draw_nitrogen_cycle_total.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  openplotmat()
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = "Nitrogen Cycle: Summary Figure") # Add a tiitl
  elpos <- coordinates (c(1, 2, 2, 2, 1), mx = 0.1, my = -0.1) # Put the coordinat
  
  straightarrow(from = elpos[1, ], to = elpos[3, ], lty = 1, lcol = 1) #N-S-01:Nitrogen fixation
  straightarrow(from = elpos[3, ], to = elpos[7, ], lty = 1, lcol = 1) #N-S-02:Ammonia oxidation
  straightarrow(from = elpos[7, ], to = elpos[8, ], lty = 1, lcol = 1) #N-S-03:Nitrite oxidation
  straightarrow(from = elpos[8, ], to = elpos[6, ], lty = 1, lcol = 1) #N-S-04:Nitrate reduction
  straightarrow(from = elpos[6, ], to = elpos[4, ], lty = 1, lcol = 1) #N-S-05:Nitrite reduction
  straightarrow(from = elpos[4, ], to = elpos[2, ], lty = 1, lcol = 1) #N-S-06:Nitric oxide reduction
  straightarrow(from = elpos[2, ], to = elpos[1, ], lty = 1, lcol = 1) #N-S-07:Nitrous oxide reduction
  splitarrow(from = elpos[c(3,6), ], to = elpos[1, ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = 1) #N-S-08:Nitrite ammonification
  straightarrow(from = elpos[8, ], to = elpos[3, ], lty = 1, lcol = 1) #N-S-09:Anammox
  
  textrect (elpos[1, ], 0.05, 0.05, lab = expression(N['2'](0)), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression(paste(N['2'],O,(+1))), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression(paste(NH['4'])^'+'* (-3)), cex = 1.5)
  textrect (elpos[4, ], 0.07, 0.05, lab = expression(NO(+2)), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* (+3)), cex = 1.5)
  textrect (elpos[7, ], 0.07, 0.05, lab = expression(paste(NO['2'])^'-'* (+3)), cex = 1.5)
  textrect (elpos[8, ], 0.07, 0.05, lab = expression(paste(NO['3'])^'-'* (+5)), cex = 1.5)
  
  # To know where to write the names (mid = c(X,Y), I ploted elpost (plot(elpos) to get approximate coordinates))
  # Make spacing between lines smaller: 
  par(lheight=0.05)
  
  textplain(mid = c(0.8, 0.8), 
            lab = c("N-S-01:Nitrogen fixation",
                    paste("Nb.Genome:",input.total$Nb.Genome[11]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[11],"%"))) 
  textplain(mid = c(0.9, 0.4), 
            lab = c("N-S-02:Ammonia oxidation",
                    paste("Nb.Genome:",input.total$Nb.Genome[12]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[12],"%"))) 
  textplain(mid = c(0.75, 0.10), 
            lab = c("N-S-03:Nitrite oxidation",
                    paste("Nb.Genome:",input.total$Nb.Genome[13]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[13],"%"))) 
  textplain(mid = c(0.45, 0.10), 
            lab = c("N-S-04:Nitrate reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[14]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[15],"%"))) 
  textplain(mid = c(0.2, 0.3), 
            lab = c("N-S-05:Nitrite reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[15]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[5],"%"))) 
  textplain(mid = c(0.2, 0.5), 
            lab = c("N-S-06:Nitric oxide reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[16]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[16],"%"))) 
  textplain(mid = c(0.45, 0.75), 
            lab = c("N-S-07:Nitrous oxide reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[17]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[17],"%"))) 
  textplain(mid = c(0.6, 0.6), 
            lab = c("N-S-09:Anammox",
                    paste("Nb.Genome:",input.total$Nb.Genome[18]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[18],"%"))) 
  textplain(mid = c(0.6, 0.2), 
            lab = c("N-S-08:Nitrite ammonification",
                    paste("Nb.Genome:",input.total$Nb.Genome[19]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[19],"%")))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


##Sulfur cycle##
drawScycle.single <- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/", as.character(name.of.MAG),"draw_sulfur_cycle_single.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  openplotmat()
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = paste("Sulfur Cycle", name.of.MAG)) # Add a tiitle
  
  elpos <- coordinates (c(1, 3, 3, 3, 1), mx = 0.1, my = -0.1) # Put the coordinate
  elpos
  straightarrow(from = elpos[1, ], to = elpos[4, ], lty = 1, lcol = R_input[24,2]) #S-S-01:Sulfide oxidation
  straightarrow(from = elpos[4, ], to = elpos[1, ], lty = 1, lcol = R_input[25,2]) #S-S-02:Sulfur reduction
  straightarrow(from = elpos[4, ], to = elpos[10, ], lty = 1, lcol = R_input[26,2]) #S-S-03:Sulfur oxidation
  straightarrow(from = elpos[10, ], to = elpos[11, ], lty = 1, lcol = R_input[27,2]) #S-S-04:Sulfite oxidation
  straightarrow(from = elpos[11, ], to = elpos[5, ], lty = 1, lcol = R_input[28,2]) #S-S-05:Sulfate reduction
  straightarrow(from = elpos[5, ], to = elpos[1, ], lty = 1, lcol = R_input[29,2]) #S-S-06:Sulfite reduction
  straightarrow(from = elpos[6, ], to = elpos[11, ], lty = 1, lcol = R_input[30,2]) #S-S-07:Thiosulfate oxidation
  splitarrow(from = elpos[6, ], to = elpos[c(1,10), ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = R_input[31,2]) #S-S-08:Thiosulfate disproportionation
  
  #https://stackoverflow.com/questions/17083362/colorize-parts-of-the-title-in-a-plot
  
  textrect (elpos[1, ], 0.07, 0.05, lab = expression(paste(H['2'],S,(-2))), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(S(0)), cex = 1.5)
  textrect (elpos[5, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'* (+4)), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(S['2']*O['3'])^'2-'* (+2)), cex = 1.5)
  textrect (elpos[10, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'* (+4)), cex = 1.5)
  textrect (elpos[11, ], 0.07, 0.05, lab = expression(paste(SO['4'])^'2-'* (+6)), cex = 1.5)
  
  par(lheight=0.08)
  
  textplain(mid = c(0.8, 0.8), 
            lab = c("S-S-01:Sulfide oxidation"))
  textplain(mid = c(0.7, 0.65), 
            lab = c("S-S-02:Sulfur reduction"))
  textplain(mid = c(1, 0.4), 
            lab = c("S-S-03:Sulfur oxidation"))
  textplain(mid = c(0.8, 0.02), 
            lab = c("S-S-04:Sulfite oxidation"))
  textplain(mid = c(0.35, 0.1), 
            lab = c("S-S-05:Sulfate reduction"))
  textplain(mid = c(0.4, 0.7), 
            lab = c("S-S-06:Sulfite reduction"))
  textplain(mid = c(0.70, 0.2), 
            lab = c("S-S-07:Thiosulfate oxidation"))
  textplain(mid = c(0.80, 0.5), 
            lab = c("S-08:Thiosulfate disproportionation"))
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}


drawScycle.total <- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection
  plot.name <- paste(plot.folder, "/draw_sulfur_cycle_total.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  openplotmat()
  par(mar = c(2, 2, 2, 2))
  openplotmat(main = "Sulfur Cycle : Summary Figure") # Add a tiitle
  
  elpos <- coordinates (c(1, 3, 3, 3, 1), mx = 0.1, my = -0.1) # Put the coordinate
  
  straightarrow(from = elpos[1, ], to = elpos[4, ], lty = 1, lcol = 1) #S-S-01:Sulfide oxidation
  straightarrow(from = elpos[4, ], to = elpos[1, ], lty = 1, lcol = 1) #S-S-02:Sulfur reduction
  straightarrow(from = elpos[4, ], to = elpos[10, ], lty = 1, lcol = 1) #S-S-03:Sulfur oxidation
  straightarrow(from = elpos[10, ], to = elpos[11, ], lty = 1, lcol = 1) #S-S-04:Sulfite oxidation
  straightarrow(from = elpos[11, ], to = elpos[5, ], lty = 1, lcol = 1) #S-S-05:Sulfate reduction
  straightarrow(from = elpos[5, ], to = elpos[1, ], lty = 1, lcol = 1) #S-S-06:Sulfite reduction
  straightarrow(from = elpos[6, ], to = elpos[11, ], lty = 1, lcol = 1) #S-S-07:Thiosulfate oxidation
  splitarrow(from = elpos[6, ], to = elpos[c(1,10), ], lty = 1, lwd = 1, dd = 0.7, arr.side = 1:2, lcol = 1) #S-S-08:Thiosulfate disproportionation
  
  #https://stackoverflow.com/questions/17083362/colorize-parts-of-the-title-in-a-plot
  
  textrect (elpos[1, ], 0.07, 0.05, lab = expression(paste(H['2'],S,(-2))), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(S(0)), cex = 1.5)
  textrect (elpos[5, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'* (+4)), cex = 1.5)
  textrect (elpos[6, ], 0.07, 0.05, lab = expression(paste(S['2']*O['3'])^'2-'* (+2)), cex = 1.5)
  textrect (elpos[10, ], 0.07, 0.05, lab = expression(paste(SO['3'])^'2-'* (+4)), cex = 1.5)
  textrect (elpos[11, ], 0.07, 0.05, lab = expression(paste(SO['4'])^'2-'* (+6)), cex = 1.5)
  

  par(lheight=0.01)
  
  textplain(mid = c(0.8, 0.8), 
            lab = c("S-S-01:Sulfide oxidation",
                    paste("Nb.Genome:",input.total$Nb.Genome[11]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[11],"%")),
            cex=0.8) 
  textplain(mid = c(0.7, 0.65), 
            lab = c("S-S-02:Sulfur reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[12]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[12],"%")),
            cex=0.8) 
  textplain(mid = c(1, 0.4), 
            lab = c("S-S-03:Sulfur oxidation",
                    paste("Nb.Genome:",input.total$Nb.Genome[13]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[13],"%")),
            cex=0.8)
  textplain(mid = c(0.8, 0.02), 
            lab = c("S-S-04:Sulfite oxidation",
                    paste("Nb.Genome:",input.total$Nb.Genome[14]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[15],"%")),
            cex=0.8) 
  textplain(mid = c(0.35, 0.1), 
            lab = c("S-S-05:Sulfate reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[15]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[5],"%")),
            cex=0.8) 
  textplain(mid = c(0.4, 0.7), 
            lab = c("S-S-06:Sulfite reduction",
                    paste("Nb.Genome:",input.total$Nb.Genome[17]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[17],"%")),
            cex=0.8) 
  textplain(mid = c(0.70, 0.2), 
            lab = c("S-S-07:Thiosulfate oxidation",
                    paste("Nb.Genome:",input.total$Nb.Genome[18]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[18],"%")),
            cex=0.8) 
  textplain(mid = c(0.80, 0.5), 
            lab = c("S-08:Thiosulfate disproportionation",
                    paste("Nb.Genome:",input.total$Nb.Genome[19]),
                    paste("Coverage:",input.total$Genome.Coverage.Percentages.Round[19],"%")),
            cex=0.8)
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
  
}

# Carbon cycle
drawCcycle.single <- function(R_input, OutputFolder){
  
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/", as.character(name.of.MAG),"draw_carbon_cycle_single.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  openplotmat()
  par(mar = c(1, 1, 1, 1))
  openplotmat(main = paste("Carbon Cycle:", name.of.MAG)) # Add a title
  elpos <- coordinates (c(1, 2, 1, 2, 1, 1), mx = 0.1, my = -0.1) # Put the coordinate
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
  
  textrect (elpos[1, ], 0.12, 0.05, lab = expression("Organic carbon"), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression("Acetate"), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression("Ethanol"), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(paste(H['2'])), cex = 1.5)
  textrect (elpos[6, ], 0.05, 0.05, lab = expression(paste(H['2'])*"O"), cex = 1.5)
  textrect (elpos[7, ], 0.05, 0.05, lab = expression(paste(CH['4'])), cex = 1.5)
  textrect (elpos[8, ], 0.05, 0.05, lab = expression(paste(CO['2'])), cex = 1.5)
  
  #Once the plot is done, export it:
  dev.off()
  cat("made plot: ", plot.name, "\n")
}

drawCcycle.total <- function(R_input, OutputFolder){
  
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/draw_carbon_cycle_total.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  openplotmat()
  par(mar = c(1, 1, 1, 1))
  openplotmat(main = "Carbon Cycle") # Add a title
  elpos <- coordinates (c(1, 2, 1, 2, 1, 1), mx = 0.1, my = -0.1) # Put the coordinate
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
  
  textrect (elpos[1, ], 0.12, 0.05, lab = expression("Organic carbon"), cex = 1.5)
  textrect (elpos[2, ], 0.07, 0.05, lab = expression("Acetate"), cex = 1.5)
  textrect (elpos[3, ], 0.07, 0.05, lab = expression("Ethanol"), cex = 1.5)
  textrect (elpos[4, ], 0.05, 0.05, lab = expression(paste(H['2'])), cex = 1.5)
  textrect (elpos[6, ], 0.05, 0.05, lab = expression(paste(H['2'])*"O"), cex = 1.5)
  textrect (elpos[7, ], 0.05, 0.05, lab = expression(paste(CH['4'])), cex = 1.5)
  textrect (elpos[8, ], 0.05, 0.05, lab = expression(paste(CO['2'])), cex = 1.5)
  
  #Once the plot is done, export it:
  unnecessary.comment <- dev.off()
  cat("made plot: ", plot.name, "\n")
}

# Other cycles
drawOthercycles.single<- function(R_input, OutputFolder){
  input <- R_input
  plot.folder <- OutputFolder
  
  #Open file connection:
  plot.name <- paste(plot.folder, "/", as.character(name.of.MAG),"draw_other_cycle_single.png", sep="")
  png(filename = plot.name, width = 10, height = 10, units = "in", res = 100)
  
  library(diagram)
  openplotmat()
  par(mar = c(1, 1, 1, 1))
  openplotmat(main = paste("Other cycles",name.of.MAG)) # Add a title
  elpos <- coordinates (c(5, 5, 2, 2), mx = 0.1, my = -0.1) # Put the coordinate
  elpos
  curvedarrow(from = elpos[2, ], to = elpos[7, ], curve = 0.1, lty = 1, lcol = 1) #O-S-01:Metal reduction
  curvedarrow(from = elpos[7, ], to = elpos[2, ], curve = 0.1, lty = 1, lcol = 1) #
  straightarrow(from = elpos[3, ], to = elpos[9, ], lty = 1, lcol = 1) #O-S-02:chlorate reduction
  straightarrow(from = elpos[5, ], to = elpos[9, ], lty = 1, lcol = 1) #O-S-02:Perchlorate reduction
  curvedarrow(from = elpos[11, ], to = elpos[13, ], curve = 0.1, lty = 1, lcol = 1) #O-S-03:Arsenate reduction
  curvedarrow(from = elpos[13, ], to = elpos[11, ], curve = 0.1, lty = 1, lcol = 1) #O-S-04:Arsenite oxidation
  straightarrow(from = elpos[12, ], to = elpos[14, ], lty = 1, lcol = 1) #C-S-05:Selenate reduction

  dev.off()
  cat("made plot: ", plot.name, "\n")

  }  


##### Use Functions

make.plot.directory(FolderPath = plots.folder.path)

print(biogeochemcycles.plots.folder)
dir.create(biogeochemcycles.plots.folder)

files <- list.files(path=R_input, pattern="*.txt", full.names=TRUE, recursive=FALSE)
file.total <- list.files(path=R_input, pattern="Total.R_input.txt", full.names=TRUE)
# Remove the total file from the "files" list of individual MAG:
files <- setdiff(files, file.total)

print(files[1])

# files are path to files not the actual files!

print(paste("There are:",length(files), "MAGs to process", sep=" "))
print(paste("There is:",length(file.total), "total summary file", sep=" "))

# Total
# in the "Total R input file, the 2nd column is the number of genomes and the 3rd column is the Percentage of them.
input.total <- read.table("~/Documents/METABOLIC_R_input/R_input_files/Total.R_input.txt", sep="\t")
colnames(input.total) <- c("Reaction","Nb.Genome","Genome.Coverage.Percentages")
# change the genome coverage percentages to actual percentages:
input.total$Genome.Coverage.Percentages <- input.total$Genome.Coverage.Percentages*100
# Round to 1 digit:
input.total$Genome.Coverage.Percentages.Round <- round(input.total$Genome.Coverage.Percentages, digits = 1)

for (i in 1:length(files)){
  input <- read.table(files[i], sep="\t")
  name.of.MAG <- as.character(files[i])
  
  name.of.MAG <- unlist(strsplit(name.of.MAG, "/"))[-1]
  
  print(name.of.MAG)
  
  colnames(input) <- c("Reaction","PresenceAbsence")
  input$PresenceAbsence <- input$PresenceAbsence + 1
    
  drawNcycle.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  drawScycle.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  drawCcycle.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  drawOthercycles.single(R_input = input, OutputFolder = biogeochemcycles.plots.folder)
  print("Next MAG")
}

# Generating summary figures:

print("Making the summary figures:")
drawNcycle.total(R_input = file.total, OutputFolder = biogeochemcycles.plots.folder)
drawScycle.total(R_input = file.total, OutputFolder = biogeochemcycles.plots.folder)
drawCcycle.total(R_input = file.total, OutputFolder = biogeochemcycles.plots.folder)






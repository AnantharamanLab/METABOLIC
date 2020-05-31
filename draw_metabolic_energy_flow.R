# Patricia Tran
# Nov 5, 2019

userprefs <- commandArgs(trailingOnly = TRUE)
Energy_flow_input <- userprefs[1] # Path to to Energy_flow_input.txt
plots.folder.path <- userprefs[2] # Name of new output folder

if (length(userprefs) > 2){
  mirror.location <- userprefs[4]
}else{
  mirror.location <- "https://cran.mtu.edu"
}

energy.plots.folder <- paste(plots.folder.path, "Energy_plot", sep = "/")

library.path <- .libPaths() 

# create directories to hold plots
make.plot.directory <- function(FolderPath){
  if (!dir.exists(FolderPath)){
    dir.create(FolderPath)
    cat("made folder: ", FolderPath, "\n")
  }
}


make.plot.directory(FolderPath = plots.folder.path)
print(energy.plots.folder)
dir.create(energy.plots.folder)

plot.folder <- energy.plots.folder

# Energy flow plot:

energy.flow <- read.table(Energy_flow_input,sep="\t", header=FALSE)
colnames(energy.flow) <- c("Taxa", "Reaction", "Freq")


#unique(energy.flow$Reaction)
# There are 31 levels, and it probably makes sense to split them by general category?


energy.flow$Category <- ifelse(grepl("C-S", energy.flow$Reaction), "Carbon", 
                               ifelse(grepl("N-S", energy.flow$Reaction), "Nitrogen",
                                      ifelse(grepl("S-S", energy.flow$Reaction), "Sulfur",
                                             ifelse(grepl("O-S", energy.flow$Reaction), "Oxygen",""))))

#str(energy.flow)
energy.flow$Category <- as.factor(energy.flow$Category)



library(ggalluvial)
print("Is this in the alluvial format?")
is_alluvia_form(as.data.frame(energy.flow), axes = 1:3, silent = TRUE)

#str(energy.flow)

library(ggthemes)
alluvial.plot <- ggplot(as.data.frame(energy.flow),
       aes(y = Freq, axis1= Taxa, axis2 = Reaction, axis3 = Category)) +
  geom_alluvium(aes(fill = Taxa), width = 1/12)+
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", infer.label = TRUE) +
# label.strata was deprecated so I changed it to infer.label
  scale_x_discrete(limits = c("Reaction", "Taxa","Category"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Reactions and taxa")+
  theme_bw()

plot.name <- paste(plot.folder,"/","network.plot.pdf", sep="")
pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
alluvial.plot
dev.off()

# Patricia Tran
# May 29, 2020

# In order to fix the Perl Excel issue, we are trying to load the csv files and modify it and save it as a Excel spreadsheet.

userprefs <- commandArgs(trailingOnly = TRUE)
Dir_with_TSV_files <- userprefs[1] # Path to the directory that contains all the .TSV files (named METABOLIC_result_worksheetX.tsv)

if (length(userprefs) > 2){
  mirror.location <- userprefs[3]
}else{
  mirror.location <- "https://cran.mtu.edu"
}

# Load all the pages:
page1 <- read.delim(paste0(Dir_with_TSV_files,"METABOLIC_result_worksheet1.tsv"), sep="\t", header=TRUE)
page2 <- read.delim(paste0(Dir_with_TSV_files,"METABOLIC_result_worksheet2.tsv"), sep="\t", header=TRUE)
page3 <- read.delim(paste0(Dir_with_TSV_files,"METABOLIC_result_worksheet3.tsv"), sep="\t", header=TRUE)
page4 <- read.delim(paste0(Dir_with_TSV_files,"METABOLIC_result_worksheet4.tsv"), sep="\t", header=TRUE)
page5 <- read.delim(paste0(Dir_with_TSV_files,"METABOLIC_result_worksheet5.tsv"), sep="\t", header=TRUE)
page6 <- read.delim(paste0(Dir_with_TSV_files,"METABOLIC_result_worksheet6.tsv"), sep="\t", header=TRUE)

#devtools::install_github("colearendt/xlsx@develop")
#library(xlsx)

#install.packages("rJava",type='source')
#library(rJava)

#+++++++++++++++++++++++++++
# xlsx.writeMultipleData
#+++++++++++++++++++++++++++++
# file : the path to the output file
# ... : a list of data to write to the workbook
xlsx.writeMultipleData <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i], row.names=FALSE)
    else write.xlsx(objects[[i]], file, sheetName = objnames[i], row.names = FALSE,
                    append = TRUE)
  }
}


## Now load it again and make soem transformations:
# First bold the header:
# install.packages("openxlsx")
library(openxlsx)

# Make a workboot with all these pages:

wb <- createWorkbook()
addWorksheet(wb, "HMMHitNum")
writeData(wb, sheet = 1, page1)
addWorksheet(wb, "FunctionHit")
writeData(wb, sheet = 2, page2)
addWorksheet(wb, "KEGGModuleHit")
writeData(wb, sheet = 3, page3)
addWorksheet(wb, "KEGGModuleStepHit")
writeData(wb, sheet = 4, page4)
addWorksheet(wb, "dbCAN2Hit")
writeData(wb, sheet = 5, page5)
addWorksheet(wb, "MEROPSHit")
writeData(wb, sheet = 6, page6)

#wb <- loadWorkbook("Metabolic-Output.xlsx")
# Create Header Style
headerStyle <- createStyle(fontSize = 11, fontColour = "#000000", halign = "center",
                           fgFill = "#4F81BD", border="TopBottom", borderColour = "#4F81BD",textDecoration= "bold")
# Apply the bold style to all the sheets:
addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:ncol(page1), gridExpand = TRUE)
addStyle(wb, sheet = 2, headerStyle, rows = 1, cols = 1:ncol(page2), gridExpand = TRUE)
addStyle(wb, sheet = 3, headerStyle, rows = 1, cols = 1:ncol(page3), gridExpand = TRUE)
addStyle(wb, sheet = 4, headerStyle, rows = 1, cols = 1:ncol(page4), gridExpand = TRUE)
addStyle(wb, sheet = 5, headerStyle, rows = 1, cols = 1:ncol(page5), gridExpand = TRUE)
addStyle(wb, sheet = 6, headerStyle, rows = 1, cols = 1:ncol(page6), gridExpand = TRUE)
# Apply this to all the sheets

# Create a style for conditional formatting of the cells:
AbsStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
PresStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")

#We want to apply this conditional formatting to only certain sheets:
conditionalFormatting(
  wb,
  sheet = "HMMHitNum",
  cols=1:ncol(page1),
  rows=1:nrow(page1),
  style = AbsStyle,
  type = "contains",
  rule = "Absent"
)

conditionalFormatting(
  wb,
  sheet = "HMMHitNum",
  cols=1:ncol(page1),
  rows=1:nrow(page1),
  style = PresStyle,
  type = "contains",
  rule = "Present"
)

# Repeat for next sheets:
conditionalFormatting(
  wb, sheet = "FunctionHit", cols=1:ncol(page2),rows=1:nrow(page2),style = AbsStyle,type = "contains",rule = "Absent")
conditionalFormatting(
  wb, sheet = "FunctionHit", cols=1:ncol(page2),rows=1:nrow(page2),style = PresStyle,type = "contains",rule = "Present")

conditionalFormatting(
  wb, sheet = "KEGGModuleHit", cols=1:ncol(page3),rows=1:nrow(page3),style = AbsStyle,type = "contains",rule = "Absent")
conditionalFormatting(
  wb, sheet = "KEGGModuleHit", cols=1:ncol(page3),rows=1:nrow(page3),style = PresStyle,type = "contains",rule = "Present")

conditionalFormatting(
  wb, sheet = "KEGGModuleStepHit", cols=1:ncol(page4),rows=1:nrow(page4),style = AbsStyle,type = "contains",rule = "Absent")
conditionalFormatting(
  wb, sheet = "KEGGModuleStepHit", cols=1:ncol(page4),rows=1:nrow(page4),style = PresStyle,type = "contains",rule = "Present")

conditionalFormatting(
  wb, sheet = "dbCAN2Hit", cols=1:ncol(page5),rows=1:nrow(page5),style = AbsStyle,type = "contains",rule = "Absent")
conditionalFormatting(
  wb, sheet = "dbCAN2Hit", cols=1:ncol(page5),rows=1:nrow(page5),style = PresStyle,type = "contains",rule = "Present")

conditionalFormatting(
  wb, sheet = "MEROPSHit", cols=1:ncol(page6),rows=1:nrow(page6),style = AbsStyle,type = "contains",rule = "Absent")
conditionalFormatting(
  wb, sheet = "MEROPSHit", cols=1:ncol(page6),rows=1:nrow(page6),style = PresStyle,type = "contains",rule = "Present")

# Make the font to be Arial, bold and size 11.
modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Arial")

# Save it and see:
saveWorkbook(wb, "METABOLIC_result.xlsx", overwrite = TRUE)

<<<<<<< HEAD
print("Finished making the Excel Spreadsheet")
=======
print("Finished making the Excel Spreadsheet")
>>>>>>> 504d316611cab29d1ab54d6a7a952a7667806a8a

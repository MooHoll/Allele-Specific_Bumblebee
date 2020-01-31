

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0030866","cortical actin cytoskeleton organization",0.032,3.2017,0.729,0.000,"cortical actin cytoskeleton organization"),
c("GO:0032506","cytokinetic process",0.192,1.5934,0.708,0.131,"cortical actin cytoskeleton organization"),
c("GO:0030029","actin filament-based process",0.398,1.3388,0.775,0.159,"cortical actin cytoskeleton organization"),
c("GO:2001135","regulation of endocytic recycling",0.002,2.1461,0.833,0.017,"regulation of endocytic recycling"),
c("GO:0055085","transmembrane transport",8.916,1.5496,0.883,0.218,"regulation of endocytic recycling"),
c("GO:0072583","clathrin-dependent endocytosis",0.016,1.6179,0.856,0.461,"regulation of endocytic recycling"),
c("GO:0045464","R8 cell fate specification",0.000,2.8439,0.346,0.091,"R8 cell fate specification"),
c("GO:0007541","sex determination, primary response to X:A ratio",0.000,2.3674,0.487,0.300,"R8 cell fate specification"),
c("GO:0050908","detection of light stimulus involved in visual perception",0.005,1.7678,0.520,0.325,"R8 cell fate specification"),
c("GO:0035019","somatic stem cell population maintenance",0.010,1.5934,0.460,0.380,"R8 cell fate specification"),
c("GO:0061382","Malpighian tubule tip cell differentiation",0.001,2.3674,0.428,0.393,"R8 cell fate specification"),
c("GO:0007349","cellularization",0.006,1.3216,0.457,0.415,"R8 cell fate specification"),
c("GO:0061101","neuroendocrine cell differentiation",0.002,2.1461,0.489,0.510,"R8 cell fate specification"),
c("GO:0061337","cardiac conduction",0.012,1.5701,0.458,0.526,"R8 cell fate specification"),
c("GO:0008407","chaeta morphogenesis",0.002,1.5481,0.417,0.555,"R8 cell fate specification"),
c("GO:0007538","primary sex determination",0.001,1.8919,0.437,0.677,"R8 cell fate specification"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_ASE_repro.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "Significantly Enriched GO Terms, p-val<0.05",
  inflate.labels = T,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 1,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCC00",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "bottom",
  title.legend = "",
  fontsize.labels = c(0,5),
  #force.print.labels=TRUE,
  
  fontsize.legend = 8
  ,palette = "Set3"
)

dev.off()

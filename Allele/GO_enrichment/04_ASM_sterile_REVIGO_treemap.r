

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
revigo.data <- rbind(c("GO:0043484","regulation of RNA splicing",0.040,2.6834,0.751,0.000,"regulation of RNA splicing"),
c("GO:0019284","L-methionine biosynthetic process from S-adenosylmethionine",0.037,1.5815,0.839,0.124,"regulation of RNA splicing"),
c("GO:0034474","U2 snRNA 3'-end processing",0.000,1.5815,0.885,0.329,"regulation of RNA splicing"),
c("GO:0000349","generation of catalytic spliceosome for first transesterification step",0.002,1.5815,0.801,0.587,"regulation of RNA splicing"),
c("GO:0000096","sulfur amino acid metabolic process",0.539,1.4230,0.854,0.599,"regulation of RNA splicing"),
c("GO:0000381","regulation of alternative mRNA splicing, via spliceosome",0.017,2.0466,0.745,0.662,"regulation of RNA splicing"),
c("GO:0045087","innate immune response",0.148,2.7629,0.750,0.000,"innate immune response"),
c("GO:0035332","positive regulation of hippo signaling",0.001,2.6989,0.695,0.201,"innate immune response"),
c("GO:2000286","receptor internalization involved in canonical Wnt signaling pathway",0.000,1.5815,0.716,0.211,"innate immune response"),
c("GO:0007166","cell surface receptor signaling pathway",0.920,1.4213,0.672,0.330,"innate immune response"),
c("GO:0034198","cellular response to amino acid starvation",0.016,1.6642,0.795,0.385,"innate immune response"),
c("GO:0046824","positive regulation of nucleocytoplasmic transport",0.029,1.6642,0.771,0.397,"innate immune response"),
c("GO:0050794","regulation of cellular process",18.840,1.4283,0.779,0.416,"innate immune response"),
c("GO:0051647","nucleus localization",0.012,1.4948,0.918,0.444,"innate immune response"),
c("GO:0023057","negative regulation of signaling",0.290,2.0798,0.671,0.523,"innate immune response"),
c("GO:0051716","cellular response to stimulus",9.561,1.4830,0.781,0.525,"innate immune response"),
c("GO:0007516","hemocyte development",0.000,1.5815,0.609,0.531,"innate immune response"),
c("GO:0034125","negative regulation of MyD88-dependent toll-like receptor signaling pathway",0.000,1.5815,0.656,0.566,"innate immune response"),
c("GO:0045879","negative regulation of smoothened signaling pathway",0.007,1.4230,0.660,0.599,"innate immune response"),
c("GO:0002314","germinal center B cell differentiation",0.001,1.5815,0.524,0.619,"innate immune response"),
c("GO:0035745","T-helper 2 cell cytokine production",0.002,1.5815,0.579,0.655,"innate immune response"),
c("GO:0010648","negative regulation of cell communication",0.290,2.0250,0.688,0.658,"innate immune response"),
c("GO:0048585","negative regulation of response to stimulus",0.344,1.7261,0.670,0.667,"innate immune response"),
c("GO:0065007","biological regulation",20.498,1.3755,0.977,0.000,"biological regulation"),
c("GO:0097150","neuronal stem cell population maintenance",0.005,2.6989,0.642,0.000,"neuronal stem cell population maintenance"),
c("GO:0030717","karyosome formation",0.002,2.4053,0.582,0.349,"neuronal stem cell population maintenance"),
c("GO:0006723","cuticle hydrocarbon biosynthetic process",0.000,1.5815,0.668,0.349,"neuronal stem cell population maintenance"),
c("GO:0046960","sensitization",0.001,1.5815,0.713,0.357,"neuronal stem cell population maintenance"),
c("GO:0097374","sensory neuron axon guidance",0.002,1.5815,0.528,0.465,"neuronal stem cell population maintenance"),
c("GO:0009653","anatomical structure morphogenesis",1.542,2.0587,0.606,0.494,"neuronal stem cell population maintenance"),
c("GO:0007391","dorsal closure",0.005,1.7854,0.590,0.507,"neuronal stem cell population maintenance"),
c("GO:0007446","imaginal disc growth",0.002,1.7660,0.591,0.561,"neuronal stem cell population maintenance"),
c("GO:2000255","negative regulation of male germ cell proliferation",0.001,1.5815,0.569,0.568,"neuronal stem cell population maintenance"),
c("GO:0098727","maintenance of cell number",0.036,1.3026,0.681,0.573,"neuronal stem cell population maintenance"),
c("GO:0072205","metanephric collecting duct development",0.001,1.5815,0.625,0.578,"neuronal stem cell population maintenance"),
c("GO:0001764","neuron migration",0.030,1.4230,0.556,0.584,"neuronal stem cell population maintenance"),
c("GO:0060465","pharynx development",0.002,1.5815,0.629,0.585,"neuronal stem cell population maintenance"),
c("GO:0042473","outer ear morphogenesis",0.002,1.5815,0.601,0.607,"neuronal stem cell population maintenance"),
c("GO:0008359","regulation of bicoid mRNA localization",0.000,1.5815,0.540,0.610,"neuronal stem cell population maintenance"),
c("GO:0003402","planar cell polarity pathway involved in axis elongation",0.001,1.5815,0.466,0.611,"neuronal stem cell population maintenance"),
c("GO:0071678","olfactory bulb axon guidance",0.000,1.5815,0.549,0.628,"neuronal stem cell population maintenance"),
c("GO:0048749","compound eye development",0.020,1.3056,0.587,0.659,"neuronal stem cell population maintenance"),
c("GO:0060441","epithelial tube branching involved in lung morphogenesis",0.006,1.5815,0.578,0.659,"neuronal stem cell population maintenance"),
c("GO:0060711","labyrinthine layer development",0.008,1.5815,0.571,0.667,"neuronal stem cell population maintenance"),
c("GO:0007300","ovarian nurse cell to oocyte transport",0.001,1.5746,0.590,0.677,"neuronal stem cell population maintenance"),
c("GO:0001578","microtubule bundle formation",0.027,1.3577,0.835,0.045,"microtubule bundle formation"),
c("GO:0018095","protein polyglutamylation",0.007,1.5815,0.870,0.081,"protein polyglutamylation"),
c("GO:0006468","protein phosphorylation",4.137,1.5351,0.852,0.370,"protein polyglutamylation"),
c("GO:0051567","histone H3-K9 methylation",0.010,1.5746,0.815,0.440,"protein polyglutamylation"),
c("GO:0038083","peptidyl-tyrosine autophosphorylation",0.011,1.5746,0.865,0.451,"protein polyglutamylation"),
c("GO:0000165","MAPK cascade",0.245,1.3275,0.627,0.514,"protein polyglutamylation"),
c("GO:0018212","peptidyl-tyrosine modification",0.215,1.4112,0.849,0.543,"protein polyglutamylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ASM_sterile_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

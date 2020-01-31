

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
revigo.data <- rbind(c("GO:0033120","positive regulation of RNA splicing",0.007,2.9883,0.666,0.000,"positive regulation of RNA splicing"),
c("GO:0043610","regulation of carbohydrate utilization",0.000,1.5589,0.871,0.109,"positive regulation of RNA splicing"),
c("GO:0045838","positive regulation of membrane potential",0.001,1.3816,0.840,0.122,"positive regulation of RNA splicing"),
c("GO:0042749","regulation of circadian sleep/wake cycle",0.005,2.1467,0.668,0.136,"positive regulation of RNA splicing"),
c("GO:0042789","mRNA transcription from RNA polymerase II promoter",0.004,2.6539,0.892,0.157,"positive regulation of RNA splicing"),
c("GO:0016071","mRNA metabolic process",0.798,1.8630,0.880,0.215,"positive regulation of RNA splicing"),
c("GO:0034654","nucleobase-containing compound biosynthetic process",14.533,1.6971,0.859,0.288,"positive regulation of RNA splicing"),
c("GO:0051253","negative regulation of RNA metabolic process",0.626,1.8416,0.714,0.312,"positive regulation of RNA splicing"),
c("GO:0021553","olfactory nerve development",0.001,1.5589,0.803,0.331,"positive regulation of RNA splicing"),
c("GO:1904056","positive regulation of cholangiocyte proliferation",0.000,1.5589,0.775,0.333,"positive regulation of RNA splicing"),
c("GO:0002128","tRNA nucleoside ribose methylation",0.011,1.5589,0.883,0.356,"positive regulation of RNA splicing"),
c("GO:0035158","regulation of tube diameter, open tracheal system",0.001,1.5589,0.715,0.358,"positive regulation of RNA splicing"),
c("GO:1904594","regulation of termination of RNA polymerase II transcription",0.000,1.5589,0.742,0.369,"positive regulation of RNA splicing"),
c("GO:0034641","cellular nitrogen compound metabolic process",34.137,1.3402,0.903,0.387,"positive regulation of RNA splicing"),
c("GO:1904193","negative regulation of cholangiocyte apoptotic process",0.000,1.5589,0.840,0.396,"positive regulation of RNA splicing"),
c("GO:0007389","pattern specification process",0.147,1.6653,0.759,0.442,"positive regulation of RNA splicing"),
c("GO:0097309","cap1 mRNA methylation",0.000,1.5589,0.895,0.448,"positive regulation of RNA splicing"),
c("GO:0090702","non-reproductive fruiting body development",0.008,1.5589,0.879,0.448,"positive regulation of RNA splicing"),
c("GO:0048752","semicircular canal morphogenesis",0.002,1.5589,0.788,0.475,"positive regulation of RNA splicing"),
c("GO:0035167","larval lymph gland hemopoiesis",0.002,1.4530,0.793,0.477,"positive regulation of RNA splicing"),
c("GO:0035886","vascular smooth muscle cell differentiation",0.004,1.5589,0.760,0.491,"positive regulation of RNA splicing"),
c("GO:0060999","positive regulation of dendritic spine development",0.004,1.8401,0.551,0.502,"positive regulation of RNA splicing"),
c("GO:0010611","regulation of cardiac muscle hypertrophy",0.006,1.7230,0.643,0.509,"positive regulation of RNA splicing"),
c("GO:0021762","substantia nigra development",0.008,1.5589,0.774,0.524,"positive regulation of RNA splicing"),
c("GO:0031529","ruffle organization",0.010,1.8401,0.859,0.527,"positive regulation of RNA splicing"),
c("GO:0051703","intraspecies interaction between organisms",0.017,1.3816,0.969,0.536,"positive regulation of RNA splicing"),
c("GO:0019219","regulation of nucleobase-containing compound metabolic process",10.258,2.0642,0.690,0.537,"positive regulation of RNA splicing"),
c("GO:1904197","positive regulation of granulosa cell proliferation",0.000,1.5589,0.777,0.540,"positive regulation of RNA splicing"),
c("GO:0051549","positive regulation of keratinocyte migration",0.002,1.5589,0.546,0.549,"positive regulation of RNA splicing"),
c("GO:0090304","nucleic acid metabolic process",21.449,1.4744,0.868,0.562,"positive regulation of RNA splicing"),
c("GO:0060566","positive regulation of DNA-templated transcription, termination",0.000,1.5589,0.673,0.563,"positive regulation of RNA splicing"),
c("GO:0009953","dorsal/ventral pattern formation",0.034,1.5996,0.767,0.565,"positive regulation of RNA splicing"),
c("GO:0003215","cardiac right ventricle morphogenesis",0.004,1.5589,0.777,0.596,"positive regulation of RNA splicing"),
c("GO:0044057","regulation of system process",0.088,1.3506,0.671,0.599,"positive regulation of RNA splicing"),
c("GO:1902373","negative regulation of mRNA catabolic process",0.002,1.8401,0.758,0.624,"positive regulation of RNA splicing"),
c("GO:0003062","regulation of heart rate by chemical signal",0.001,1.5589,0.702,0.648,"positive regulation of RNA splicing"),
c("GO:0007628","adult walking behavior",0.008,1.9784,0.768,0.655,"positive regulation of RNA splicing"),
c("GO:0060421","positive regulation of heart growth",0.005,1.5589,0.591,0.667,"positive regulation of RNA splicing"),
c("GO:0048869","cellular developmental process",1.896,1.4755,0.792,0.670,"positive regulation of RNA splicing"),
c("GO:0016480","negative regulation of transcription from RNA polymerase III promoter",0.009,1.5589,0.774,0.683,"positive regulation of RNA splicing"),
c("GO:0002087","regulation of respiratory gaseous exchange by neurological system process",0.002,1.5589,0.698,0.687,"positive regulation of RNA splicing"),
c("GO:0030719","P granule organization",0.001,1.5589,0.715,0.688,"positive regulation of RNA splicing"),
c("GO:2000806","positive regulation of termination of RNA polymerase II transcription, poly(A)-coupled",0.000,1.5589,0.671,0.690,"positive regulation of RNA splicing"),
c("GO:0051931","regulation of sensory perception",0.003,1.5589,0.704,0.691,"positive regulation of RNA splicing"),
c("GO:0021510","spinal cord development",0.021,1.3167,0.764,0.697,"positive regulation of RNA splicing"),
c("GO:1900029","positive regulation of ruffle assembly",0.001,1.5589,0.684,0.699,"positive regulation of RNA splicing"),
c("GO:0040007","growth",0.317,1.4681,0.988,0.000,"growth"),
c("GO:0098655","cation transmembrane transport",2.290,1.8723,0.849,0.000,"cation transmembrane transport"),
c("GO:0099003","vesicle-mediated transport in synapse",0.036,1.4992,0.872,0.244,"cation transmembrane transport"),
c("GO:0051028","mRNA transport",0.075,1.7456,0.848,0.261,"cation transmembrane transport"),
c("GO:0019090","mitochondrial rRNA export from mitochondrion",0.000,1.5589,0.810,0.646,"cation transmembrane transport"),
c("GO:0071229","cellular response to acid chemical",0.066,2.8174,0.863,0.019,"cellular response to acid chemical"),
c("GO:0023041","neuronal signal transduction",0.001,1.5589,0.790,0.223,"cellular response to acid chemical"),
c("GO:0031668","cellular response to extracellular stimulus",0.436,2.6112,0.822,0.333,"cellular response to acid chemical"),
c("GO:0046928","regulation of neurotransmitter secretion",0.012,2.5451,0.580,0.337,"cellular response to acid chemical"),
c("GO:0009914","hormone transport",0.074,1.3543,0.718,0.402,"cellular response to acid chemical"),
c("GO:0042542","response to hydrogen peroxide",0.029,2.1392,0.891,0.437,"cellular response to acid chemical"),
c("GO:1902474","positive regulation of protein localization to synapse",0.001,1.5589,0.667,0.474,"cellular response to acid chemical"),
c("GO:0010890","positive regulation of sequestering of triglyceride",0.001,1.5589,0.647,0.485,"cellular response to acid chemical"),
c("GO:0006606","protein import into nucleus",0.102,2.4091,0.752,0.505,"cellular response to acid chemical"),
c("GO:0042493","response to drug",0.266,1.3180,0.894,0.510,"cellular response to acid chemical"),
c("GO:0098942","retrograde trans-synaptic signaling by trans-synaptic protein complex",0.000,1.5589,0.871,0.516,"cellular response to acid chemical"),
c("GO:0031887","lipid particle transport along microtubule",0.000,1.5589,0.783,0.534,"cellular response to acid chemical"),
c("GO:0007637","proboscis extension reflex",0.000,1.5589,0.820,0.585,"cellular response to acid chemical"),
c("GO:1900244","positive regulation of synaptic vesicle endocytosis",0.001,2.3607,0.548,0.593,"cellular response to acid chemical"),
c("GO:1901950","dense core granule transport",0.000,1.5589,0.830,0.597,"cellular response to acid chemical"),
c("GO:1990253","cellular response to leucine starvation",0.001,2.6539,0.857,0.649,"cellular response to acid chemical"),
c("GO:0051656","establishment of organelle localization",0.180,2.0604,0.796,0.666,"cellular response to acid chemical"),
c("GO:0034765","regulation of ion transmembrane transport",0.197,1.9210,0.661,0.683,"cellular response to acid chemical"),
c("GO:0023061","signal release",0.092,1.4769,0.736,0.684,"cellular response to acid chemical"),
c("GO:0032252","secretory granule localization",0.001,1.5589,0.833,0.685,"cellular response to acid chemical"),
c("GO:0043933","macromolecular complex subunit organization",2.371,2.6964,0.868,0.029,"macromolecular complex subunit organization"),
c("GO:0072553","terminal button organization",0.003,1.6216,0.861,0.342,"macromolecular complex subunit organization"),
c("GO:0000022","mitotic spindle elongation",0.006,2.1467,0.864,0.364,"macromolecular complex subunit organization"),
c("GO:0030034","microvillar actin bundle assembly",0.000,1.5589,0.872,0.404,"macromolecular complex subunit organization"),
c("GO:0016569","covalent chromatin modification",0.424,2.2089,0.840,0.529,"macromolecular complex subunit organization"),
c("GO:0051262","protein tetramerization",0.044,1.4530,0.872,0.554,"macromolecular complex subunit organization"),
c("GO:0051276","chromosome organization",1.477,1.3028,0.867,0.612,"macromolecular complex subunit organization"),
c("GO:0045840","positive regulation of mitotic nuclear division",0.017,1.5324,0.652,0.625,"macromolecular complex subunit organization"),
c("GO:1905606","regulation of presynapse assembly",0.001,1.5589,0.618,0.647,"macromolecular complex subunit organization"),
c("GO:0043984","histone H4-K16 acetylation",0.006,1.5324,0.868,0.655,"macromolecular complex subunit organization"),
c("GO:1901360","organic cyclic compound metabolic process",30.324,1.7448,0.962,0.042,"organic cyclic compound metabolism"),
c("GO:0006725","cellular aromatic compound metabolic process",29.628,1.9141,0.933,0.066,"cellular aromatic compound metabolism"),
c("GO:0046483","heterocycle metabolic process",29.664,1.7276,0.933,0.245,"cellular aromatic compound metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ASM_both_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

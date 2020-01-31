

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
revigo.data <- rbind(c("GO:0043279","response to alkaloid",0.017,1.5960,0.945,0.000,"response to alkaloid"),
c("GO:0098942","retrograde trans-synaptic signaling by trans-synaptic protein complex",0.000,3.1230,0.816,0.000,"retrograde trans-synaptic signaling by trans-synaptic protein complex"),
c("GO:0023041","neuronal signal transduction",0.001,3.1230,0.804,0.186,"retrograde trans-synaptic signaling by trans-synaptic protein complex"),
c("GO:0032230","positive regulation of synaptic transmission, GABAergic",0.002,2.3455,0.659,0.479,"retrograde trans-synaptic signaling by trans-synaptic protein complex"),
c("GO:0099565","chemical synaptic transmission, postsynaptic",0.066,1.7114,0.672,0.688,"retrograde trans-synaptic signaling by trans-synaptic protein complex"),
c("GO:1902474","positive regulation of protein localization to synapse",0.001,3.1230,0.687,0.000,"positive regulation of protein localization to synapse"),
c("GO:1905606","regulation of presynapse assembly",0.001,3.1230,0.452,0.113,"positive regulation of protein localization to synapse"),
c("GO:0048789","cytoskeletal matrix organization at active zone",0.000,2.4246,0.808,0.200,"positive regulation of protein localization to synapse"),
c("GO:0043063","intercellular bridge organization",0.001,2.1698,0.774,0.208,"positive regulation of protein localization to synapse"),
c("GO:0045838","positive regulation of membrane potential",0.001,1.9209,0.797,0.284,"positive regulation of protein localization to synapse"),
c("GO:0006610","ribosomal protein import into nucleus",0.004,2.3455,0.778,0.326,"positive regulation of protein localization to synapse"),
c("GO:0046960","sensitization",0.001,2.8221,0.694,0.327,"positive regulation of protein localization to synapse"),
c("GO:0051290","protein heterotetramerization",0.004,1.8947,0.783,0.328,"positive regulation of protein localization to synapse"),
c("GO:0044091","membrane biogenesis",0.096,1.8467,0.782,0.352,"positive regulation of protein localization to synapse"),
c("GO:0007485","imaginal disc-derived male genitalia development",0.001,1.9786,0.660,0.371,"positive regulation of protein localization to synapse"),
c("GO:0008258","head involution",0.002,1.9488,0.687,0.373,"positive regulation of protein localization to synapse"),
c("GO:0035277","spiracle morphogenesis, open tracheal system",0.001,1.9488,0.698,0.381,"positive regulation of protein localization to synapse"),
c("GO:0016200","synaptic target attraction",0.000,2.5213,0.648,0.399,"positive regulation of protein localization to synapse"),
c("GO:0007349","cellularization",0.006,1.4952,0.680,0.403,"positive regulation of protein localization to synapse"),
c("GO:2000809","positive regulation of synaptic vesicle clustering",0.001,2.8221,0.641,0.424,"positive regulation of protein localization to synapse"),
c("GO:0043576","regulation of respiratory gaseous exchange",0.004,2.5213,0.612,0.457,"positive regulation of protein localization to synapse"),
c("GO:2000310","regulation of N-methyl-D-aspartate selective glutamate receptor activity",0.002,2.8221,0.649,0.470,"positive regulation of protein localization to synapse"),
c("GO:1900029","positive regulation of ruffle assembly",0.001,3.1230,0.587,0.471,"positive regulation of protein localization to synapse"),
c("GO:0048468","cell development",0.573,1.4097,0.639,0.481,"positive regulation of protein localization to synapse"),
c("GO:0044087","regulation of cellular component biogenesis",0.404,1.6068,0.665,0.481,"positive regulation of protein localization to synapse"),
c("GO:0031529","ruffle organization",0.010,2.0453,0.710,0.494,"positive regulation of protein localization to synapse"),
c("GO:0010469","regulation of receptor activity",0.025,1.8467,0.729,0.497,"positive regulation of protein localization to synapse"),
c("GO:0007435","salivary gland morphogenesis",0.016,1.3978,0.658,0.515,"positive regulation of protein localization to synapse"),
c("GO:0031644","regulation of neurological system process",0.008,1.6221,0.625,0.516,"positive regulation of protein localization to synapse"),
c("GO:0010841","positive regulation of circadian sleep/wake cycle, wakefulness",0.000,2.8221,0.569,0.521,"positive regulation of protein localization to synapse"),
c("GO:0008340","determination of adult lifespan",0.020,1.3600,0.669,0.524,"positive regulation of protein localization to synapse"),
c("GO:0007215","glutamate receptor signaling pathway",0.096,1.5714,0.739,0.526,"positive regulation of protein localization to synapse"),
c("GO:0070192","chromosome organization involved in meiotic cell cycle",0.031,1.6221,0.690,0.539,"positive regulation of protein localization to synapse"),
c("GO:0097061","dendritic spine organization",0.009,1.8700,0.711,0.550,"positive regulation of protein localization to synapse"),
c("GO:0032409","regulation of transporter activity",0.039,1.4952,0.695,0.562,"positive regulation of protein localization to synapse"),
c("GO:0006886","intracellular protein transport",1.199,1.3920,0.803,0.563,"positive regulation of protein localization to synapse"),
c("GO:0051590","positive regulation of neurotransmitter transport",0.003,2.0107,0.631,0.574,"positive regulation of protein localization to synapse"),
c("GO:0061002","negative regulation of dendritic spine morphogenesis",0.001,2.5213,0.493,0.583,"positive regulation of protein localization to synapse"),
c("GO:0097114","NMDA glutamate receptor clustering",0.001,2.8221,0.597,0.603,"positive regulation of protein localization to synapse"),
c("GO:0030833","regulation of actin filament polymerization",0.123,1.3900,0.545,0.608,"positive regulation of protein localization to synapse"),
c("GO:0007301","female germline ring canal formation",0.000,2.3455,0.549,0.623,"positive regulation of protein localization to synapse"),
c("GO:0099084","postsynaptic specialization organization",0.002,2.3455,0.679,0.639,"positive regulation of protein localization to synapse"),
c("GO:1990709","presynaptic active zone organization",0.001,2.2787,0.697,0.645,"positive regulation of protein localization to synapse"),
c("GO:0030717","karyosome formation",0.002,2.2208,0.586,0.646,"positive regulation of protein localization to synapse"),
c("GO:1901631","positive regulation of presynaptic membrane organization",0.000,2.8221,0.573,0.647,"positive regulation of protein localization to synapse"),
c("GO:0007622","rhythmic behavior",0.012,1.4058,0.912,0.671,"positive regulation of protein localization to synapse"),
c("GO:0007619","courtship behavior",0.005,1.3389,0.690,0.673,"positive regulation of protein localization to synapse"),
c("GO:1904062","regulation of cation transmembrane transport",0.046,1.3978,0.692,0.676,"positive regulation of protein localization to synapse"),
c("GO:0072553","terminal button organization",0.003,1.8245,0.683,0.681,"positive regulation of protein localization to synapse"),
c("GO:1904861","excitatory synapse assembly",0.001,2.5213,0.529,0.684,"positive regulation of protein localization to synapse"),
c("GO:1900244","positive regulation of synaptic vesicle endocytosis",0.001,2.5213,0.515,0.685,"positive regulation of protein localization to synapse"),
c("GO:0016080","synaptic vesicle targeting",0.001,2.2208,0.588,0.692,"positive regulation of protein localization to synapse"),
c("GO:0007158","neuron cell-cell adhesion",0.002,2.0453,0.875,0.032,"neuron cell-cell adhesion"),
c("GO:0007157","heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules",0.008,1.9209,0.924,0.597,"neuron cell-cell adhesion"),
c("GO:0016339","calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules",0.004,1.7834,0.926,0.699,"neuron cell-cell adhesion"),
c("GO:0042023","DNA endoreduplication",0.007,2.1242,0.861,0.075,"DNA endoreduplication"),
c("GO:0038083","peptidyl-tyrosine autophosphorylation",0.011,2.0107,0.948,0.075,"peptidyl-tyrosine autophosphorylation"),
c("GO:0018212","peptidyl-tyrosine modification",0.215,1.3458,0.945,0.543,"peptidyl-tyrosine autophosphorylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_ASE_and_ASM.pdf", width=16, height=9 ) # width and height are in inches

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

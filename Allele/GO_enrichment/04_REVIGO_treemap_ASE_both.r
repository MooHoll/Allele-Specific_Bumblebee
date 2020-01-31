

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
revigo.data <- rbind(c("GO:0010501","RNA secondary structure unwinding",0.025,1.5420,0.949,0.000,"RNA secondary structure unwinding"),
c("GO:0023052","signaling",6.765,1.3769,0.990,0.000,"signaling"),
c("GO:0035114","imaginal disc-derived appendage morphogenesis",0.021,2.9029,0.667,0.000,"imaginal disc-derived appendage morphogenesis"),
c("GO:0014806","smooth muscle hyperplasia",0.000,2.0058,0.714,0.338,"imaginal disc-derived appendage morphogenesis"),
c("GO:0032764","negative regulation of mast cell cytokine production",0.000,2.0058,0.562,0.370,"imaginal disc-derived appendage morphogenesis"),
c("GO:0060288","formation of a compartment boundary",0.000,1.5329,0.732,0.380,"imaginal disc-derived appendage morphogenesis"),
c("GO:0018990","ecdysis, chitin-based cuticle",0.002,1.4100,0.723,0.403,"imaginal disc-derived appendage morphogenesis"),
c("GO:0035810","positive regulation of urine volume",0.002,1.5329,0.655,0.413,"imaginal disc-derived appendage morphogenesis"),
c("GO:0031622","positive regulation of fever generation",0.001,2.0058,0.513,0.418,"imaginal disc-derived appendage morphogenesis"),
c("GO:0008258","head involution",0.002,2.0333,0.697,0.430,"imaginal disc-derived appendage morphogenesis"),
c("GO:1905606","regulation of presynapse assembly",0.001,2.0058,0.538,0.433,"imaginal disc-derived appendage morphogenesis"),
c("GO:0021942","radial glia guided migration of Purkinje cell",0.001,1.5329,0.653,0.435,"imaginal disc-derived appendage morphogenesis"),
c("GO:0048872","homeostasis of number of cells",0.086,1.9813,0.807,0.438,"imaginal disc-derived appendage morphogenesis"),
c("GO:1902903","regulation of supramolecular fiber organization",0.166,1.9328,0.703,0.445,"imaginal disc-derived appendage morphogenesis"),
c("GO:0061484","hematopoietic stem cell homeostasis",0.002,1.5329,0.839,0.451,"imaginal disc-derived appendage morphogenesis"),
c("GO:0046849","bone remodeling",0.015,2.3875,0.694,0.456,"imaginal disc-derived appendage morphogenesis"),
c("GO:0060078","regulation of postsynaptic membrane potential",0.077,1.5727,0.807,0.459,"imaginal disc-derived appendage morphogenesis"),
c("GO:0007292","female gamete generation",0.060,2.8435,0.664,0.460,"imaginal disc-derived appendage morphogenesis"),
c("GO:1901223","negative regulation of NIK/NF-kappaB signaling",0.002,2.0058,0.679,0.478,"imaginal disc-derived appendage morphogenesis"),
c("GO:0048563","post-embryonic animal organ morphogenesis",0.026,2.4799,0.646,0.492,"imaginal disc-derived appendage morphogenesis"),
c("GO:0048736","appendage development",0.062,2.2755,0.662,0.524,"imaginal disc-derived appendage morphogenesis"),
c("GO:0007552","metamorphosis",0.031,2.1284,0.673,0.530,"imaginal disc-derived appendage morphogenesis"),
c("GO:0048813","dendrite morphogenesis",0.034,2.2432,0.573,0.533,"imaginal disc-derived appendage morphogenesis"),
c("GO:0016200","synaptic target attraction",0.000,1.4100,0.645,0.544,"imaginal disc-derived appendage morphogenesis"),
c("GO:0010642","negative regulation of platelet-derived growth factor receptor signaling pathway",0.002,1.3152,0.667,0.547,"imaginal disc-derived appendage morphogenesis"),
c("GO:0003012","muscle system process",0.079,1.3984,0.709,0.549,"imaginal disc-derived appendage morphogenesis"),
c("GO:0007500","mesodermal cell fate determination",0.000,1.5329,0.681,0.550,"imaginal disc-derived appendage morphogenesis"),
c("GO:0048815","hermaphrodite genitalia morphogenesis",0.000,1.7069,0.701,0.567,"imaginal disc-derived appendage morphogenesis"),
c("GO:0051491","positive regulation of filopodium assembly",0.005,2.2267,0.686,0.569,"imaginal disc-derived appendage morphogenesis"),
c("GO:0014827","intestine smooth muscle contraction",0.001,1.7069,0.760,0.570,"imaginal disc-derived appendage morphogenesis"),
c("GO:0048749","compound eye development",0.020,1.9221,0.668,0.574,"imaginal disc-derived appendage morphogenesis"),
c("GO:0045080","positive regulation of chemokine biosynthetic process",0.002,2.0058,0.563,0.578,"imaginal disc-derived appendage morphogenesis"),
c("GO:0001826","inner cell mass cell differentiation",0.001,1.7069,0.676,0.589,"imaginal disc-derived appendage morphogenesis"),
c("GO:0051241","negative regulation of multicellular organismal process",0.211,1.5921,0.506,0.593,"imaginal disc-derived appendage morphogenesis"),
c("GO:0007444","imaginal disc development",0.034,1.7197,0.646,0.594,"imaginal disc-derived appendage morphogenesis"),
c("GO:0002246","wound healing involved in inflammatory response",0.001,1.7069,0.881,0.598,"imaginal disc-derived appendage morphogenesis"),
c("GO:0033004","negative regulation of mast cell activation",0.002,1.7069,0.747,0.601,"imaginal disc-derived appendage morphogenesis"),
c("GO:0097114","NMDA glutamate receptor clustering",0.001,2.0058,0.749,0.603,"imaginal disc-derived appendage morphogenesis"),
c("GO:0000710","meiotic mismatch repair",0.006,1.5329,0.699,0.609,"imaginal disc-derived appendage morphogenesis"),
c("GO:0021799","cerebral cortex radially oriented cell migration",0.006,1.3152,0.621,0.616,"imaginal disc-derived appendage morphogenesis"),
c("GO:0043576","regulation of respiratory gaseous exchange",0.004,1.4100,0.626,0.619,"imaginal disc-derived appendage morphogenesis"),
c("GO:0061458","reproductive system development",0.155,1.4862,0.641,0.622,"imaginal disc-derived appendage morphogenesis"),
c("GO:0045035","sensory organ precursor cell division",0.001,1.3152,0.688,0.623,"imaginal disc-derived appendage morphogenesis"),
c("GO:0001778","plasma membrane repair",0.002,1.5329,0.757,0.626,"imaginal disc-derived appendage morphogenesis"),
c("GO:0009967","positive regulation of signal transduction",0.329,1.9041,0.589,0.627,"imaginal disc-derived appendage morphogenesis"),
c("GO:0009791","post-embryonic development",0.163,1.3290,0.645,0.636,"imaginal disc-derived appendage morphogenesis"),
c("GO:0099084","postsynaptic specialization organization",0.002,1.4100,0.791,0.639,"imaginal disc-derived appendage morphogenesis"),
c("GO:1990709","presynaptic active zone organization",0.001,1.3152,0.805,0.645,"imaginal disc-derived appendage morphogenesis"),
c("GO:1901631","positive regulation of presynaptic membrane organization",0.000,1.7069,0.679,0.647,"imaginal disc-derived appendage morphogenesis"),
c("GO:0051493","regulation of cytoskeleton organization",0.232,1.3405,0.724,0.654,"imaginal disc-derived appendage morphogenesis"),
c("GO:0007313","maternal specification of dorsal/ventral axis, oocyte, soma encoded",0.001,1.5329,0.671,0.660,"imaginal disc-derived appendage morphogenesis"),
c("GO:0060279","positive regulation of ovulation",0.001,2.0058,0.592,0.672,"imaginal disc-derived appendage morphogenesis"),
c("GO:0032230","positive regulation of synaptic transmission, GABAergic",0.002,1.5329,0.684,0.676,"imaginal disc-derived appendage morphogenesis"),
c("GO:0044089","positive regulation of cellular component biogenesis",0.193,1.3119,0.700,0.682,"imaginal disc-derived appendage morphogenesis"),
c("GO:1904861","excitatory synapse assembly",0.001,1.7069,0.604,0.684,"imaginal disc-derived appendage morphogenesis"),
c("GO:0007300","ovarian nurse cell to oocyte transport",0.001,1.9264,0.638,0.693,"imaginal disc-derived appendage morphogenesis"),
c("GO:0001659","temperature homeostasis",0.006,1.3152,0.615,0.694,"imaginal disc-derived appendage morphogenesis"),
c("GO:0043950","positive regulation of cAMP-mediated signaling",0.003,1.4100,0.671,0.695,"imaginal disc-derived appendage morphogenesis"),
c("GO:0032012","regulation of ARF protein signal transduction",0.040,1.3152,0.678,0.696,"imaginal disc-derived appendage morphogenesis"),
c("GO:0042745","circadian sleep/wake cycle",0.005,2.0927,0.952,0.000,"circadian sleep/wake cycle"),
c("GO:0046960","sensitization",0.001,1.7069,0.708,0.581,"circadian sleep/wake cycle"),
c("GO:0007635","chemosensory behavior",0.017,1.9376,0.801,0.689,"circadian sleep/wake cycle"),
c("GO:0050896","response to stimulus",12.210,1.9881,0.990,0.000,"response to stimulus"),
c("GO:0006915","apoptotic process",0.406,2.2829,0.777,0.058,"apoptotic process"),
c("GO:0051790","short-chain fatty acid biosynthetic process",0.001,2.0058,0.885,0.106,"apoptotic process"),
c("GO:0006788","heme oxidation",0.015,2.0058,0.854,0.131,"apoptotic process"),
c("GO:0017157","regulation of exocytosis",0.044,2.0414,0.676,0.141,"apoptotic process"),
c("GO:1904706","negative regulation of vascular smooth muscle cell proliferation",0.002,2.0058,0.771,0.145,"apoptotic process"),
c("GO:0008219","cell death",0.458,2.3049,0.852,0.172,"apoptotic process"),
c("GO:0072488","ammonium transmembrane transport",0.086,1.7069,0.909,0.193,"apoptotic process"),
c("GO:0044205","'de novo' UMP biosynthetic process",0.221,1.5329,0.835,0.217,"apoptotic process"),
c("GO:2000117","negative regulation of cysteine-type endopeptidase activity",0.018,1.9264,0.748,0.401,"apoptotic process"),
c("GO:0051187","cofactor catabolic process",0.062,1.6733,0.951,0.410,"apoptotic process"),
c("GO:0090119","vesicle-mediated cholesterol transport",0.000,2.0058,0.837,0.432,"apoptotic process"),
c("GO:0051792","medium-chain fatty acid biosynthetic process",0.001,2.0058,0.881,0.432,"apoptotic process"),
c("GO:0031954","positive regulation of protein autophosphorylation",0.005,1.3152,0.760,0.449,"apoptotic process"),
c("GO:1903356","positive regulation of distal tip cell migration",0.000,2.0058,0.711,0.490,"apoptotic process"),
c("GO:1902474","positive regulation of protein localization to synapse",0.001,2.0058,0.750,0.506,"apoptotic process"),
c("GO:2000809","positive regulation of synaptic vesicle clustering",0.001,1.7069,0.731,0.543,"apoptotic process"),
c("GO:0033015","tetrapyrrole catabolic process",0.025,1.7069,0.942,0.556,"apoptotic process"),
c("GO:2000310","regulation of N-methyl-D-aspartate selective glutamate receptor activity",0.002,2.0058,0.641,0.565,"apoptotic process"),
c("GO:0033002","muscle cell proliferation",0.027,1.3346,0.909,0.575,"apoptotic process"),
c("GO:0006909","phagocytosis",0.051,1.3843,0.824,0.584,"apoptotic process"),
c("GO:1990927","calcium ion regulated lysosome exocytosis",0.000,1.7069,0.805,0.616,"apoptotic process"),
c("GO:0050432","catecholamine secretion",0.006,1.7069,0.790,0.626,"apoptotic process"),
c("GO:0060455","negative regulation of gastric acid secretion",0.000,2.0058,0.556,0.630,"apoptotic process"),
c("GO:1903298","negative regulation of hypoxia-induced intrinsic apoptotic signaling pathway",0.001,2.0058,0.618,0.632,"apoptotic process"),
c("GO:0070092","regulation of glucagon secretion",0.001,1.7069,0.646,0.678,"apoptotic process"),
c("GO:1900244","positive regulation of synaptic vesicle endocytosis",0.001,1.4100,0.650,0.685,"apoptotic process"),
c("GO:0010656","negative regulation of muscle cell apoptotic process",0.004,1.3152,0.718,0.691,"apoptotic process"),
c("GO:1904901","positive regulation of myosin II filament organization",0.000,2.0058,0.729,0.099,"positive regulation of myosin II filament organization"),
c("GO:0038168","epidermal growth factor receptor signaling pathway via I-kappaB kinase/NF-kappaB cascade",0.000,2.0058,0.771,0.103,"positive regulation of myosin II filament organization"),
c("GO:0034395","regulation of transcription from RNA polymerase II promoter in response to iron",0.001,2.0058,0.780,0.173,"positive regulation of myosin II filament organization"),
c("GO:0098942","retrograde trans-synaptic signaling by trans-synaptic protein complex",0.000,2.0058,0.836,0.178,"positive regulation of myosin II filament organization"),
c("GO:0070935","3'-UTR-mediated mRNA stabilization",0.003,1.5329,0.825,0.179,"positive regulation of myosin II filament organization"),
c("GO:0051775","response to redox state",0.016,1.5329,0.895,0.182,"positive regulation of myosin II filament organization"),
c("GO:0023041","neuronal signal transduction",0.001,2.0058,0.780,0.205,"positive regulation of myosin II filament organization"),
c("GO:1990255","subsynaptic reticulum organization",0.000,2.0058,0.905,0.243,"positive regulation of myosin II filament organization"),
c("GO:0043974","histone H3-K27 acetylation",0.001,1.5329,0.890,0.290,"positive regulation of myosin II filament organization"),
c("GO:0007166","cell surface receptor signaling pathway",0.920,1.9289,0.684,0.308,"positive regulation of myosin II filament organization"),
c("GO:0043570","maintenance of DNA repeat elements",0.064,1.7069,0.864,0.314,"positive regulation of myosin II filament organization"),
c("GO:0043279","response to alkaloid",0.017,1.4055,0.878,0.335,"positive regulation of myosin II filament organization"),
c("GO:0072719","cellular response to cisplatin",0.001,1.4100,0.868,0.392,"positive regulation of myosin II filament organization"),
c("GO:0060052","neurofilament cytoskeleton organization",0.002,1.7069,0.824,0.401,"positive regulation of myosin II filament organization"),
c("GO:0007040","lysosome organization",0.019,1.6048,0.880,0.403,"positive regulation of myosin II filament organization"),
c("GO:0007188","adenylate cyclase-modulating G-protein coupled receptor signaling pathway",0.043,1.6382,0.735,0.408,"positive regulation of myosin II filament organization"),
c("GO:0048789","cytoskeletal matrix organization at active zone",0.000,1.5329,0.896,0.418,"positive regulation of myosin II filament organization"),
c("GO:0043619","regulation of transcription from RNA polymerase II promoter in response to oxidative stress",0.006,1.5329,0.756,0.430,"positive regulation of myosin II filament organization"),
c("GO:1990926","short-term synaptic potentiation",0.000,2.0058,0.723,0.455,"positive regulation of myosin II filament organization"),
c("GO:0071243","cellular response to arsenic-containing substance",0.004,1.3152,0.859,0.461,"positive regulation of myosin II filament organization"),
c("GO:0007215","glutamate receptor signaling pathway",0.096,1.4843,0.701,0.467,"positive regulation of myosin II filament organization"),
c("GO:0010646","regulation of cell communication",0.929,1.5408,0.757,0.521,"positive regulation of myosin II filament organization"),
c("GO:0023051","regulation of signaling",0.934,1.5529,0.766,0.536,"positive regulation of myosin II filament organization"),
c("GO:0007169","transmembrane receptor protein tyrosine kinase signaling pathway",0.167,1.3876,0.690,0.697,"positive regulation of myosin II filament organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_ASE_both.pdf", width=16, height=9 ) # width and height are in inches

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

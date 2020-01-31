

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
revigo.data <- rbind(c("GO:0006043","glucosamine catabolic process",0.000,1.3166,0.974,0.000,"glucosamine catabolism"),
c("GO:0006998","nuclear envelope organization",0.016,2.2224,0.856,0.000,"nuclear envelope organization"),
c("GO:0060260","regulation of transcription initiation from RNA polymerase II promoter",0.023,1.6576,0.755,0.284,"nuclear envelope organization"),
c("GO:0006366","transcription from RNA polymerase II promoter",1.430,1.3016,0.920,0.317,"nuclear envelope organization"),
c("GO:0030035","microspike assembly",0.001,1.3166,0.882,0.361,"nuclear envelope organization"),
c("GO:1905090","negative regulation of parkin-mediated mitophagy in response to mitochondrial depolarization",0.000,1.3166,0.681,0.381,"nuclear envelope organization"),
c("GO:1901741","positive regulation of myoblast fusion",0.003,1.3166,0.581,0.437,"nuclear envelope organization"),
c("GO:0042795","snRNA transcription from RNA polymerase II promoter",0.003,1.3835,0.916,0.448,"nuclear envelope organization"),
c("GO:0010606","positive regulation of cytoplasmic mRNA processing body assembly",0.003,1.5158,0.724,0.587,"nuclear envelope organization"),
c("GO:1903077","negative regulation of protein localization to plasma membrane",0.004,1.3166,0.629,0.660,"nuclear envelope organization"),
c("GO:0007608","sensory perception of smell",0.268,3.2671,0.681,0.000,"sensory perception of smell"),
c("GO:0007318","pole plasm protein localization",0.001,2.1732,0.613,0.405,"sensory perception of smell"),
c("GO:0001967","suckling behavior",0.002,1.3166,0.967,0.424,"sensory perception of smell"),
c("GO:0048875","chemical homeostasis within a tissue",0.002,1.3166,0.633,0.437,"sensory perception of smell"),
c("GO:0046552","photoreceptor cell fate commitment",0.003,2.7223,0.672,0.443,"sensory perception of smell"),
c("GO:0072160","nephron tubule epithelial cell differentiation",0.003,1.3166,0.658,0.448,"sensory perception of smell"),
c("GO:0034383","low-density lipoprotein particle clearance",0.004,1.3166,0.649,0.451,"sensory perception of smell"),
c("GO:2000036","regulation of stem cell population maintenance",0.005,1.8861,0.614,0.461,"sensory perception of smell"),
c("GO:0042066","perineurial glial growth",0.000,1.3166,0.695,0.467,"sensory perception of smell"),
c("GO:0032941","secretion by tissue",0.012,1.3835,0.668,0.489,"sensory perception of smell"),
c("GO:0051482","positive regulation of cytosolic calcium ion concentration involved in phospholipase C-activating G-protein coupled signaling pathway",0.004,1.3166,0.707,0.497,"sensory perception of smell"),
c("GO:0007129","synapsis",0.018,1.8861,0.769,0.505,"sensory perception of smell"),
c("GO:0002164","larval development",0.023,2.5296,0.667,0.511,"sensory perception of smell"),
c("GO:0035882","defecation rhythm",0.000,1.3166,0.748,0.520,"sensory perception of smell"),
c("GO:0021885","forebrain cell migration",0.013,1.8861,0.604,0.525,"sensory perception of smell"),
c("GO:0035107","appendage morphogenesis",0.055,1.5816,0.665,0.544,"sensory perception of smell"),
c("GO:0043403","skeletal muscle tissue regeneration",0.005,1.3166,0.682,0.558,"sensory perception of smell"),
c("GO:0003231","cardiac ventricle development",0.025,1.6576,0.655,0.559,"sensory perception of smell"),
c("GO:2000011","regulation of adaxial/abaxial pattern formation",0.001,1.3166,0.634,0.563,"sensory perception of smell"),
c("GO:0090188","negative regulation of pancreatic juice secretion",0.001,1.3166,0.514,0.574,"sensory perception of smell"),
c("GO:0060348","bone development",0.038,1.5682,0.654,0.596,"sensory perception of smell"),
c("GO:0002244","hematopoietic progenitor cell differentiation",0.026,1.3425,0.617,0.597,"sensory perception of smell"),
c("GO:0044770","cell cycle phase transition",0.207,1.6779,0.833,0.609,"sensory perception of smell"),
c("GO:0048569","post-embryonic animal organ development",0.036,1.3037,0.657,0.610,"sensory perception of smell"),
c("GO:0022605","oogenesis stage",0.001,1.3166,0.678,0.631,"sensory perception of smell"),
c("GO:0008406","gonad development",0.043,1.4711,0.619,0.639,"sensory perception of smell"),
c("GO:1901207","regulation of heart looping",0.000,1.3166,0.605,0.646,"sensory perception of smell"),
c("GO:0048729","tissue morphogenesis",0.180,1.3039,0.730,0.652,"sensory perception of smell"),
c("GO:0033058","directional locomotion",0.001,1.3166,0.968,0.000,"directional locomotion"),
c("GO:0043268","positive regulation of potassium ion transport",0.006,2.1732,0.678,0.000,"positive regulation of potassium ion transport"),
c("GO:1900044","regulation of protein K63-linked ubiquitination",0.001,1.3166,0.791,0.127,"positive regulation of potassium ion transport"),
c("GO:1902402","signal transduction involved in mitotic DNA damage checkpoint",0.008,2.1281,0.631,0.138,"positive regulation of potassium ion transport"),
c("GO:0051642","centrosome localization",0.007,1.5682,0.895,0.143,"positive regulation of potassium ion transport"),
c("GO:1903326","regulation of tRNA metabolic process",0.002,1.3166,0.823,0.173,"positive regulation of potassium ion transport"),
c("GO:0015931","nucleobase-containing compound transport",0.198,1.4338,0.896,0.180,"positive regulation of potassium ion transport"),
c("GO:0018002","N-terminal peptidyl-glutamic acid acetylation",0.001,1.3166,0.932,0.210,"positive regulation of potassium ion transport"),
c("GO:0007216","G-protein coupled glutamate receptor signaling pathway",0.005,1.3166,0.753,0.261,"positive regulation of potassium ion transport"),
c("GO:0086009","membrane repolarization",0.008,1.8861,0.747,0.290,"positive regulation of potassium ion transport"),
c("GO:0043401","steroid hormone mediated signaling pathway",0.102,1.7131,0.682,0.310,"positive regulation of potassium ion transport"),
c("GO:0048583","regulation of response to stimulus",1.120,1.3149,0.747,0.324,"positive regulation of potassium ion transport"),
c("GO:0071500","cellular response to nitrosative stress",0.002,1.3166,0.852,0.334,"positive regulation of potassium ion transport"),
c("GO:2000237","positive regulation of tRNA processing",0.000,1.3166,0.775,0.338,"positive regulation of potassium ion transport"),
c("GO:2000152","regulation of ubiquitin-specific protease activity",0.001,1.3166,0.810,0.364,"positive regulation of potassium ion transport"),
c("GO:1902946","protein localization to early endosome",0.001,1.3166,0.896,0.366,"positive regulation of potassium ion transport"),
c("GO:0009408","response to heat",0.166,1.3653,0.835,0.371,"positive regulation of potassium ion transport"),
c("GO:1902373","negative regulation of mRNA catabolic process",0.002,1.3835,0.764,0.390,"positive regulation of potassium ion transport"),
c("GO:0051235","maintenance of location",0.129,1.4882,0.756,0.406,"positive regulation of potassium ion transport"),
c("GO:0043489","RNA stabilization",0.008,1.3835,0.773,0.407,"positive regulation of potassium ion transport"),
c("GO:0042770","signal transduction in response to DNA damage",0.023,1.3425,0.721,0.410,"positive regulation of potassium ion transport"),
c("GO:1903019","negative regulation of glycoprotein metabolic process",0.003,1.3166,0.763,0.413,"positive regulation of potassium ion transport"),
c("GO:1990168","protein K33-linked deubiquitination",0.001,1.3166,0.924,0.415,"positive regulation of potassium ion transport"),
c("GO:0006408","snRNA export from nucleus",0.000,1.3166,0.850,0.424,"positive regulation of potassium ion transport"),
c("GO:0032057","negative regulation of translational initiation in response to stress",0.000,1.3166,0.716,0.449,"positive regulation of potassium ion transport"),
c("GO:0043123","positive regulation of I-kappaB kinase/NF-kappaB signaling",0.037,1.7562,0.633,0.452,"positive regulation of potassium ion transport"),
c("GO:0072344","rescue of stalled ribosome",0.001,1.3166,0.806,0.469,"positive regulation of potassium ion transport"),
c("GO:0046000","positive regulation of ecdysteroid secretion",0.000,1.3166,0.551,0.481,"positive regulation of potassium ion transport"),
c("GO:0015862","uridine transport",0.000,1.3166,0.877,0.497,"positive regulation of potassium ion transport"),
c("GO:1903781","positive regulation of cardiac conduction",0.000,1.3166,0.576,0.506,"positive regulation of potassium ion transport"),
c("GO:0010883","regulation of lipid storage",0.010,2.1281,0.714,0.517,"positive regulation of potassium ion transport"),
c("GO:0015801","aromatic amino acid transport",0.013,1.3166,0.848,0.523,"positive regulation of potassium ion transport"),
c("GO:0035666","TRIF-dependent toll-like receptor signaling pathway",0.000,1.3166,0.676,0.529,"positive regulation of potassium ion transport"),
c("GO:0035523","protein K29-linked deubiquitination",0.001,1.3166,0.923,0.537,"positive regulation of potassium ion transport"),
c("GO:0070327","thyroid hormone transport",0.004,1.3166,0.756,0.545,"positive regulation of potassium ion transport"),
c("GO:2000400","positive regulation of thymocyte aggregation",0.002,1.3166,0.711,0.554,"positive regulation of potassium ion transport"),
c("GO:0061393","positive regulation of transcription from RNA polymerase II promoter in response to osmotic stress",0.001,1.3166,0.678,0.556,"positive regulation of potassium ion transport"),
c("GO:0036245","cellular response to menadione",0.002,1.3166,0.850,0.562,"positive regulation of potassium ion transport"),
c("GO:0042327","positive regulation of phosphorylation",0.245,1.6525,0.696,0.586,"positive regulation of potassium ion transport"),
c("GO:2001137","positive regulation of endocytic recycling",0.001,1.3166,0.710,0.588,"positive regulation of potassium ion transport"),
c("GO:0086023","adrenergic receptor signaling pathway involved in heart process",0.001,1.3166,0.583,0.592,"positive regulation of potassium ion transport"),
c("GO:0034340","response to type I interferon",0.008,1.3166,0.806,0.597,"positive regulation of potassium ion transport"),
c("GO:0007089","traversing start control point of mitotic cell cycle",0.003,1.3166,0.689,0.598,"positive regulation of potassium ion transport"),
c("GO:2001234","negative regulation of apoptotic signaling pathway",0.042,1.6489,0.631,0.602,"positive regulation of potassium ion transport"),
c("GO:0001975","response to amphetamine",0.003,1.3166,0.861,0.609,"positive regulation of potassium ion transport"),
c("GO:0048387","negative regulation of retinoic acid receptor signaling pathway",0.002,1.3166,0.696,0.631,"positive regulation of potassium ion transport"),
c("GO:0048523","negative regulation of cellular process",1.830,1.3289,0.703,0.652,"positive regulation of potassium ion transport"),
c("GO:0035456","response to interferon-beta",0.004,1.3166,0.865,0.672,"positive regulation of potassium ion transport"),
c("GO:0099625","ventricular cardiac muscle cell membrane repolarization",0.004,1.8861,0.736,0.676,"positive regulation of potassium ion transport"),
c("GO:0046628","positive regulation of insulin receptor signaling pathway",0.002,1.3166,0.651,0.679,"positive regulation of potassium ion transport"),
c("GO:0098901","regulation of cardiac muscle cell action potential",0.005,1.6781,0.763,0.682,"positive regulation of potassium ion transport"),
c("GO:0033182","regulation of histone ubiquitination",0.005,1.3166,0.752,0.690,"positive regulation of potassium ion transport"),
c("GO:0033083","regulation of immature T cell proliferation",0.002,1.3166,0.780,0.692,"positive regulation of potassium ion transport"),
c("GO:0034764","positive regulation of transmembrane transport",0.023,2.1222,0.669,0.695,"positive regulation of potassium ion transport"),
c("GO:0008626","granzyme-mediated apoptotic signaling pathway",0.001,1.3166,0.750,0.695,"positive regulation of potassium ion transport"),
c("GO:0036150","phosphatidylserine acyl-chain remodeling",0.000,1.3166,0.903,0.085,"phosphatidylserine acyl-chain remodeling"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ASM_repro_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

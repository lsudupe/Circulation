### SCRIPT:  Read variables to check in different areas CIRCULATION proyect

## 01.10.23 Laura Sudupe , git @lsudupe

###### Libraries
library(openxlsx)
library("Seurat")
library("GSEABase")
library("ggplot2")
library(UCell)
library(RColorBrewer)
library(dplyr)
set.seed(123)


###### Read data
#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "Circulation/objects/area/enrich/")
DIR_RES <-  file.path(DIR_ROOT, "Circulation/results/new_area_features/")

combined <- readRDS("./Circulation/objects/area/new_areas.RDS")

# separate data
control <- subset(combined, subset = sample == "control")
dpi3 <- subset(combined, subset = sample == "dpi3")
dpi5_female <- subset(combined, subset = sample == "dpi5_female")
dpi5_male <- subset(combined, subset = sample == "dpi5_male")

# list
objects <- list(control = control, dpi3 = dpi3, dpi5_female = dpi5_female, dpi5_male = dpi5_male)

# features
###### Read the gene sets
#NEW
five.g <- c("Cthrc1", "Sfrp2", "Postn", "Fn1", "Ptn", "Col1a1", "Ddah1", 
"Tagln", "Acta2", "Lox", "Col1a2", "Col3a1", "Csrp2", "Col5a2", "Timp1", "Ctgf", "Ltbp2", "Tnc", 
"1500015O10Rik", "Comp", "Palld", "Tpm2", "Rflnb", "Mdk", "Thbs4", "Actn1", "Aspn", "Cilp", "Wisp2",
"Mfap4", "Wisp1", "Col8a1", "Thbs1", "Col5a1", "Sparc", "Tgfb3", "H19", "Fmod", "Loxl3", "Marcksl1",
"C1qtnf6", "Col14a1", "Col12a1", "Fstl1", "Lbh", "Vamp5", "Sfrp1", "Fkbp11", "Runx1", "Col16a1", "Cyr61",
"Myl9", "Serpinh1", "Cnn2", "Thbs2", "Ppic", "Rbp1", "Tspan6", "P4hb", "Mfap2", "Mmp14", "Fhl2", "Pdgfrl",
"Serf1", "Pmepa1", "Plac8", "Ssr2", "P4ha3", "Rrbp1", "Enpp1", "Sec61b", "Plod2", "Col8a2", "Fbn1", "Fkbp10", 
"Rcn1", "Crlf1", "Emilin1", "Cd200", "Lgals1", "Maged1", "Tnfrsf12a", "Synpo", "Foxp1", "Eef1b2", "Srpx2", 
"Serf2", "Ckap4", "Rcn3", "Sh3bgrl3", "Pdlim3", "Txndc5", "Vmp1", "Bgn", "Gpx7", "Dpysl3", "Efemp2", "Rasl11b",
"Loxl1", "Dok1", "Chpf", "Nrep", "Ppib", "Cdh2", "Sdc1", "Cpe", "Fat1", "Picalm", "Rpl12", "Clec11a", 
"Nr4a2", "Ttll5", "Rgs4", "Col11a1", "Phldb2", "Rpsa", "Rpl41", "Unc5b", "Tns3", "Myh9", "Glipr1", "Tcea3",
"P3h1", "Mical2", "Mlec", "Pxdn", "Siva1", "Ccl9", "Nedd9", "Kif26b", "Tmem176b", "Pdia5", "Fscn1",
"Sox4", "Ecscr", "Ddit4l", "Wwp2", "Rpl15", "Adarb1", "Kcnj15", "Bdnf", "Col27a1", "Dkk3", "2310022B05Rik",
"Ncam1", "Nkd2", "Ada", "Rplp0", "Cercam", "Odc1", "Rpl39", "Pdzrn3", "Rab31", "Neo1", "Vcl", 
"Hif1a", "4931406P16Rik", "Fbn2", "Dnm3os", "Tln2", "Gsto1", "Panx1", "Slc20a1", "Shisa4", "Rai14", "Acot7",
"Pafah1b3", "Cdh11", "Sertad4", "Pabpc4", "Adam12", "Kctd11", "Plscr2", "Gm15867", "Nek6", "Nt5dc2", 
"Rnf149", "Tmem98", "Cdkn2b", "Azin2", "Unc119", "Chsy1", "Piezo2", "Sipa1l1", "Galnt10", "Cx3cl1", "Nuak1", 
"Myo1d", "Uck2", "Fzd2", "Adgrl1", "Pycr1", "Gm17501", "Limd2", "Ctxn1", "Fmnl3", "Fgfrl1", "Myl12a", "Gxylt2",
"Ptprs", "Colgalt1", "Myl6", "Lhfpl2", "Tubb2b", "Ost4", "Frzb", "Olfml2b", "Tyms", "Slc39a14", "Hdlbp", "Apbb2", 
"Pabpc1", "Luzp1", "Surf4", "Smad7", "Fkbp14", "Ak1", "Dut", "Nbea", "Gng8", "Mpzl1", "Arf5", "Ass1", "Itga11", 
"Ubfd1", "Bmp1", "Spp1", "Slc1a3", "Rpl23", "Ssr4", "Cnn3", "Rack1", "Bsg", "Fam198b", "Mxra7", "Lrrc59", "Glrx5",
"Slc44a2", "Rpl27a", "Furin", "Tmsb4x", "Egr2", "Rpl28", "Ext1", "Gucy1a1", "Crtap", "1810058I24Rik", "Ltbp3", 
"Adamts2", "Tmem258", "Stmn1", "Rps20", "Aebp1", "Slc16a3", "Npm1", "Praf2", "Rps10", "Rgs3", "Rps15a", "Eef1g",
"Plk2", "Pls3", "Fndc1", "1110008F13Rik", "Arhgap31", "Erp29", "Rpl31", "Sdc3", "0610012G03Rik", "Itsn1",
"Slc38a10", "Sec61g", "Rps17", "Dcun1d5", "Itgb1", "Atp5g2", "Tmem176a", "Srm", "Eid1", "Plod3", "Rpl32")

eight.g <- c("Timp1", "Stmn1", "Cks2", "H2afz", "Spp1", "Cenpa", "Acta2", "Birc5", "Selenoh", "Cdc20", "Hmgb2", 
"Tpm2", "Ran", "Ccnb2", "Mif", "Ube2s", "Pclaf", "Cdca3", "Tagln", "Ube2c", "Hnrnpa1", "Cthrc1", "Cks1b", "Tubb6", 
"Tagln2", "Slc25a5", "Pgk1", "Rbm3", "Cdca8", "Actg1", "Eif5a", "Csrp2", "Nme1", "Fkbp11", "Lgals1", "Ranbp1", 
"Cxcl5", "Tpm4", "Spc24", "Pgam1", "Ccnb1", "Hmgn2", "Ldha", "Pkm", "Tpi1", "Cdk1", "Knstrn", "Dut", "Lrrc59", 
"Ppia", "Lockd", "Tpx2", "H2afx", "Cdkn3", "Ccl2", "Ccna2", "Tyms", "Ccl9", "Snrpf", "Eno1", "Cfl1", "Rpsa",
"Txn1", "Nap1l1", "Snrpd1", "Sec61b", "Rpl12", "Tk1", "Racgap1", "Smc4", "Tpm1", "Nucks1", "Tnc", "Pfn1", "Hspe1",
"Eef1g", "Eef1b2", "Ssr2", "Cnn2", "Ddx39", "Cenpw", "Prelid1", "Cenpm", "Tcp1", "Pimreg", "Fscn1", "Ckap4", "Lsm5", 
"Snrpe", "Cdk4", "Loxl3", "Cox5a", "Pdap1", "Cenpf", "Rps27l", "Tacc3", "Actn1", "Hmgb3", "Tnfrsf12a", "Lox", "Rplp1", 
"Eloc", "Banf1", "Ppa1", "Sec61g", "Ckap2", "Atp5g2", "Rpl41", "Rack1", "Hmmr", "Pbk", "Rrm2", "Ass1", "Rps10",
"Rps17", "Kif20a", "Top2a", "Plk1", "Tmem97", "Dtymk", "Anp32e", "Hint1", "Slc16a3", "Cenpe", "Myl6", "Mki67", 
"Dctpp1", "Wisp1", "Snrpa1", "Nudcd2", "Rrm1", "Smc2", "Nuf2", "Sdc1", "Pmf1", "Mad2l1", "Tuba1c", "Erh", "Rps3",
"Anapc15", "Mxd3", "Prc1", "Rpl15", "Kif22", "Ckap2l", "Ecscr", "Marcksl1", "Rpl32", "Rps13", "Kif23", "Spc25", 
"Haus7", "Rangap1", "Tubb2b", "Fhl2", "Gmnn", "Steap1", "Bcat1", "Acot7", "Palld", "Cmc2", "Gsto1", "Gins2",
"Nasp", "Asf1b", "Jpt2", "Ctps", "Timm8a1", "Hist1h4i", "Kif2c", "Thyn1", "Anln", "Glipr1", "Aurka", "Rfc5", 
"Mtap", "Dnajc9", "Cenpq", "Shmt2", "Ppih", "Depdc1a", "Lmnb1", "Bzw2", "Nudt5", "Pabpc4", "Asns", "Mcm5",
"Mthfd2", "Dbf4", "Ada", "Lig1", "Ska1", "Ripk3", "Ppil1", "Rfc4", "C330027C09Rik", "Diaph3", "Nek6", "Cd44",
"Shcbp1", "Cenpn", "Ndc80", "Haus4", "Nup37", "Cep55", "P4ha3", "Limd2", "Bub1b", "Rcc1", "Eno1b", "Mcm7", 
"Dlgap5", "Nt5dc2", "Uck2", "Adam12", "Trip13", "Psat1", "Hells", "Mycbp", "Runx1", "Dhfr", "Cenpp", "Cdc25c",
"Incenp", "Uhrf1", "Melk", "Cenpk", "Cenph", "Alg8", "Ptprn", "Ezh2", "Tcf19", "Mns1", "Tex30", "Mcm2",
"Psrc1", "Reep4", "Nudt21", "Pfkl", "Ubl4a", "Hspd1", "Tubb5", "Myl9", "Col12a1", "S100a11", "Npm1", "Rps27a")

eleven.g <- c("User", "Postn", "Cthrc1", "Col1a1", "Fn1", "Timp1", "Col3a1", 
              "Acta2", "Lgals1", "Csrp2", "Tagln", "Tpm2", "Col1a2", "Sparc", "Fstl1", "Ptn", 
              "Myl6", "Tmsb4x", "Lox", "Ltbp2", "Col5a2", "Tnc", "Cilp", "Serf2", "Rbp1", "Ctgf",
              "Rpl12", "Serpinh1", "Ddah1", "Rps17", "Rpsa", "Ppib", "Aspn", "Rpl31", "Rps26", 
              "Mfap4", "Cald1", "Igfbp7", "Rps20", "Rpl23", "Tagln2", "Calr", "Rps8", "Rpl19", 
              "S100a11", "Gm42418", "Rpl27a", "Rpl41", "Rpl28", "Mdk", "Actn1", "Palld", "Tmsb10", 
              "Rrbp1", "Rpl22", "Stmn1", "Col8a1", "Rpl39", "Npm1", "Bgn", "Loxl3", "Mfap5",
              "Myl12a", "Sec61g", "Elob", "Hsp90ab1", "Rps10", "Eln", "Arf5", "Fkbp11", "Eef1b2",
              "Rcn3", "Mmp14", "Snx7", "Fam96b", "Marcksl1", "Tpm4", "Pdia6", "Rps27l", "Col5a1", 
              "Tgfb3", "Mettl9", "Sec61b", "Mif", "Rps19", "Col16a1", "Sh3bgrl3", "Hspa5", "Mrpl18", 
              "Bet1", "Ccnl2", "Vdac1", "Tm9sf3", "Anxa2", "Meox1", "Cct2", "Actg1", "Enpp1",
              "Nupr1", "Kctd10", "Dazap1", "Ssr2", "Nptn", "Dap", "Cct8", "Ywhah", "Mgat2", 
              "Dst", "Bzw1", "Hsp90b1", "Rpl38", "Zfp207", "Mrpl51", "Hspd1", "Vasp", "Msn", 
              "Myl9", "Cct3", "Nars", "Tomm40", "Tubb6", "Gnai3", "Gcsh", "Tpm1", "Fads3",
              "Caprin1", "Plscr3", "Tial1", "Ufm1", "Eif4b", "Maged1", "Srsf6", "Fam104a",
              "Rnf10", "Fbl", "Nop56", "Eif2s1", "Ppic", "Tmem208", "Hspe1", "Lamb2", "Bax", 
              "Srsf7", "Hspa4", "Sin3b", "Atp1a1", "Raly", "Gpr180", "Parva", "Hnrnph1", 
              "Creld2", "Hnrnpf", "Eif3c", "Nono", "Crtap", "Ube2l3", "P3h4", "Eif3a", 
              "Ywhag", "Kdelr2", "Nudcd2", "Pofut2", "Tcp1", "Manf", "Mpzl1", "Psmd2", 
              "Cbx3", "Pdia4", "H2afz", "Atp5e", "Ubxn4", "Hnrnpab", "Ptbp1", "Colgalt1",
              "Iqgap1", "Cyr61", "Arap1", "Psmb5", "Syncrip", "Map1b", "Canx", "Uhrf2", 
              "Itgb5", "Eif5b", "Pdia3", "Sfpq", "Snrpd1", "Spcs3", "Dad1", "Ap2a2", "Kpnb1", 
              "Pole3", "Hist1h1c", "Arfgap1", "Ankrd11", "Gars", "Actn4", "Lman1", "Eif4g1",
              "Plod3", "Rnf181", "Surf4", "Sf3b3", "Hist1h2bc", "Olfml2b", "Cyb5r3", "Ptms", 
              "Idh2", "C1qtnf6", "Igf2r", "Atp5j2", "Eif4g2", "Cdk11b", "Fbn1", "Ost4",
              "Smc4", "Sdc1", "Hnrnpm", "Pcbp1", "Cdv3", "Ppia", "Adamts2", "Nap1l1", 
              "Fbln2", "Psmb4", "Adcy7", "Pdlim7", "Tm4sf1", "Bola2", "Txnl4a", "Thbs3", 
              "Khdrbs1", "Sh3glb1", "Ndufa3", "Mbd3", "Cox5a", "Wdr1", "Tardbp", "Tnfrsf12a",
              "Fxyd6", "Ptprs", "Sec61a1", "Copb2", "Tceal8", "Odc1", "Plec", "Ssr1", "Fhl2", 
              "Bsg", "Kdelr3", "Actr3", "Anp32b", "H2afj", "Ubl4a", "Ubap2l", "Ncl", "Gmppb", 
              "Arpp19", "Pabpc1", "Copg1", "Arhgef40", "Aldh1a2", "Mrpl20", "Pdgfrb", "Dctpp1", 
              "Taf1d", "Luzp1", "Gspt1", "Sfrp2", "Ranbp1", "Emilin1", "Aebp1", "Tubb5", "Eef1g", 
              "Ppa1", "Oaf", "Vapa", "Nme1", "Fat1", "Mrpl52", "Eny2", "Emp1", "Ab")

genes <- read.xlsx("./Circulation/data/Silvia_clean.xlsx", colNames = TRUE)
verde <- na.omit(genes$verde)
morado <- na.omit(genes$morado)

fb.genes <- read.csv("./Circulation/results/DE/top200_tokio_FB.csv")
fb.g <- na.omit(fb.genes$gene)

nine <- read.csv("./Circulation/results/genes/top100nine.csv")
nine.g <- na.omit(nine$x)

genes.velocity <- read.csv("./Circulation/data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("./Circulation/data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
DA <- as.vector(genes.velocity$Dynamics_1)
DB <- as.vector(genes.velocity$Dynamics_2)
DB[DB == ""] <- NA 
DB <- DB[!is.na(DB)]

enrich <- c()
###############################################333
for (i in 1:length(objects)){
  ##matrix
  a <- objects[[i]]
  DefaultAssay(a) <- "SCT"
  
  ## Check intersect genes B and object
  v.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], verde))
  m.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], morado))
  fb.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], fb.g))
  nine.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], nine.g))
  b.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], B))
  da.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], DA))
  db.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], DB))
  five.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], five.g))
  eight.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], eight.g))
  eleven.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], eleven.g))
  
  ## Add UCellScore
  a <- AddModuleScore_UCell(a, features = v.genes, name = "_verde")
  a <- AddModuleScore_UCell(a, features = m.genes, name = "_morado")
  a <- AddModuleScore_UCell(a, features = m.genes, name = "_fb")
  a <- AddModuleScore_UCell(a, features = nine.genes, name = "_nine")
  a <- AddModuleScore_UCell(a, features = b.genes, name = "_b")
  a <- AddModuleScore_UCell(a, features = da.genes, name = "_da")
  a <- AddModuleScore_UCell(a, features = db.genes, name = "_db")
  a <- AddModuleScore_UCell(a, features = db.genes, name = "_five")
  a <- AddModuleScore_UCell(a, features = db.genes, name = "_eight")
  a <- AddModuleScore_UCell(a, features = db.genes, name = "_eleven")
  
  
  ## Calculate ratio verde/morado
  v <- as.vector(a$signature_1_verde)
  m <- as.vector(a$signature_1_morado)
  ratio <- log10(m/v)
  is.na(ratio) <-sapply(ratio, is.infinite)
  ratio[is.na(ratio)] = 0
  a$ratio <- ratio
  
  ## Calculate ratio da/db
  da <- as.vector(a$signature_1_da)
  db <- as.vector(a$signature_1_db)
  ratio.dadb <- log10(da/db)
  is.na(ratio.dadb) <-sapply(ratio.dadb, is.infinite)
  ratio.dadb[is.na(ratio.dadb)] = 0
  a$ratio.dadb <- ratio.dadb
  
  ## Regress out FB values
  meta <- a@meta.data
  lm <- lm(meta$signature_1_b ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_b"]] <- residuals
  
  lm <- lm(meta$signature_1_da ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_da"]] <- residuals
  
  lm <- lm(meta$signature_1_db ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_db"]] <- residuals
  
  lm <- lm(meta$signature_1_nine ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_nine"]] <- residuals
  
  lm <- lm(meta$signature_1_verde ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_verde"]] <- residuals
  
  lm <- lm(meta$signature_1_morado ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_morado"]] <- residuals
  
  lm <- lm(meta$signature_1_five ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_five"]] <- residuals
  
  lm <- lm(meta$signature_1_eight ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_eight"]] <- residuals
  
  lm <- lm(meta$signature_1_eleven ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_eleven"]] <- residuals
  
  ## Save
  saveRDS(a,file = paste0(DIR_DATA, names(objects[i]),".enrich.rds"))
  enrich[[length(enrich) + 1]] <- a
  
}
names(enrich) <- c("control", "dpi3","dpi5_female","dpi5_male")

#####################seg.PLOTS#####################
## Create a custom color scale
#nombres <- c("IZ","BZ_1","BZ_2","BZ_3","BZ_2_3","BZ_4","RZ")
#colors <- c("coral4", "skyblue", "red", "darkgreen", "yellow","khaki2", "darkolivegreen", "gray89")
#names(colors)  <- nombres

nombres <- c("IZ","BZ1","BZ2","RZ")
colors <- c("red", "skyblue", "darkgreen", "yellow")
names(colors)  <- nombres

for (i in 1:length(enrich)){
  ## Object
  a <- enrich[[i]]
  
  ## Order area as factor
  a@meta.data[["new_area"]] <- factor(a@meta.data[["new_area"]],
                                  levels = c("IZ","BZ1","BZ2","RZ"),ordered = TRUE)
  
  ## boxplot
  # b
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_b.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_b, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_b_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_b, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  error <- "aaaaa"
  
  # da
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_da.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_da, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_da_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_da, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  # db
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_db.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_db, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_db_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_db, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  # nine
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_nine.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_nine, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_nine_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_nine, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  # verde
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_verde.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_verde, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_verde_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_verde, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  # morado
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_morado.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_morado, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_morado_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_morado, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_ratio.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=ratio, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  
  # five
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_five.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_five, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_five_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_five, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  # eight
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_eight.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_eight, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_eight_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_eight, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  
  # eleven
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_eleven.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=signature_1_eleven, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_eleven_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_eleven, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  
  # ratio
  pdf(file.path(paste0(DIR_RES, names(enrich[i]),"/", names(enrich[i]), "boxplot_ratio.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=ratio, y= new_area, fill=new_area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  
  
  ## spatial plots
  library(BuenColors)
  color <- jdb_palette("brewer_spectra")
  # b
  b <- c(min(a@meta.data[["signature_1_b"]]), max(a@meta.data[["signature_1_b"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_b"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_b.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_b"]]), max(a@meta.data[["residuals_b"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_b"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "sparial_enrich_b_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # da
  b <- c(min(a@meta.data[["signature_1_da"]]), max(a@meta.data[["signature_1_da"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_da"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_da.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_da"]]), max(a@meta.data[["residuals_da"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_da"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "sparial_enrich_da_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # db
  b <- c(min(a@meta.data[["signature_1_db"]]), max(a@meta.data[["signature_1_db"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_db"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_db.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_db"]]), max(a@meta.data[["residuals_db"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_da"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "sparial_enrich_db_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  
  # nine
  b <- c(min(a@meta.data[["signature_1_nine"]]), max(a@meta.data[["signature_1_nine"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_nine"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_nine.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_nine"]]), max(a@meta.data[["residuals_nine"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_nine"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "sparial_enrich_nine_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # verde
  b <- c(min(a@meta.data[["signature_1_verde"]]), max(a@meta.data[["signature_1_verde"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_verde"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_verde.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_verde"]]), max(a@meta.data[["residuals_verde"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_verde"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "sparial_enrich_verde_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # morado
  b <- c(min(a@meta.data[["signature_1_morado"]]), max(a@meta.data[["signature_1_morado"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_morado"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_morado.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_morado"]]), max(a@meta.data[["residuals_morado"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_morado"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_morado_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # five
  b <- c(min(a@meta.data[["signature_1_five"]]), max(a@meta.data[["signature_1_five"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_five"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_five.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_five"]]), max(a@meta.data[["residuals_five"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_five"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_five_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # eight
  b <- c(min(a@meta.data[["signature_1_eight"]]), max(a@meta.data[["signature_1_eight"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_eight"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_eight.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_eight"]]), max(a@meta.data[["residuals_eight"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_eight"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_eight_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # eleven
  b <- c(min(a@meta.data[["signature_1_eleven"]]), max(a@meta.data[["signature_1_eleven"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_eleven"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_eleven.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_eleven"]]), max(a@meta.data[["residuals_eleven"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_eleven"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_enrich_eleven_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # ratio
  b <- c(min(a@meta.data[["ratio"]]), max(a@meta.data[["ratio"]]))
  p1 <- SpatialFeaturePlot(a, features = c("ratio"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_ratio.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # ratio red and blue
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  
  l <- c(min(a@meta.data[["ratio"]]), max(a@meta.data[["ratio"]]))
  p1 <- SpatialFeaturePlot(a, features = c("ratio"), combine = FALSE, ncol = 2)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re),
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 na.value = "grey98",
                                 limits = l)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0(DIR_RES,names(enrich[i]), "/", paste(names(enrich[i]), "spatial_ratio_blue_red.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
}

###################################################
nombres <- c("IZ","BZ1","BZ2","RZ")
colors <- c("#000033", "#0066CC", "#33CC33","#FFF000")
names(colors) <- nombres

p <- SpatialDimPlot(combined, group.by = "new_area",  combine = TRUE, ncol = 2, cols = colors, alpha =  0.65)

pdf(paste0(DIR_RES, "spatial_new_areas.pdf"))
print(p)
dev.off()






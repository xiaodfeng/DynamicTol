# This file contains functions used for plotting
# The functions are sorted alphabetically

#' Bland-altman plot
#' @export
F_Bland <- function(x, y, name) {
  Stat <- blandr.statistics(x, y)
  DT <- data.table("Means" = Stat$means, "Differences" = Stat$differences)
  g_plot <- ggplot(DT, aes(x = Means, y = Differences)) +
    geom_point(size = 0.8, alpha = 0.7, color = "red") +
    xlim(0, qvalueCutoff) + # ylim(-0.6,0.6)+
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = Stat$upperLOA, linetype = "dashed", color = "green") +
    geom_hline(yintercept = Stat$bias, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = Stat$lowerLOA, linetype = "dashed", color = "brown") +
    labs(title = name) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"), plot.title = element_text(size = 17, face = "bold")
    )
  return(g_plot)
}

F_Distribution <- function(data) {
  data <- data[data > 0]
  weibull <- fitdist(data, "weibull", method = "mle") # for positive
  gamma <- fitdist(data, "gamma", method = "mle")
  norm <- fitdist(data, "norm", method = "mle")
  logis <- fitdist(data, "logis", method = "mle")
  plot.legend <- c("weibull", "gamma", "norm", "logis") # for positive
  # plot.legend <- c("gamma", "norm",'logis') #for neg
  List <- list(weibull, gamma, norm, logis) # for pos
  # List<-list(gamma, norm,logis) # for Neg
  par(mfrow = c(2, 2))
  denscomp(List,
    legendtext = plot.legend,
    xlab = "Scores",
    xlegend = "topleft"
  )
  cdfcomp(List,
    legendtext = plot.legend,
    xlab = "Scores",
    xlegend = "topleft"
  )
  qqcomp(List, legendtext = plot.legend)
  ppcomp(List, legendtext = plot.legend)
  # print(gofstat(List, fitnames = plot.legend))
  # print(summary(weibull))
  # print(summary(gamma))
  # print(summary(norm))
  # print(summary(logis))
}
F_DensityDotProduct <-function(Score, adjust = 2) {
  # F_DensityDotProduct(Score, adjust = 2) https://stackoverflow.com/questions/60879339/setting-width-and-height-of-a-single-panel-in-multi-panel-plot-in-cowplotplot
  Score.dp<-rbind(data.table('Score'=Score[,dpF],'Legend'='Fixed'),data.table('Score'=Score[,dpD],'Legend'='Dynamic'))
  g_dp <- ggplot(Score.dp)  + ylab('Count') + xlab('dot product') + theme_classic() + theme(legend.position="none") +
    geom_density(aes(x = Score, y = after_stat(count),colour = Legend), adjust = adjust) +
    scale_colour_manual( name="Legend",values = c("Fixed" = "Green","Dynamic" = "Blue"))
  
  Score.dp.decoy<-rbind(data.table('Score'=Score[,dpF.decoy],'Legend'='Fixed'),data.table('Score'=Score[,dpD.decoy],'Legend'='Dynamic'))
  g_dp.decoy <- ggplot(Score.dp.decoy)  + ylab('Count') + xlab('decoy') + theme_classic() + theme(legend.position="none") +
    geom_density(aes(x = Score, y = after_stat(count),colour = Legend), adjust = adjust) +
    scale_colour_manual( name="Legend",values = c("Fixed" = "Green","Dynamic" = "Blue"))
  
  Score.xcorr<-rbind(data.table('Score'=Score[,dpF.xcorr],'Legend'='Fixed'),data.table('Score'=Score[,dpD.xcorr],'Legend'='Dynamic'))
  g_dp.xcorr <- ggplot(Score.xcorr)  + ylab('Count') + xlab('xcorr') + theme_classic() + theme(legend.position="none") +
    geom_density(aes(x = Score, y = after_stat(count),colour = Legend), adjust = adjust) +
    scale_colour_manual( name="Legend",values = c("Fixed" = "Green","Dynamic" = "Blue"))
  # add legend
  pL <- ggplot(Score.xcorr)  +  geom_density(aes(x = Score, y = after_stat(count),colour = Legend))+
    scale_colour_manual( name="Legend",values = c("Fixed" = "Green","Dynamic" = "Blue"))
  l <- cowplot::get_legend(pL)
  cowplot::plot_grid(g_dp, g_dp.decoy, g_dp.xcorr,l,align="hv")
}
F_DensityScores <-function(DT, adjust = 2) {
  g_dpc <- ggplot(DT)  + ylab('Count') + xlab('dpc') + theme_classic() + theme(legend.position='top') +
    geom_density(aes(x = dpc, y = after_stat(count),colour = Legend), adjust = adjust) +
    scale_colour_manual(name="",values = ColorValues)
  # g_dpc
  g_dpc.decoy <- ggplot(DT)  + ylab('Count') + xlab('decoy') + theme_classic() + theme(legend.position='top') +
    geom_density(aes(x = dpc.decoy, y = after_stat(count),colour = Legend), adjust = adjust) +
    scale_colour_manual(name="",values = ColorValues)
  # g_dpc.decoy
  # g_dpc.xcorr <- ggplot(DT)  + ylab('Count') + xlab('dpc.xcorr') + theme_classic() + theme(legend.position="none") +
  #   geom_density(aes(x = dpc.xcorr, y = after_stat(count),colour = Legend), adjust = adjust) +
  #   scale_colour_manual(name="Legend",values = ColorValues)
  # g_dpc.xcorr
  # add legend
  # pL <- ggplot(DT)  +  geom_density(aes(x = dpc, y = after_stat(count),colour = Legend))+
  #   scale_colour_manual( name="Legend",values = ColorValues)
  # l <- get_legend(pL)
  print(plot_grid(g_dpc, g_dpc.decoy,align="v",nrow=2))
}

F_ExtractPlot <- function(MS2, PlotName) {
  name <- paste(PlotName, ".jpeg", seq = "")
  jpeg(file = name)
  ggplot(MS2, aes(mz, intensity)) +
    geom_col(colour = "black", size = 0.5) +
    labs(title = PlotName) #+xlim(0,1000)
  dev.off()
}
F_ExportPlotly <- function(fig, name = "plot", format = "svg") {
  # https://github.com/plotly/orca/issues/298#issue-582674958
  debug <- verbose <- safe <- F
  b <- plotly:::plotly_build(fig)
  plotlyjs <- plotly:::plotlyjsBundle(b)
  plotlyjs_path <- file.path(plotlyjs$src$file, plotlyjs$script)
  if (!is.null(plotlyjs$package)) {
    plotlyjs_path <- system.file(plotlyjs_path, package = plotlyjs$package)
  }
  tmp <- tempfile(fileext = ".json")
  cat(plotly:::to_JSON(b$x[c("data", "layout")]), file = tmp)
  args <- c(
    "graph", tmp, "-o", name, "--format",
    format, "--plotlyjs", plotlyjs_path, if (debug) "--debug",
    if (verbose) "--verbose", if (safe) "--safe-mode"
  )
  base::system(paste("orca", paste(args, collapse = " ")))
}
F_FDRCutoff <- function(DT, FDRCutoff, FileName, qvalueCutoff = 0.35) {
  print("Dynamic construction")
  TarDynamic <- DT[mztol == "NA" & Database == "Target"] %>% setnames(., "dpc", "Thres")
  DecDynamic <- DT[mztol == "NA" & Database == "Decoy"] %>% setnames(., "dpc.decoy", "Thres")
  Dynamic <- rbind(TarDynamic, DecDynamic, fill = TRUE)
  # Dynamic <- Dynamic[Thres!=1]
  print("Calculate the estimated qValue based on Pep accumulation")
  setorder(Dynamic, -Thres) # decreasing of the score according to dpc
  LambdaDynamic <- F_getPEPFromScoreLambda(
    Dynamic[Database == "Target"]$Thres,
    Dynamic[Database == "Decoy"]$Thres, paste0("FDR cutoff prediction DynamicCurve")
  )
  Dynamic[, pep := sapply(Dynamic$Thres, LambdaDynamic[[1]])]
  Dynamic[, SumPep := cumsum(pep)]
  Dynamic[, Length := seq_len(length(pep))]
  Dynamic[, FDR_Pep := SumPep / Length]
  Dynamic[, qValue_Pep := rev(cummin(rev(FDR_Pep)))]
  print("Calculate the estimated qValue based on target decoy")
  Dynamic[, decoy_hit := grepl("Decoy", Dynamic$Database)]
  Dynamic[, NDecoy := cumsum(decoy_hit)]
  Dynamic[, NTarget := Length - NDecoy]
  Dynamic[, FDR_DecoyTarget := NDecoy / NTarget]
  Dynamic[, qValue_DecoyTarget := rev(cummin(rev(FDR_DecoyTarget)))]
  print("Calculate the actual qValue based on inchkey 14")
  Dynamic[, NegInch := grepl("False", Dynamic$outcomeInch) & grepl("Target", Dynamic$Database)]
  Dynamic[, NNegInch := cumsum(NegInch)]
  Dynamic[, FDR_Inch := NNegInch / NTarget]
  Dynamic[, qValue_Inch := rev(cummin(rev(FDR_Inch)))]
  print("Calculate the actual qValue based on similarity matrix")
  Dynamic[, NegMatr := grepl("False", Dynamic$outcomeMatr) & grepl("Target", Dynamic$Database)]
  Dynamic[, NNegMatr := cumsum(NegMatr)]
  Dynamic[, FDR_Matr := NNegMatr / NTarget]
  Dynamic[, qValue_Matr := rev(cummin(rev(FDR_Matr)))]
  Colu <- c(
    "Thres", "outcomeInch", "outcomeMatr", "Database", "pep", "FDR_Pep", "qValue_Pep",
    "FDR_DecoyTarget", "qValue_DecoyTarget",
    "NNegInch", "NTarget", "FDR_Inch", "qValue_Inch", "FDR_Matr", "qValue_Matr",
    "SumPep", "Length", "decoy_hit", "NDecoy"
  )
  setcolorder(Dynamic, c(Colu, colnames(Dynamic)[!(colnames(Dynamic) %in% Colu)]))
  ## Calculate the identification rate
  DynamicUni <- unique(Dynamic[Database == "Target"], by = "query_inchikey14")
  DynamicUni[, NegInch := grepl("False", DynamicUni$outcomeInch)]
  DynamicUni[, NNegInch := cumsum(NegInch)]
  DynamicUni[, PosInch := grepl("True", DynamicUni$outcomeInch)]
  DynamicUni[, NPosInch := cumsum(PosInch)]
  DynamicUni[, RatioNegInch := NNegInch / 500]
  DynamicUni[, RatioPosInch := NPosInch / 500]
  # write.csv(Dynamic,'Dynamic.csv')
  ## Density plot
  Pep_Thres <- round(min(Dynamic[qValue_Pep < FDRCutoff]$Thres), digits = 2)
  print(paste("Pep_Thres", Pep_Thres))
  FDR_Thres <- round(min(Dynamic[qValue_DecoyTarget <= FDRCutoff]$Thres), digits = 2)
  print(paste("FDR_Thres", FDR_Thres))
  # Score.dpc<-rbind(data.table('Score'=Dynamic[,dpc],'Legend'='Target'),data.table('Score'=Dynamic[,dpc.decoy.Dec],'Legend'='Decoy'))
  svg(paste("Fig. 1e DensityTargetDecoyFDR", FDRCutoff, ".svg")) # ,width = 500, height = 500
  g_density <- ggplot(Dynamic) +
    ylab("Count") +
    xlab("Score") +
    theme_classic() +
    theme(legend.position = "top") +
    geom_density(aes(x = Thres, y = after_stat(count), colour = Database), adjust = 2) +
    # geom_density(data=Dynamic[Database=='Target'],aes(x = Thres, y = after_stat(count),colour = outcomeMatr), adjust = 2) +
    # geom_vline(xintercept = FDR_Thres, linetype="dotted",color = "orange", size=1.5)+
    # annotate(geom = "label", x = FDR_Thres, y = 1000, label = FDR_Thres,color = "orange")+
    geom_vline(xintercept = Pep_Thres, linetype = "dotted", color = "red", size = 1.5) +
    annotate(geom = "label", x = Pep_Thres, y = 1000, label = Pep_Thres, color = "red") +
    scale_colour_manual(name = "Legend", values = c("Target" = "blue", "Decoy" = "green"))
  print(g_density)
  dev.off()
  ## ThresVSQvalue
  ColorValues <<- c("Estimated q-value" = "red", "qValue_DecoyTarget" = "orange")
  svg(file = paste("Fig. S9 ThresVSQvalue", ".svg"))
  g_ThresVSQvalue <- ggplot(data = Dynamic, aes(x = Thres)) +
    # geom_point(aes(y=qValue_Matr,colour="True q-value"),size=1)+
    # geom_point(aes(y=qValue_InchFull,colour="qValue_InchFull"),size=1)+
    geom_point(aes(y = qValue_Pep, colour = "Estimated q-value"), size = 1) +
    geom_point(aes(y = qValue_DecoyTarget, colour = "qValue_DecoyTarget"), size = 1) +
    xlab("Thres") +
    ylab("Estimated q-value") + # labs(title="Peak detection mass tolerance in Da")+
    theme_classic() +
    theme(
      axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"), plot.title = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 12), legend.title = element_blank(),
      legend.justification = c(0.8, 0.3), legend.position = "top"
    ) +
    scale_colour_manual(name = "Legend", values = ColorValues)
  print(g_ThresVSQvalue)
  dev.off()
  ## Identification rate
  ColorValues <<- c("In library" = "purple", "Not in library" = "grey")
  svg(file = paste("Fig. S10 Identification rate", ".svg"))
  g_IdentificationRate <- ggplot(data = DynamicUni, aes(x = qValue_Pep)) +
    geom_point(aes(y = RatioPosInch, colour = "In library"), size = 1) +
    geom_point(aes(y = RatioNegInch, colour = "Not in library"), size = 1) +
    xlab("Estimated q-value") +
    ylab("Identification Rate") + # labs(title="Peak detection mass tolerance in Da")+
    theme_classic() +
    theme(
      axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"), plot.title = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 12), legend.title = element_blank(),
      legend.justification = c(0.8, 0.3), legend.position = "top"
    ) +
    scale_colour_manual(name = "Legend", values = ColorValues)
  print(g_IdentificationRate)
  dev.off()
  print(paste("qValue 0.05 NPosInch", max(DynamicUni[qValue_Pep < 0.05]$NPosInch)))
  print(paste("qValue 0.05 NNegInch", max(DynamicUni[qValue_Pep < 0.05]$NNegInch)))
  print(paste("qValue 0.01 NPosInch", max(DynamicUni[qValue_Pep < 0.01]$NPosInch)))
  print(paste("qValue 0.01 NNegInch", max(DynamicUni[qValue_Pep < 0.01]$NNegInch)))
  ## QvalueEstVSReal
  ColorValues <<- c("True q-value" = "grey", "Estimated q-value" = "red")
  svg(file = paste("Fig. 7a QvalueEstVSReal", ".svg"))
  DynamicCut <- Dynamic[qValue_Matr < qvalueCutoff]
  g_EstVSReal <- ggplot(data = DynamicCut, aes(x = qValue_Matr)) +
    geom_line(aes(y = qValue_Matr, colour = "True q-value"), size = 1) +
    # geom_point(aes(y=qValue_InchFull,colour="qValue_InchFull"),size=1)+
    geom_point(aes(y = qValue_Pep, colour = "Estimated q-value"), size = 1) +
    # geom_point(aes(y=qValue_DecoyTarget,colour="qValue_DecoyTarget"),size=1)+
    xlab("True q-value") +
    ylab("Estimated q-value") +
    theme_classic() +
    xlim(0, qvalueCutoff) +
    ylim(0, qvalueCutoff) +
    theme(
      axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"), plot.title = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 12), legend.title = element_blank(),
      legend.justification = c(0.8, 0.3), legend.position = "top"
    ) +
    scale_colour_manual(name = "Legend", values = ColorValues)
  print(g_EstVSReal)
  dev.off()
  ## Bland plot
  svg(file = paste("Fig. 7b Bland", ".svg"))
  
  g_Bland <- F_Bland(DynamicCut$qValue_Matr, DynamicCut$qValue_Pep, "True q-value - Estimated q-value")
  print(g_Bland)
  dev.off()
  ## Output plots together
  png(file = paste0("Combined FDR figures", FileName, ".png"), width = 1000, height = 1000)
  print(plot_grid(g_density, g_ThresVSQvalue, g_EstVSReal, g_Bland))
  # png(file = paste0('DensityQvaluePlots',FileName , '.png'),width = 1000,height = 500)
  # print(plot_grid(g_ThresVSQvalue,g_CountVSQvalue))
  dev.off()
  return(Dynamic)
}
F_gghistogram <- function(z, name, BinNum = 30) {
  ## Extract the scores seperately
  DynamicTarget <- z[mztol == "NA" & Database == "Target", ]
  DynamicDecoy <- z[mztol == "NA" & Database == "Decoy", ]
  F0.005Target <- z[mztol == 0.005 & Database == "Target", ]
  F0.005Decoy <- z[mztol == 0.005 & Database == "Decoy", ]
  F0.028Target <- z[mztol == 0.028 & Database == "Target", ]
  F0.028Decoy <- z[mztol == 0.028 & Database == "Decoy", ]
  F0.050Target <- z[mztol == 0.050 & Database == "Target", ]
  F0.050Decoy <- z[mztol == 0.050 & Database == "Decoy", ]
  png(paste0(name, "_His.png"), width = 1000, height = 500) # ,res = 300
  # layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 2, ncol = 4, byrow = TRUE))
  g_DynamicTarget <- gghistogram(DynamicTarget, ylim = c(0, 300), main = "DynamicTarget", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_F0.005Target <- gghistogram(F0.005Target, ylim = c(0, 300), main = "F0.005Target", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_F0.028Target <- gghistogram(F0.028Target, ylim = c(0, 300), main = "F0.028Target", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_F0.050Target <- gghistogram(F0.050Target, ylim = c(0, 300), main = "F0.050Target", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_DynamicDecoy <- gghistogram(DynamicDecoy, ylim = c(0, 300), main = "DynamicDecoy", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_F0.005Decoy <- gghistogram(F0.005Decoy, ylim = c(0, 300), main = "F0.005Decoy", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_F0.028Decoy <- gghistogram(F0.028Decoy, ylim = c(0, 300), main = "F0.028Decoy", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  g_F0.050Decoy <- gghistogram(F0.050Decoy, ylim = c(0, 300), main = "F0.050Decoy", x = "score", add = "mean", bins = BinNum, color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800"))
  print(plot_grid(nrow = 2, g_DynamicTarget, g_F0.005Target, g_F0.028Target, g_F0.050Target, g_DynamicDecoy, g_F0.005Decoy, g_F0.028Decoy, g_F0.050Decoy))
  dev.off()
}


F_PepRegression <- function(z, name, BinNum = 30) {
  ## Calculate the Pep score
  LambdaDynamic <- F_getPEPFromScoreLambda(z$dpD, z$dpD.decoy, paste0(name, "DynamicCurve"))
  z[, dpD.pep := sapply(z$dpD, LambdaDynamic[[1]])]
  z[, dpD.qval := sapply(z$dpD, LambdaDynamic[[2]])]
  LambdaFixed <- F_getPEPFromScoreLambda(z$dpF, z$dpF.decoy, paste0(name, "FixedCurve"))
  z[, dpF.pep := sapply(z$dpF, LambdaFixed[[1]])]
  z[, dpF.qval := sapply(z$dpF, LambdaFixed[[2]])]
  Lambda0.028 <- F_getPEPFromScoreLambda(z$dp0.028, z$dp0.028.decoy, paste0(name, "Fixed0.028Curve"))
  z[, dp0.028.pep := sapply(z$dp0.028, Lambda0.028[[1]])]
  z[, dp0.028.qval := sapply(z$dp0.028, Lambda0.028[[2]])]
  Lambda0.050 <- F_getPEPFromScoreLambda(z$dp0.050, z$dp0.050.decoy, paste0(name, "Fixed0.050Curve"))
  z[, dp0.050.pep := sapply(z$dp0.050, Lambda0.050[[1]])]
  z[, dp0.050.qval := sapply(z$dp0.050, Lambda0.050[[2]])]
  # z$dpD.pep.neg<-1-z$dpD.pep
  # z$dpF.pep.neg<-1-z$dpF.pep
  ## Histogram plot with mean lines
  # Change outline and fill colors by groups
  # Use custom color palette
  svg(file = paste0(name, "His_dpD", ".svg"))
  plot(gghistogram(z,
    x = "dpD", add = "mean", bins = BinNum,
    color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800")
  ))
  dev.off()

  svg(file = paste0(name, "His_dpF", ".svg"))
  plot(gghistogram(z,
    x = "dpF", add = "mean", bins = BinNum,
    color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800")
  ))
  dev.off()

  svg(file = paste0(name, "His_dpD_decoy", ".svg"))
  plot(gghistogram(z,
    x = "dpD.decoy", add = "mean", bins = BinNum,
    color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800")
  ))
  dev.off()

  svg(file = paste0(name, "His_dpF_decoy", ".svg"))
  plot(gghistogram(z,
    x = "dpF.decoy", add = "mean", bins = BinNum,
    color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800")
  ))
  dev.off()

  svg(file = paste0(name, "His_dp0.028", ".svg"))
  plot(gghistogram(z,
    x = "dp0.028", add = "mean", bins = BinNum,
    color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800")
  ))
  dev.off()
  svg(file = paste0(name, "His_dp0.050", ".svg"))
  plot(gghistogram(z,
    x = "dp0.050", add = "mean", bins = BinNum,
    color = "outcome", fill = "outcome", palette = c("#00AFBB", "#E7B800")
  ))
  dev.off()
  # F_PepRegression(z=TopN, name=paste0('Top',TopCutoff),BinNum=50)
}
F_PlotcombinePeaks <- function(MetaL, inchi, Single, Sensus, Color) {
  Selected <- MetaL[inchikey_14_precursor_mz == inchi]
  ## Extract MS2 spectra related to each meta id
  L <- list()
  for (j in 1:nrow(Selected)) {
    # j <- 1
    SelectedMS2 <- as.data.table(MSMS[library_spectra_meta_id == Selected[j, ]$id])
    SelectedMS2$normalized <- 100 * SelectedMS2$i / max(SelectedMS2$i)
    L[[j]] <- SelectedMS2[, c("mz", "normalized")] %>%
      setnames(., "mz", "intensity") %>%
      as.matrix(.)
  }
  Intersect <- Sensus[inchikey_14_precursor_mz == inchi]
  Intersect$i <- 100 * Intersect$i / max(Intersect$i)
  # Intersect <- combinePeaks(L,ppm = 10, peaks = 'intersect',minProp = 0.5)
  xlimV <- c(min(Intersect[, 1]) - 10, max(Intersect[, 1]) + 10)
  ylimV <- c(0, 100)
  # Union <- combinePeaks(L,ppm = 10, peaks = 'union',minProp = 0.5)
  svg(paste0("PlotcombinePeaks_", inchi, ".svg")) # ,width = 1000, height = 1100
  par(mfrow = c(3, 3), mar = c(5, 2, 2, 1)) #
  plot(L[[1]][, 1], L[[1]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 1", xlim = xlimV, ylim = ylimV)
  plot(L[[2]][, 1], L[[2]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 2", xlim = xlimV, ylim = ylimV)
  plot(L[[3]][, 1], L[[3]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 3", xlim = xlimV, ylim = ylimV)
  plot(L[[4]][, 1], L[[4]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 4", xlim = xlimV, ylim = ylimV)
  plot(L[[5]][, 1], L[[5]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 5", xlim = xlimV, ylim = ylimV)
  plot(L[[6]][, 1], L[[6]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 6", xlim = xlimV, ylim = ylimV)
  plot(L[[7]][, 1], L[[7]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 7", xlim = xlimV, ylim = ylimV)
  plot(L[[8]][, 1], L[[8]][, 2], type = "h", col = "black", xlab = "mz", ylab = "intensity", main = "Scan 8", xlim = xlimV, ylim = ylimV)
  # plot(L[[6]][, 1], L[[6]][, 2], type = "h", col = "black")
  # plot(Union[, 1], Union[, 2], type = "h", col = Color)
  # title(paste0(inchi,'_Union'))
  Intersect <- as.matrix(Intersect)
  plot(Intersect[, 1], Intersect[, 2], type = "h", col = Color, xlab = "mz", ylab = "intensity", main = "Consensus", xlim = xlimV, ylim = ylimV)
  # title(paste0(inchi,'_Consensus'))
  dev.off()
}

F_plot.roc <- function(z, name, BinNum = 30) {
  ## Extract the scores seperately
  DynamicTarget <- z[mztol == "NA" & Database == "Target", ]
  DynamicDecoy <- z[mztol == "NA" & Database == "Decoy", ]
  F0.005Target <- z[mztol == 0.005 & Database == "Target", ]
  F0.005Decoy <- z[mztol == 0.005 & Database == "Decoy", ]
  F0.028Target <- z[mztol == 0.028 & Database == "Target", ]
  F0.028Decoy <- z[mztol == 0.028 & Database == "Decoy", ]
  F0.050Target <- z[mztol == 0.050 & Database == "Target", ]
  F0.050Decoy <- z[mztol == 0.050 & Database == "Decoy", ]
  scorenames <- c("Dynamic", "F0.005", "F0.028", "F0.050")
  # colors <- RColorBrewer::brewer.pal(length(scorenames), "Dark2")
  # names(colors) <- scorenames
  ## Dynamic
  roc.Dynamic <- plot.roc(DynamicTarget$outcome, DynamicTarget$score,
    legacy.axes = T,
    col = ColorValues["Dynamic"], main = name
  )
  auc.Dynamic <- round(roc.Dynamic$auc, digits = 3)
  l.Dynamic <- paste("Dynamic AUC", auc.Dynamic)
  ## Fixed 0.005
  roc.F0.005 <- lines.roc(F0.005Target$outcome, F0.005Target$score, col = ColorValues["F0.005"])
  auc.F0.005 <- round(roc.F0.005$auc, digits = 3)
  l.F0.005 <- paste("F0.005 AUC", auc.F0.005)
  ## Fixed 0.028
  roc.F0.028 <- lines.roc(F0.028Target$outcome, F0.028Target$score, col = ColorValues["F0.028"])
  auc.F0.028 <- round(roc.F0.028$auc, digits = 3)
  l.F0.028 <- paste("F0.028 AUC", auc.F0.028)
  ## Fixed 0.050
  roc.F0.050 <- lines.roc(F0.050Target$outcome, F0.050Target$score, col = ColorValues["F0.050"])
  auc.F0.050 <- round(roc.F0.050$auc, digits = 3)
  l.F0.050 <- paste("F0.050 AUC", auc.F0.050)
  legend("bottomright",
    lwd = 2, col = ColorValues[scorenames],
    legend = c(l.Dynamic, l.F0.005, l.F0.028, l.F0.050)
  )
  Results <- data.table(
    "name" = name, "auc.F0.005" = auc.F0.005, "auc.F0.028" = auc.F0.028,
    "auc.F0.050" = auc.F0.050, "auc.Dynamic" = auc.Dynamic
  )
  return(Results)
}

F_PLotMir <- function(Top, Bottom, labelTitle, Query_inchkey, Library_inchkey, labelTop, labelBottom) {
  spec.top <- data.table("mz" = Top$mz, "intensity" = Top$i)
  spec.bottom <- data.table("mz" = Bottom$mz, "intensity" = Bottom$i)
  b <- 0
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 2])
  top_tmp$normalized <- round((top_tmp$intensity / max(top_tmp$intensity)) * 100, digits = 2)
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized) # data frame for plotting spectrum
  top <- subset(top_plot, top_plot$intensity > b) # data frame for similarity score calculation
  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[, 2])
  bottom_tmp$normalized <- round((bottom_tmp$intensity / max(bottom_tmp$intensity)) * 100, digits = 2)
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized) # data frame for plotting spectrum
  bottom <- subset(bottom_plot, bottom_plot$intensity > b) # data frame for similarity score calculation
  ## Dynamic matching
  alignment <- F_DynamicMatching(top, bottom)
  Score <- MSsim(dplyr::filter(alignment, intensity.bottom > 0)) %>% formatC(., format = "f", digits = 3) # round(.,digits = 4)
  print(Score)
  Equal <- setequal(top, bottom)
  ## plot the head to tail target
  svg(paste(
    "Equal", Equal, "Score", Score, "Precursor mz", labelTitle,
    "Query_inchkey", Query_inchkey, "Library_inchkey", Library_inchkey,
    ".svg"
  )) # ,width = 1000, height = 1100
  mzRange <- range(alignment$mz)
  # xlim <- c(mzRange[1] - 25, mzRange[2] + 25)
  xlim <- c(50, 150) # for plotting specific figure
  plot.new()
  plot.window(xlim = xlim, ylim = c(-125, 125))
  ticks <- c(-100, -50, 0, 50, 100)
  # COL <- adjustcolor(c("blue"), alpha.f = 1)
  for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]), col = "blue", lwd = 1) ## Add line
  for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]), col = "red")
  axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity")
  axis(1, pos = -125)
  lines(xlim, c(0, 0))
  rect(xlim[1], -125, xlim[2], 125)
  mtext("m/z", side = 1, line = 2)
  mtext("intensity (%)", side = 2, line = 2)
  title(paste("Score", Score, "Precursor m/z", labelTitle))
  text(mean(xlim), 100, labelTop)
  text(mean(xlim), -100, labelBottom)
  dev.off()
  return(Equal)
}
F_PLotMirOfftarget <- function(DT, MSMS) {
  for (id in c(1:nrow(DT))) {
    # id <- 1
    print(id)
    Query_inchkey <- DT[id, ]$query_inchikey14
    Query_precursor <- DT[id, ]$query_precursor_mz
    Query_id <- DT[id, ]$query_qpid
    Query_accession <- DT[id, ]$query_accession
    Query_name <- DT[id, ]$query_name
    Query_spectra <- MSMS[library_spectra_meta_id == Query_id]
    Library_inchkey <- DT[id, ]$library_inchikey14
    Library_precursor <- DT[id, ]$library_precursor_mz
    Library_id <- DT[id, ]$library_lpid
    Library_accession <- DT[id, ]$library_accession
    Library_name <- DT[id, ]$library_entry_name
    Library_spectra <- MSMS[library_spectra_meta_id == Library_id]
    EqualValue <- F_PLotMir(Query_spectra, Library_spectra,
      labelTitle = Query_precursor,
      Query_inchkey = Query_inchkey, Library_inchkey = Library_inchkey,
      labelTop = paste(Query_name, Query_inchkey, sep = " , "),
      labelBottom = paste(Library_name, Library_inchkey, sep = " , ")
    )
    print(EqualValue)
    DT[id, Equal := EqualValue]
  }
  return(DT)
}
F_PlotMSMS <- function(alignment, bottom_plot, top_plot) {
  ## generate plot
  ## calculation based on the common
  mzRange <- range(alignment$mz)
  xlim <- c(abs(mzRange[1] - 25), mzRange[2] + 25)
  # xlim <- c(1,100)
  plot.new()
  plot.window(xlim = xlim, ylim = c(-125, 125))
  ticks <- c(-100, -50, 0, 50, 100)
  for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]), col = "grey") ## Add line
  for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]), col = "grey")
  axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity", cex.axis = 1.5)
  axis(1, pos = -125, cex.axis = 1.5)
  lines(xlim, c(0, 0))
  rect(xlim[1], -125, xlim[2], 125)
  mtext("m/z", side = 1, line = 2, cex = 2)
  mtext("intensity (%)", side = 2, line = 2, cex = 2)
  text(mean(xlim), 100, "Query", col = "blue")
  text(mean(xlim), -100, "Library", col = "red")
  # text(mean(xlim), 100, paste("dp", round(dp, digits = 3), "rdp", round(rdp, digits = 3)))
  ## Add the mz value if there is a match
  ## The mz value in the alignment is from the bottom
  align_plot <- data.table(subset(alignment, alignment$intensity.top > 0))
  align_plot[, Topmz := top_plot$mz] #
  align_plot[, Da := formatC(abs(Topmz - mz), format = "e", digits = 2)] #
  align_plot[, Roundmz := formatC(mz, format = "f", digits = 2)]
  align_plot <- align_plot[intensity.bottom > 0] # only use the common for plotting the align
  if (length(align_plot$mz) > 0) {
    for (i in 1:length(align_plot$mz)) {
      text(align_plot$mz[i], align_plot$intensity.top[i] + 15, align_plot$Da[i], col = "blue", cex = 0.5, srt = 90)
      text(align_plot$mz[i], -align_plot$intensity.bottom[i] - 15, align_plot$Roundmz[i], col = "red", cex = 0.5, srt = 90)
      lines(rep(align_plot$mz[i], 2), c(0, align_plot$intensity.top[i]), col = "blue")
      lines(rep(align_plot$mz[i], 2), c(0, -align_plot$intensity.bottom[i]), col = "red")
    }
  }
}

F_PlotMSMSOriginal <- function(alignment, top_plot, bottom_plot, Lwd = 1, label, dp) {
  ## To plot the original MSMS head to tail figure, for demonstration purpose
  ## calculation based on the common
  mzRange <- range(alignment$mz)
  xlim <- c(mzRange[1] - 25, mzRange[2] + 25)
  plot.new()
  plot.window(xlim = xlim, ylim = c(-125, 125))
  ticks <- c(-100, -50, 0, 50, 100)
  # COL <- adjustcolor(c("blue"), alpha.f = 1)
  for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]), col = "blue", lwd = Lwd) ## Add line
  for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]), col = "red")
  axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity")
  axis(1, pos = -125)
  lines(xlim, c(0, 0))
  rect(xlim[1], -125, xlim[2], 125)
  mtext("m/z", side = 1, line = 2)
  mtext("intensity (%)", side = 2, line = 2)
  text(mean(xlim), 100, label)
  # text(mean(xlim), 80, paste("dp", formatC(dp, format = "e", digits = 3)))

  # svg(paste('Decoy Offset = dynamic range', '.svg')) #,width = 1000, height = 1100
  # F_PlotMSMSOriginal(alignment = AlignFixedDecoy.rbind,top_plot=decoy_plot.rbind,bottom_plot=bottom_plot,
  #                    Lwd=0.005, label = paste('Decoy Offset = dynamic range'),dp = dpFDecoy.rbind)
  # dev.off()
}



F_PlotMSMSXcorr <- function(alignment, top_plot, bottom_plot, label, dp) {
  ## To plot the original MSMS head to tail figure, for demonstration purpose
  ## calculation based on the common
  mzRange <- range(alignment$mz)
  xlim <- c(mzRange[1] - 25, mzRange[2] + 25)
  plot.new()
  plot.window(xlim = xlim, ylim = c(-125, 125))
  ticks <- c(-100, -50, 0, 50, 100)
  for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]), col = "blue") ## Add line
  for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]), col = "red")
  axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity")
  axis(1, pos = -125)
  lines(xlim, c(0, 0))
  rect(xlim[1], -125, xlim[2], 125)
  mtext("m/z", side = 1, line = 2)
  mtext("intensity (%)", side = 2, line = 2)
  text(mean(xlim), 100, label)
  text(mean(xlim), 80, paste("dp", round(dp, digits = 4)))
  ## Add the mz value if there is a match
  align_plot <- data.table(subset(alignment, alignment$intensity.top > 0))
  align_plot[, Topmz := top$mz]
  align_plot[, Da := formatC(abs(Topmz - mz), format = "e", digits = 2)]
  align_plot[, Roundmz := formatC(mz, format = "f", digits = 2)]
  align_plot <- align_plot[intensity.bottom > 0]
  if (length(align_plot$mz) > 0) {
    for (i in 1:length(align_plot$mz)) {
      text(align_plot$mz[i], -align_plot$intensity.bottom[i] - 15, align_plot$Roundmz[i], col = "red", cex = 0.5, srt = 90)
      text(align_plot$mz[i], align_plot$intensity.top[i] + 15, align_plot$Da[i], col = "blue", cex = 0.5, srt = 90)
    }
  }
}
F_ROC <- function(Sensus,Name){
  png(paste(Name,'Sensus Inch.VS.Matr xcorr.VS.dpc','ROC.png'),width = 1000, height = 1000) # ,res = 300
  layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
  # png(paste(Name,'Single.VS.Sensus Inch.VS.Matr xcorr','ROC.png'),width = 1500, height = 1500) # ,res = 300
  # layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
  # Single[,outcome:=outcomeInch] %>% .[,score:=dpc.xcorr] %>% F_plot.roc(z=., name=paste0('SingleInch'))
  # Single[,outcome:=outcomeMatr] %>% .[,score:=dpc.xcorr] %>% F_plot.roc(z=., name=paste0('SingleMatr'))
  Sensus[,outcome:=outcomeInch] %>% .[,score:=dpc.xcorr] %>% F_plot.roc(z=., name=paste0('SensusInchXcorr'))
  Sensus[,outcome:=outcomeMatr] %>% .[,score:=dpc.xcorr] %>% F_plot.roc(z=., name=paste0('SensusMatrXcorr'))
  Sensus[,outcome:=outcomeInch] %>% .[,score:=dpc] %>% F_plot.roc(z=., name=paste0('SensusInchDpc'))
  Sensus[,outcome:=outcomeMatr] %>% .[,score:=dpc] %>% F_plot.roc(z=., name=paste0('SensusMatrDpc'))
  dev.off()
}

F_simplot <- function(mat, FileName) {
  # matAbove <- mat[1:Cutoff,1:Cutoff]
  # View(matAbove)
  svg(paste("Heatmap similarity cutoff", Cutoff, FileName, ".svg", seq = "")) # ,width = 13,height=11
  plot(DOSE::simplot(mat, xlab = "InchiKey", ylab = "InchiKey", font.size = 3, labs = F, color.high = "red"))
  dev.off()
}
F_SpectrumLoop <- function(Spectra, FT, mztol) {
  # ScoreDT <- F_SpectrumLoop(Spectra = Spectra, FT = Test, mztol = 0.005)
  # ScoreDT <- F_SpectrumLoop(Spectra = Spectra, FT = filteredFeature, mztol = 0.005)
  # names(FT)
  ScoreDT <- apply(FT, 1, F_SpectrumSingle)
  ScoreDT <- data.table(plyr::ldply(ScoreDT, data.frame))
  return(ScoreDT)
}
F_SpectrumSingle <- function(x) {
  ft <- x[1]
  mz <- as.numeric(x[2])
  rt <- as.numeric(x[5])
  int <- as.numeric(x[14])
  # print(x)
  ex_FT <- Spectra[mcols(Spectra)$feature_id == ft]
  if (length(ex_FT@listData) < 2) {
    ScoreDT <- data.table()
  } else {
    mgf_FT1 <- data.table("mz" = ex_FT[[1]]@mz, "intensity" = ex_FT[[1]]@intensity)
    mgf_FT2 <- data.table("mz" = ex_FT[[2]]@mz, "intensity" = ex_FT[[2]]@intensity)
    ScoreDT <- F_SpectrumSimilarity(spec.top = mgf_FT1, spec.bottom = mgf_FT2, Tol = Tol, mz = mz, rt = rt, int = int, ft = ft)
  }
  return(ScoreDT)
}

F_SpectrumSimilarity <- function(spec.top, spec.bottom, mztol = 0.005, mz, rt, int, ft, b = 0, x.threshold = 0) {

  ## For testing purpose
  # feature_id <- "FT0003"
  # ex_FT <- Spectra[mcols(Spectra)$feature_id == feature_id]
  # ex_FT
  # spec.top <-data.table("mz" = ex_FT[[1]]@mz, "intensity" = ex_FT[[1]]@intensity) # mix1 as top
  # spec.bottom <-data.table("mz" = ex_FT[[2]]@mz, "intensity" = ex_FT[[2]]@intensity) # mix2 as bottom
  ## Target
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 2])
  top_tmp$normalized <- (top_tmp$intensity / max(top_tmp$intensity)) * 100
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized) # data frame for plotting spectrum
  top <- subset(top_plot, top_plot$intensity >= b) # data frame for similarity score calculation
  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[, 2])
  bottom_tmp$normalized <- (bottom_tmp$intensity / max(bottom_tmp$intensity)) * 100
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized) # data frame for plotting spectrum
  bottom <- subset(bottom_plot, bottom_plot$intensity >= b) # data frame for similarity score calculation
  AlignFixed <- F_FixedMatching(top = top, bottom = bottom)
  AlignDynamic <- F_DynamicMatching(top = top, bottom = bottom)
  dpF <- MSsim(AlignFixed[AlignFixed[, 3] > 0, ])
  # rdpF <- MSsim(AlignFixed[AlignFixed[, 3] > 0, ])
  dpD <- MSsim(AlignDynamic[AlignDynamic[, 3] > 0, ])
  # rdpD <- MSsim(AlignDynamic[AlignDynamic[, 3] > 0, ])
  ## Decoy
  decoy.rbind <- data.frame()
  PMMTmin <- F_CalPMMT(min(top$mz), MRP, RefMZ = 200)
  PMMTmax <- F_CalPMMT(max(top$mz), MRP, RefMZ = 200)
  HalfLenth <- 75
  for (id in c(seq(-HalfLenth, -1), seq(1, HalfLenth))) {
    decoy <- top
    decoy$mz <- top$mz + 1.5 * (id / abs(id)) * PMMTmax + id * PMMTmin
    # decoy$mz <- top$mz + ((id/abs(id))*PMMTmax + id*PMMTmin)/10 # for testing purpose
    decoy.rbind <- rbind(decoy.rbind, decoy)
  }
  setorder(decoy.rbind, mz)
  AlignFixed.decoy <- F_FixedMatching(top = decoy.rbind, bottom = bottom)
  AlignDynamic.decoy <- F_DynamicMatching(top = decoy.rbind, bottom = bottom)
  # write.csv(AlignFixed.decoy,'AlignFixed.decoy.csv')
  # write.csv(AlignDynamic.decoy,'AlignDynamic.decoy.csv')
  # write.csv(decoy.rbind,'decoy.rbind.csv')
  # MSsim(AlignFixed.decoy)
  # MSsim(AlignDynamic.decoy)
  # MSsim(AlignFixed.decoy[AlignFixed.decoy[, 3] > 0, ])
  # MSsim(AlignDynamic.decoy[AlignDynamic.decoy[, 3] > 0, ])
  ## Merge the decoy according to the unique mz of aligned results
  MergeFixed <- F_MergeDecoy(AlignFixed.decoy[AlignFixed.decoy[, 3] > 0, ])
  MergeDynamic <- F_MergeDecoy(AlignDynamic.decoy[AlignDynamic.decoy[, 3] > 0, ])
  ## Calculate the dot product and xcorr score
  dpF.decoy <- MSsim(MergeFixed)
  dpF.xcorr <- dpF - dpF.decoy
  # rdpF.decoy <- MSsim(AlignFixed.decoy[AlignFixed.decoy[, 3] > 0, ])

  dpD.decoy <- MSsim(MergeDynamic)
  dpD.xcorr <- dpD - dpD.decoy
  # rdpD.decoy <- MSsim(AlignDynamic.decoy[AlignDynamic.decoy[, 3] > 0, ])
  ## Plot
  png(paste(
    "xcorrD_F", round((dpD.xcorr - dpF.xcorr), digits = 3), "xcorrD", round(dpD.xcorr, digits = 3), "xcorrF", round(dpF.xcorr, digits = 3),
    "mz", round(mz, digits = 3), "rt", round(rt, digits = 0), "int", formatC(int, format = "e", digits = 2), mztol, "Da", ".png"
  ), width = 1000, height = 1100, res = 150)
  par(mfcol = c(2, 1), mar = c(4.5, 3.5, 1, 1)) #
  F_PlotMSMS(alignment = AlignFixed, top_plot = top_plot, top = top, bottom_plot = bottom_plot, bottom = bottom)
  title(paste("fixed tol, dp", round(dpF, digits = 3), "decoy", formatC(dpF.decoy, format = "e", digits = 3), "xcorr", round(dpF.xcorr, digits = 3)), outer = F)

  F_PlotMSMS(alignment = AlignDynamic, top_plot = top_plot, top = top, bottom_plot = bottom_plot, bottom = bottom)
  title(paste("dynamic tol, dp", round(dpD, digits = 3), "decoy", formatC(dpD.decoy, format = "e", digits = 3), "xcorr", round(dpD.xcorr, digits = 3)), outer = F)

  dev.off()
  ScoreDT <- data.table(
    "Row.names" = ft,
    "dpF" = dpF, "dpF.decoy" = dpF.decoy, "dpF.xcorr" = dpF.xcorr,
    "dpD" = dpD, "dpD.decoy" = dpD.decoy, "dpD.xcorr" = dpD.xcorr
  )
  return(ScoreDT)
}


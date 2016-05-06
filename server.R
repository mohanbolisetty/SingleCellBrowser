
# LIBRARY -----------------------------------------------------------------

library(shiny)
library(ggplot2)
library(Rtsne)
library(DT)
library(edgeR)
library(pheatmap)
#library(vioplot)
library(threejs)
library(plotly)
library(sm)

shinyServer(function(input, output) {
  
  set.seed(1)


# INPUTS ------------------------------------------------------------------

  readqc<-read.csv('Data/QC_PASS.csv',row.names=1,stringsAsFactors = TRUE,as.is=T)
  counts<-read.csv('Data/Counts.csv',row.names=1,stringsAsFactors = TRUE, as.is=T, check.names=F)
  featuredata<-read.csv('Data/FeatureData.csv',row.names=1,stringsAsFactors = TRUE, as.is=T)
  
  libsize<-colSums(counts)
  cpm<-1e6*(sweep(counts,2,libsize,"/"))
  log2cpm<-log2(cpm+1)

  saved_plots_and_tables <- reactiveValues(pca_plt = NULL,
                                           heatmap = NULL,
                                           maplt = NULL,
                                           genes_plt = NULL,
                                           qc_plt = NULL,
                                           violin_plt = NULL,
                                           volcano_plt = NULL,
                                           volcano_table = NULL,
                                           mv_plt = NULL,
                                           scatter_plt = NULL,
                                           scatter_var_plt = NULL,
                                           scatter_table = NULL,
                                           qc_plt = NULL,
                                           cond_dens_plt = NULL,
                                           samp_dens_plt = NULL,
                                           sample_table = NULL,
                                           kallisto_table = NULL,
                                           hm_plt = NULL,
                                           bs_var_plt = NULL,
                                           grp1=NULL,
                                           grp2=NULL)


# REACTIVE ----------------------------------------------------------------
  react_qc<-reactive({
    readqc_values<-readqc
    readqc_values$PASS.FAIL<-'FAIL'
    readqc_values[which(readqc_values$Genes_Detected > input$genes & 
                        readqc_values$Percent.Mitochondria < input$mitocontent &
                          readqc_values$Trimmed_filtered_reads > 300000 &
                          readqc_values$Percent_Mapping > 10.),'PASS.FAIL']<-'PASS'
    return(readqc_values)
    })
  
  filteredValues <- reactive({
    readqc<-react_qc()
    pass_samples<-rownames(readqc[which(readqc$PASS.FAIL=='PASS'),])
    filtered_genes<-log2cpm[rowSums(log2cpm > input$cpm_threshold) > input$cells_threshold,]
    filtered_genes<-filtered_genes[,pass_samples]
    return(filtered_genes)
  })
  
  filteredCounts <- reactive({
    filtered_genes<-filteredValues()
    counts<-counts[rownames(filtered_genes),colnames(filtered_genes)]
    return(counts)
  })


# QC_PLOTS ----------------------------------------------------------------

      output$qc_plot <- renderPlot({
    values<-react_qc()
    if (input$hm_qc == 'Total Reads'){
      p1<-ggplot(values,aes(x=1:nrow(values),y=Total_reads))+
        geom_bar(stat='identity')+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="top",
              legend.text=element_text(size=14),
              legend.title = element_blank())+
        xlab('Samples')
      saved_plots_and_tables$qc_plt<-p1
      p1
    }
    
    else if (input$hm_qc == 'Mapped Reads'){
      p1<-ggplot(values,aes(x=1:nrow(values),y=Mapped_reads))+
        geom_bar(stat='identity')+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="top",
              legend.text=element_text(size=14),
              legend.title = element_blank())+
        xlab('Samples')
      saved_plots_and_tables$qc_plt<-p1
      p1
    }
    else if (input$hm_qc == 'Duplication Rates'){
      p1<-ggplot(values,aes(x=1:nrow(values),y=Percent_Duplication))+
        geom_bar(stat='identity')+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="top",
              legend.text=element_text(size=14),
              legend.title = element_blank())+
        xlab('Samples')
      saved_plots_and_tables$qc_plt<-p1
      p1
    }
    else if (input$hm_qc == 'Percent Target'){
      p1<-ggplot(values,aes(x=1:nrow(values),y=Percent_Target))+
        geom_bar(stat='identity')+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="top",
              legend.text=element_text(size=14),
              legend.title = element_blank())+
        xlab('Samples')
      saved_plots_and_tables$qc_plt<-p1
      p1
    }
    
  })
  

# SECONDARY_QC ------------------------------------------------------------

  output$genes_plot <- renderPlot({
    values<-react_qc()
    p1<-ggplot(values,aes(x=Genes_Detected,y=Percent.Mitochondria))+
      geom_point(aes(color=factor(PASS.FAIL)))+
      theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
            axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
            strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16),
            legend.position="top",
            legend.text=element_text(size=14),
            legend.title = element_blank())+
      xlab('Genes Detected')+
      ylab('Percent Mitochondria')
    saved_plots_and_tables$genes_plt<-p1
    p1
    })

# PCA_TSNE ----------------------------------------------------------------

  pca.values<-reactive({
    filtereddata<-filteredValues()
    filtered.pca<-prcomp(t(filtereddata))
    return(filtered.pca)
  })
  
  tsne.values<-reactive({
    filtereddata<-filteredValues()
    rtsne_in<-as.matrix(t(filtereddata))
    if (nrow(rtsne_in) -1  < 3 * 30){
      rtsne_out<-Rtsne(rtsne_in,dims=3,perplexity = 10)
      return(as.data.frame(rtsne_out$Y))
    }
    else{
      rtsne_out<-Rtsne(rtsne_in,dims=3)
      return(as.data.frame(rtsne_out$Y))
    }
  })
  
  output$pca_plt<-renderPlot({
    pca.data<-pca.values()
    pca.frame<-as.data.frame(pca.data$x)
    x<-colnames(pca.frame)[as.integer(input$pc_x)]
    y<-colnames(pca.frame)[as.integer(input$pc_y)]
    
    p1<-ggplot(pca.frame,aes_string(x=x,y=y))+
      geom_point(aes(color='red'))+
      theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
            axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
            strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16),
            legend.position="none")+
      xlab(paste('PC',input$pc_x,sep=''))+
      ylab(paste('PC',input$pc_y,sep=''))
    saved_plots_and_tables$pca_plt<-p1
    p1
  })
  
  output$plt_pc_var<-renderPlot({
    pca.data<-pca.values()
    loadings<-as.data.frame(pca.data$rotation)
    if (input$sample == ''){
      plot(pca.data,type='b',main="PCA",xlab='Principal Components',
           cex.names=2)
    }
    else if(!is.null(input$sample)){
      ensgid<-rownames(featuredata[which(featuredata$Associated.Gene.Name==input$sample),])[1]
      barplot(as.matrix(loadings[ensgid,1:as.integer(input$pc_count)]),
              xlab=NULL)
    }
  })
  
  output$plt_pc_loadings<-renderPlot({
    pca.data<-pca.values()
    loadings<-as.data.frame(pca.data$rotation)
    pc_value=paste('PC',input$pc_input,sep='')
    sorteddata<-loadings[order(loadings[pc_value],decreasing=TRUE),]
    subset<-sorteddata[1:10,pc_value]
    barplot(subset,names.arg=featuredata[rownames(sorteddata)[1:10],'Associated.Gene.Name'],cex.names=0.8)
  })
  
  output$tsne_plt<-renderPlotly({
    tsne.data<-tsne.values()
    p <- plot_ly(tsne.data, x = V1, y = V2, z = V3, type = "scatter3d",
                 mode = 'markers',marker=list(size=3,symbol='circle',line=list(width=0))
    )
    layout(p)
    
#     dim_x<-paste('V',input$dim1,sep='')
#     dim_y<-paste('V',input$dim2,sep='')
# 
#      p1<-ggplot(tsne.data,aes_string(x=dim_x,y=dim_y))+
#        geom_point(aes(color='red'))+
#        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
#              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
#              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
#              axis.title.y = element_text(face="bold", size=16),
#              legend.position="none")+
#        xlab(paste('T-SNE',input$dim1,sep=''))+
#        ylab(paste('T-SNE',input$dim2,sep=''))
#      p1
  })

  
  
  
  
# DGE_SELECT_GROUPS ------------------------------------------------------------

    output$dge_plot1<-renderPlot({
    if (input$cluster == 'PCA'){
      pca.data<-pca.values()
      pca.frame<-as.data.frame(pca.data$x)
      x<-colnames(pca.frame)[as.integer(input$dge_dim1)]
      y<-colnames(pca.frame)[as.integer(input$dge_dim2)]
      
      p1<-ggplot(pca.frame,aes_string(x=x,y=y))+
        geom_point(aes(color='red'))+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="none")+
        xlab(paste('PC',input$dge_dim1,sep=''))+
        ylab(paste('PC',input$dge_dim2,sep=''))
      p1
    }
    else if(input$cluster == 'TSNE'){
      tsne.data<-tsne.values()
      dim_x<-paste('V',input$dge_dim1,sep='')
      dim_y<-paste('V',input$dge_dim2,sep='')
      
      p1<-ggplot(tsne.data,aes_string(x=dim_x,y=dim_y))+
        geom_point(aes(color='red'))+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="none")+
        xlab(paste('T-SNE',input$dge_dim1,sep=''))+
        ylab(paste('T-SNE',input$dge_dim2,sep=''))
      p1
    }
  })
  
  output$cluster1 <- DT::renderDataTable({
    if (input$cluster == 'PCA'){
      pca.data<-pca.values()
      pca.frame<-as.data.frame(pca.data$x)
      DT::datatable(brushedPoints(pca.frame[,1:5],input$b1),options = list(orderClasses = TRUE,lengthMenu = c(5, 30, 50), pageLength = 10))
    }
    else if(input$cluster == 'TSNE'){
      tsne.data<-tsne.values()
      DT::datatable(brushedPoints(tsne.data[,1:3],input$b1),options = list(orderClasses = TRUE,lengthMenu = c(5, 30, 50), pageLength = 5))
    }
  })  

  output$dge_plot2<-renderPlot({
    if (input$cluster == 'PCA'){
      pca.data<-pca.values()
      pca.frame<-as.data.frame(pca.data$x)
      x<-colnames(pca.frame)[as.integer(input$dge_dim1)]
      y<-colnames(pca.frame)[as.integer(input$dge_dim2)]
      
      p1<-ggplot(pca.frame,aes_string(x=x,y=y))+
        geom_point(aes(color='red'))+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="none")+
        xlab(paste('PC',input$dge_dim1,sep=''))+
        ylab(paste('PC',input$dge_dim2,sep=''))
      p1
    }
    else if(input$cluster == 'TSNE'){
      tsne.data<-tsne.values()
      dim_x<-paste('V',input$dge_dim1,sep='')
      dim_y<-paste('V',input$dge_dim2,sep='')
      
      p1<-ggplot(tsne.data,aes_string(x=dim_x,y=dim_y))+
        geom_point(aes(color='red'))+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="none")+
        xlab(paste('T-SNE',input$dge_dim1,sep=''))+
        ylab(paste('T-SNE',input$dge_dim2,sep=''))
      p1
    }
  })
  
  output$cluster2 <- DT::renderDataTable({
    if (input$cluster == 'PCA'){
      pca.data<-pca.values()
      pca.frame<-as.data.frame(pca.data$x)
      DT::datatable(brushedPoints(pca.frame[,1:5],input$b2),options = list(orderClasses = TRUE,lengthMenu = c(5, 30, 50), pageLength = 10))
    }
    else if(input$cluster == 'TSNE'){
      tsne.data<-tsne.values()
      DT::datatable(brushedPoints(tsne.data[,1:3],input$b2),options = list(orderClasses = TRUE,lengthMenu = c(5, 30, 50), pageLength = 5))
    }
  })  

# DGE_CALCULATIONS --------------------------------------------------------
    
  dge<-reactive({
    data<-filteredCounts()
    if (input$cluster == 'PCA'){
      pca.data<-pca.values()
      pca.frame<-as.data.frame(pca.data$x)
      
      grp1<-rownames(brushedPoints(pca.frame[,1:5],input$b1))
      grp2<-rownames(brushedPoints(pca.frame[,1:5],input$b2))

      saved_plots_and_tables$grp1<-grp1
      saved_plots_and_tables$grp2<-grp2
      
      if (length(grp1)>1 & length(grp2)>1){
        cnts<-cbind(data[,grp1],data[,grp2])
        grp.ids<-c(rep(0,length(grp1)),rep(1,length(grp2)))

        cds <- DGEList(counts = cnts, group = factor(grp.ids,levels=c(0,1)))
        cds <- estimateCommonDisp(cds)
        cds <- estimateTagwiseDisp(cds)
        de.tgw <- exactTest(cds)
        diff.dat <- topTags(de.tgw,n=nrow(cds$counts))
        res <- diff.dat$table
        write.csv(res,file='Data/diffGenes.csv')
        return(res)
      }
    }
    else if(input$cluster == 'TSNE'){
      tsne.data<-tsne.values()
      
      grp1<-rownames(brushedPoints(tsne.data[,1:3],input$b1))
      grp2<-rownames(brushedPoints(tsne.data[,1:3],input$b2))
      
      saved_plots_and_tables$grp1<-grp1
      saved_plots_and_tables$grp2<-grp2
      
      cnts<-cbind(data[,grp1],data[,grp2])
      grp.ids<-c(rep(0,length(grp1)),rep(1,length(grp2)))
      if (length(grp1)>1 & length(grp2)>1){
        cds <- DGEList(counts = cnts, group = factor(grp.ids,levels=c(0,1)))
        cds <- estimateCommonDisp(cds)
        cds <- estimateTagwiseDisp(cds)
        de.tgw <- exactTest(cds)
        diff.dat <- topTags(de.tgw,n=nrow(cds$counts))
        res <- diff.dat$table
        write.csv(res,file='Data/diffGenes.csv')
        return(res)
      }
      }
    
    })
  
  output$dge<-DT::renderDataTable({
    top.genes<-dge()
    top.genes$Associated.Gene.Name<-featuredata[rownames(top.genes),'Associated.Gene.Name']
    if (dim(top.genes)[1]>1){
      DT::datatable(top.genes, options = list(orderClasses = TRUE,lengthMenu = c(10, 30, 50), pageLength = 10))
    }
  })

# MA_PLOT -----------------------------------------------------------------
  output$maplot<-renderPlot({
    top.genes<-read.csv('Data/diffGenes.csv',row.names=1)
    p1<-ggplot(top.genes,aes(x=logCPM,y=logFC,colour = FDR<0.05))+
      geom_point()+
      scale_colour_manual(name = 'FDR<0.05', values = setNames(c('#e41a1c','#4daf4a'),c(T, F))) +
      theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
            axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
            strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16),
            legend.position="top",
            legend.text=element_text(size=14),
            legend.title=element_text(face='bold',size=14))+
      xlab('log2(CPM)')+
      ylab('log2(FC)')
    p1
  })
  
  output$viewdge<-renderPrint({
    top.genes<-read.csv('Data/diffGenes.csv',row.names=1)
    top.genes$Associated.Gene.Name<-featuredata[rownames(top.genes),'Associated.Gene.Name']
    nearPoints(top.genes,input$maplot_click)
  })

# HEATMAP -----------------------------------------------------------------
  
  output$heatmap<-renderPlot({
    top.genes<-read.csv('Data/diffGenes.csv',row.names=1,as.is=T)
    log2cpm_filtered<-filteredValues()
    
    grp1<-saved_plots_and_tables$grp1
    grp2<-saved_plots_and_tables$grp2
    
    #cat(stderr(),'test',grp1)
    
    samples<-c(grp1,grp2)
    #sample_color<-c(rep(0,length(grp1)),rep(1,length(grp2)))
    
    plotdata<-log2cpm[
      rownames(top.genes[which(top.genes$FDR<input$fdr_heatmap),]),
      samples]
    
    h<-pheatmap(as.matrix(plotdata),cluster_rows=TRUE,cluster_cols=TRUE,
                scale='row',fontsize_row=6,labels_col = colnames(plotdata),
                show_rownames = FALSE,clustering_distance_rows='correlation',
                clustering_distance_cols='correlation'
                )
    saved_plots_and_tables$heatmap<-h
    h
  })

# VIOLIN_PLOTS ------------------------------------------------------------  
  
  output$violin<-renderPlot({
    log2cpm_filtered<-filteredValues()
    
    grp1<-saved_plots_and_tables$grp1
    grp2<-saved_plots_and_tables$grp2
    
    geneid<-rownames(
      featuredata[which(featuredata$Associated.Gene.Name==input$violin_gene),]
      )[1]
    
    if(!is.null(grp1) & !is.null(grp2) & !is.null(geneid)){
      my.vioplot(t(log2cpm_filtered[geneid,grp1]),t(log2cpm_filtered[geneid,grp2]),
                  names=c('Grp1','Grp2'),col=c("#e41a1c", "#377eb8"))
      title(main=input$violin_gene,ylab='log2CPM')
    }
  })
  
  
  
  
# DOWNLOADS ---------------------------------------------------------------
    output$download_heatmap_plt <- downloadHandler(
    filename = function() { "heatmap.pdf" },
    content = function(file) {
      ggsave(file, saved_plots_and_tables$heatmap, width = 8, height = 10,dpi=300)
  })

  output$download_genes_plt <- downloadHandler(
    filename = function() { "genes_vs_mito.pdf" },
    content = function(file) {
      ggsave(file, saved_plots_and_tables$genes_plt, width = 8, height = 10,dpi=300)
    })
  
  output$download_qc_plt <- downloadHandler(
    filename = function() { paste(input$hm_qc,"QC_Plot.pdf" ,sep='_')},
    content = function(file) {
      ggsave(file, saved_plots_and_tables$qc_plt, width = 8, height = 10,dpi=300)
    })

  output$download_pca_plt <- downloadHandler(
    filename = function() {"PCA_Plot.pdf"},
    content = function(file) {
      ggsave(file, saved_plots_and_tables$pca_plt, width = 8, height = 10,dpi=300)
    })
  
  output$download_violin_plt <- downloadHandler(
    filename = function() {paste(input$violin_gene,"Plot.pdf",sep='_')},
    content = function(file) {
      pdf(file, width = 8, height = 10)
      
      log2cpm_filtered<-filteredValues()
      grp1<-saved_plots_and_tables$grp1
      grp2<-saved_plots_and_tables$grp2
      geneid<-rownames(
        featuredata[which(featuredata$Associated.Gene.Name==input$violin_gene),]
      )[1]
      
      if(!is.null(grp1) & !is.null(grp2) & !is.null(geneid)){
        my.vioplot(t(log2cpm_filtered[geneid,grp1]),t(log2cpm_filtered[geneid,grp2]),
                   names=c('Grp1','Grp2'),col=c("#e41a1c", "#377eb8"))
        title(main=input$violin_gene,ylab='log2CPM')
      }
      dev.off()
    })
  
})

# VIOPLOT_SOURCE ----------------------------------------------------------

my.vioplot<-function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
          horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
          lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
          at, add = FALSE, wex = 1, drawRect = TRUE) 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

library(shiny)
library(ggplot2)
library(Rtsne)
library(DT)
library(edgeR)

shinyServer(function(input, output) {
  
  set.seed(1)

  readqc<-read.csv('Data/QC_PASS.csv',row.names=1,stringsAsFactors = TRUE,as.is=T)
  counts<-read.csv('Data/Counts.csv',row.names=1,stringsAsFactors = TRUE, as.is=T, check.names=F)
  featuredata<-read.csv('Data/FeatureData.csv',row.names=1,stringsAsFactors = TRUE)
  
  libsize<-colSums(counts)
  cpm<-1e6*(sweep(counts,2,libsize,"/"))
  log2cpm<-log2(cpm+1)
  
  react_qc<-reactive({
    readqc_values<-readqc
    readqc_values$PASS.FAIL<-'FAIL'
    readqc_values[which(readqc_values$Genes_Detected > input$genes & 
                        readqc_values$Percent.Mitochondria < input$mitocontent),'PASS.FAIL']<-'PASS'
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
      print(p1)
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
      print(p1)
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
      print(p1)
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
      print(p1)
    }
    
  })
  
  
  output$click_info <- renderPrint({
    values=react_qc()
    nearPoints(values, input$qcplot_info, addDist = TRUE)
  })
  
  output$genes_plot <- renderPlot({
#    readqc_values<-readqc
#    readqc_values$PASS.FAIL<-'FAIL'
#    readqc_values[which(readqc_values$Genes_Detected > input$genes | 
#                          readqc_values$Percent.Mitochondria < input$mitocontent),'PASS.FAIL']<-'PASS'
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
    print(p1)
    #plot(readqc$Genes_Detected,readqc$Percent.Mitochondria,pch=16,col='red')
    })
  
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
    barplot(subset,names.arg=rownames(sorteddata)[1:10],cex.names=0.5)
  })
  
  output$tsne_plt<-renderPlot({
    tsne.data<-tsne.values()
    dim_x<-paste('V',input$dim1,sep='')
    dim_y<-paste('V',input$dim2,sep='')
    
    p1<-ggplot(tsne.data,aes_string(x=dim_x,y=dim_y))+
      geom_point(aes(color='red'))+
      theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
            axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
            strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16),
            legend.position="none")+
      xlab(paste('T-SNE',input$dim1,sep=''))+
      ylab(paste('T-SNE',input$dim2,sep=''))
    print(p1)
  })

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
  
  dge<-reactive({
    data<-filteredCounts()
    if (input$cluster == 'PCA'){
      pca.data<-pca.values()
      pca.frame<-as.data.frame(pca.data$x)
      
      grp1<-rownames(brushedPoints(pca.frame[,1:5],input$b1))
      grp2<-rownames(brushedPoints(pca.frame[,1:5],input$b2))
      
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
    if (dim(top.genes)[1]>1){
      DT::datatable(top.genes, options = list(orderClasses = TRUE,lengthMenu = c(10, 30, 50), pageLength = 10))
    }
  })
  
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
  
})

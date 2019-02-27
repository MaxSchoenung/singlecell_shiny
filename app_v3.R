library(ggplot2)
library(data.table)
library(shiny)
library(shinythemes)
library(shinybusy)
library(Seurat)

ui <- fluidPage(

  # Titel -------------------------------------------------------------------

  theme = shinythemes::shinytheme("cerulean"),
  titlePanel("Inflammation scRNA-seq"),
  br(),


  # Sidebar -----------------------------------------------------------------

  sidebarLayout(

    sidebarPanel(
      selectInput('dataset',label = "Choose the Dataset",choices = c("Mutants"="mut","Wildtype"="wt")),
      actionButton("loadData", "Load the Dataset"),
      textInput('gene', 'Gene ID Input', "Ccl3"),
      actionButton("go", "Show UMAP")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Cluster",
                 br(),
                 add_busy_spinner(spin = "fading-circle"),
                 column(12,plotOutput("plot.cluster"))),
        tabPanel("Gene Statistics",
                 br(),
                 selectInput('annotations.stat',label = "Annotation",choices = c("Genotype"=6,"Compartment"=5,"Cell Type"=10)),
                 br(),
                 column(12,plotOutput("bar.stats"))),
        tabPanel("Statistics per Genotype",
                 br(),
                 column(12,plotOutput("geno.stats"))),
        tabPanel("Subset Plots",
                 uiOutput("checkbox.subset"),
                 uiOutput("checkbox.subset2"),
                 br(),
                 column(12,plotOutput("subset.umap"),
                 plotOutput("subset.boxplot"))),
        tabPanel("Cluster Assignment",
                 br(),
                 selectInput('annotations',label = "Annotation",choices = c("Genotype"=6,"Compartment"=5,"Cell Type"=10,
                                                                            "Number of Genes"=2,
                                                                            "Number of UMI"=3)),
                 br(),
                 column(12,plotOutput("plot.stats")))
        
      )
    )
  ))


# Define the server code
options(shiny.maxRequestSize=30000*1024^2)
server <- function(input, output) {
  
  data <- eventReactive(input$loadData, {
    if(input$dataset=="mut"){
      readRDS("~/Documents/Data/10x_inflammation/AML06DFG_withoutBatchCorrection_AMA07022019.rds")
    }else if(input$dataset=="wt"){
      readRDS("~/Documents/Data/10x_inflammation/Wildtype_data_A06_AMA_06022019.rds")
    }
  })
  
  p <- eventReactive(input$go, {
    selected.gene <- which(rownames(data()@assays$RNA@scale.data)==input$gene)
    ggplot2::ggplot(data.frame(data()@reductions$umap@cell.embeddings),
                    aes(UMAP_1,UMAP_2,color=as.numeric(data()@assays$RNA@scale.data[selected.gene,])))+
      geom_point(alpha=0.6,size=0.5)+
      theme_classic()+
      scale_color_continuous(low="gray80",high="navy",
                             name=input$gene)
  })


  output$plot.cluster <- renderPlot({
    p()
  })
  
  output$plot.stats <- renderPlot({
    ggplot2::ggplot(data.frame(data()@reductions$umap@cell.embeddings),
                    aes(UMAP_1,UMAP_2,color=data()@meta.data[as.numeric(input$annotations)][,1]))+
      geom_point(alpha=0.3,size=0.5)+
      theme_classic()+
      theme(legend.title = element_blank())
  })
  
  output$bar.stats <- renderPlot({
    selected.gene <- which(rownames(data()@assays$RNA@scale.data)==input$gene)
    data.plot <- data.frame("ID"=factor(data()@meta.data[as.numeric(input$annotations.stat)][,1]),
                            "Expression"=as.numeric(data()@assays$RNA@scale.data[selected.gene,]))
    ggplot2::ggplot(data.plot,aes(ID,Expression,fill=ID))+
      theme_classic()+
      scale_fill_discrete(name="")+
      geom_violin()+
      geom_boxplot(aes(fill=NULL),width=.1,outlier.colour = NA)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  output$geno.stats <- renderPlot({
    selected.gene <- which(rownames(data()@assays$RNA@scale.data)==input$gene)
    data.plot <- data.frame("Genotype"=factor(data()@meta.data[as.numeric(6)][,1]),
                            "CellType"=factor(data()@meta.data[as.numeric(10)][,1]),
                            "Expression"=as.numeric(data()@assays$RNA@scale.data[selected.gene,]))
    ggplot2::ggplot(data.plot,aes(CellType,Expression,fill=Genotype))+
      theme_classic()+
      scale_fill_discrete(name="")+
      #geom_violin()+
      geom_boxplot(outlier.colour = NA)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  output$checkbox.subset <- renderUI({
    if(!is.null(data())){
      checkboxGroupInput("checkbox.subset.choose","Choose Cell Types",choices = levels(data()@meta.data$RNA_snn_res.3),inline=T)
      }
  })
  
  
  output$checkbox.subset2 <- renderUI({
    if(!is.null(data())){
      checkboxGroupInput("checkbox.subset.geno","Choose Genotype",choices = levels(factor(data()@meta.data$Genotype)),inline=T)
    }
  })
  
  output$subset.boxplot <- renderPlot({
    selected.gene <- which(rownames(data()@assays$RNA@scale.data)==input$gene)
    data.plot <- data.frame("Genotype"=factor(data()@meta.data[as.numeric(6)][,1]),
                            "CellType"=factor(data()@meta.data[as.numeric(10)][,1]),
                            "Expression"=as.numeric(data()@assays$RNA@scale.data[selected.gene,]))
    data.plot <- data.plot[which(as.character(data.plot$Genotype)%in%c(input$checkbox.subset.geno)&data.plot$CellType%in%c(input$checkbox.subset.choose)),]
    ggplot2::ggplot(data.plot,aes(CellType,Expression,fill=Genotype))+
      theme_classic()+
      scale_fill_discrete(name="")+
      #geom_violin()+
      geom_boxplot(outlier.colour = NA)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  output$subset.umap <- renderPlot({
    selected.gene <- which(rownames(data()@assays$RNA@scale.data)==input$gene)
    data.plot <- data.frame("Genotype"=factor(data()@meta.data[as.numeric(6)][,1]),
                            "CellType"=factor(data()@meta.data[as.numeric(10)][,1]),
                            "Expression"=as.numeric(data()@assays$RNA@scale.data[selected.gene,]))
    umap <- data.frame(data()@reductions$umap@cell.embeddings)
    umap <- umap[which(as.character(data.plot$Genotype)%in%c(input$checkbox.subset.geno)&data.plot$CellType%in%c(input$checkbox.subset.choose)),]
    data.plot <- data.plot[which(as.character(data.plot$Genotype)%in%c(input$checkbox.subset.geno)&data.plot$CellType%in%c(input$checkbox.subset.choose)),]
    
    ggplot2::ggplot(umap,aes(UMAP_1,UMAP_2,color=data.plot$Expression))+
      geom_point(alpha=0.6,size=0.5)+
      theme_classic()+
      scale_color_continuous(low="gray80",high="navy",
                             name=input$gene)
  })
  


}

# Run the application
shinyApp(ui = ui, server = server)





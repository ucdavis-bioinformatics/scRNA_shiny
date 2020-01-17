library(shiny)
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)


aggregate <- readRDS('Zhang_Seurat_UCD_blood.rds')
#aggregate
genes = aggregate@assays$RNA
meta_nums <- colnames(dplyr::select_if(aggregate@meta.data, is.numeric))
meta_cats <- colnames(dplyr::select_if(aggregate@meta.data, is.factor))
pcs <- list('PC_1','PC_2','PC_3','PC_4','PC_5','PC_6','PC_7','PC_8','PC_9')
agg_cats <- colnames(dplyr::select_if(aggregate@meta.data, is.factor))


server = function(input, output, session){
  outVar = reactive({
    if (input$dataset == 'Genes'){mydata=row.names(genes)}
    else if (input$dataset == 'Numeric Metadata') {mydata=meta_nums}
    else if (input$dataset == 'PCs') {mydata=pcs}
    mydata
  })
  
  outVar2 = reactive({
    if (input$dataset_single == 'Genes'){mydata=row.names(genes)}
    else if (input$dataset_single == 'Numeric Metadata') {mydata=meta_nums}
    else if (input$dataset_single == 'PCs') {mydata=pcs}
    mydata
  })
  
  observe({
    updateSelectInput(session, "numeric",
                      choices = outVar()
  )})
  
  observe({
    updateSelectInput(session, "numeric2",
                      choices = outVar()
  )})
  
  observe({
    updateSelectInput(session, "numeric_b",
                      choices = row.names(genes)
  )})
  
  observe({
    updateSelectInput(session, "numeric_single",
                      choices = outVar2()
    )})
  
  # formulaText <- reactive({
  #   paste("Marker Gene ~", input$numeric)
  # })

  # Marker Plot Double
  output$MarkerGenePlot <- renderPlot({
    FeaturePlot(
      aggregate,
      c(input$numeric, input$numeric2), blend=TRUE
    )
  })

  # Marker Plot Single
  output$MarkerGenePlotSingle <- renderPlot({
    FeaturePlot(
      aggregate,
      c(input$numeric_single)
    )
  })
  
  # Double Feature Categorical Feature Plot
  output$CategoricalPlot <- renderPlot({
    DimPlot(object = aggregate, group.by=input$categorical, pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
  })

  # Single Feature Categorical Feature Plot
  output$CategoricalPlotSingle <- renderPlot({
    DimPlot(object = aggregate, group.by=input$categorical_single, pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
  })
  
  # Double Feature Violin Plot
  output$ViolinPlot <- renderPlot({
    Idents(aggregate) <- input$categorical
    VlnPlot(object =  aggregate, features = c(input$numeric, input$numeric2), pt.size = 0.05)
  })

  # Single Feature Violin Plot
  output$ViolinPlotSingle <- renderPlot({
    Idents(aggregate) <- input$categorical_single
    VlnPlot(object =  aggregate, features = c(input$numeric_single), pt.size = 0.05)
  })
  
  # Marker Set Plot
  output$MarkerSet <- renderPlot({
    Idents(aggregate) <- input$categorical_b
    print(input$numeric_b)
    markers = input$numeric_b
    expr.cutoff = 3
    widedat <- FetchData(aggregate, markers)
    #print(widedat)
    widedat$Cluster <- Idents(aggregate)
    longdat <- gather(widedat, key = "Gene", value = "Expression", -Cluster)
    longdat$Is.Expressed <- ifelse(longdat$Expression > expr.cutoff, 1, 0)
    longdat$Cluster <- factor(longdat$Cluster)
    longdat$Gene <- factor(longdat$Gene, levels = markers)

    # Need to summarize into average expression, pct expressed (which is also an average)
    plotdat <- group_by(longdat, Gene, Cluster) %>% summarize(`Percentage of Expressed Cells` = mean(Is.Expressed), `Mean Expression` = mean(Expression))
    ggplot(plotdat, aes(x = Gene, y = Cluster)) +
      geom_point(aes(size = `Percentage of Expressed Cells`, col = `Mean Expression`)) +
      labs(size = "Percentage\nof Expressed\nCells", col = "Mean\nExpression", x = NULL) +
      scale_color_gradient(low = "grey", high = "slateblue4") + theme_grey(base_size = 15) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # }, height = 1000, width = 900 )
  }, height = 1000)

}


ui <- fluidPage(
  
  titlePanel("scRNA Seurat Analysis"), 
  sidebarLayout(
    sidebarPanel(width = 12,
                 tabsetPanel(
                   tabPanel("Documentation", value=-999,
                            mainPanel(width = 12,
                              br(),
                              #includeMarkdown("docs/testing.md")
                              includeMarkdown("README.md")
                            )
                   ),
                   
                   tabPanel("Double Marker", value=2,
                                br(),
                                div(style="display: inline-block;vertical-align:top; width: 24%;",
                                    selectInput("dataset", "Numeric Analysis Type:",
                                            c('Genes', 'Numeric Metadata','PCs'))),
                                div(style="display: inline-block;vertical-align:top; width: 24%;",
                                    selectInput("categorical", "Identity:",
                                            c(meta_cats))),
                                div(style="display: inline-block;vertical-align:top; width: 24%;",
                                    selectInput("numeric", "Primary Numeric:", "")),
  
                                div(style="display: inline-block;vertical-align:top; width: 24%;",
                                    selectInput('numeric2', 'Secondary Numeric', "")),
                            
                              mainPanel(width = 12,
                                 br(),
                                 br(),
                                 #h3(textOutput("caption")),
                                 plotOutput("MarkerGenePlot"),
                                 plotOutput("ViolinPlot"),
                                 plotOutput("CategoricalPlot")
                              )
                          ),
                   tabPanel("Single Marker", value=3,
                            br(),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("dataset_single", "Numeric Analysis Type:",
                                        c('Genes', 'Numeric Metadata','PCs'))),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("categorical_single", "Identity:",
                                            c(meta_cats))),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("numeric_single", "Primary Numeric:", "")),

                            mainPanel(width = 12,
                                      br(),
                                      br(),
                                      #h3(textOutput("caption")),
                                      plotOutput("MarkerGenePlotSingle"),
                                      plotOutput("ViolinPlotSingle"),
                                      plotOutput("CategoricalPlotSingle")
                            )
                   ),
                  tabPanel("Marker Set (Grid)", value=4,
                            br(),
                            selectInput("categorical_b", "Identity:",
                                        c(agg_cats)),
                            selectInput("numeric_b", "Primary Numeric:", "", multiple=TRUE),
                             mainPanel(width = 12,
                                br(),
                                br(),
                                plotOutput("MarkerSet")
                             )
                          ),
                  id = "tabselected"
                 )
    ),
    mainPanel(width = 12)
  )
)


shinyApp(ui, server)



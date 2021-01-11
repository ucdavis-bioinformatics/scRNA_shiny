library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(markdown)
library(tidyr)

# Some initial setup:
# take in the file, get list of genes, get metadata numbers and categories, get pcs 1-9, and factors..
aggregate <- readRDS('../../shiny.rds')
genes = aggregate@assays$RNA
reductions <- attributes(aggregate@reductions)
meta_nums <- colnames(dplyr::select_if(aggregate@meta.data, is.numeric))
meta_cats <- c(colnames(dplyr::select_if(aggregate@meta.data, is.character)), colnames(dplyr::select_if(aggregate@meta.data, is.factor)))
pcs <- list('PC_1','PC_2','PC_3','PC_4','PC_5','PC_6','PC_7','PC_8','PC_9')
use.pcs <- 1:29
agg_cats <- colnames(dplyr::select_if(aggregate@meta.data, is.factor))
# TODO get reduction types as a list to choose from

# Main function of the program
server = function(input, output, session){
  
  # update values based on input from ui
  outVar_double = reactive({
    if (input$dataset == 'Genes'){mydata=row.names(genes)}
    else if (input$dataset == 'Numeric Metadata') {mydata=meta_nums}
    else if (input$dataset == 'PCs') {mydata=pcs}
    mydata
  })
  
  # update values based on input from ui
  outVar_single = reactive({
    if (input$dataset_single == 'Genes'){mydata=row.names(genes)}
    else if (input$dataset_single == 'Numeric Metadata') {mydata=meta_nums}
    else if (input$dataset_single == 'PCs') {mydata=pcs}
    mydata
  })
  
  # update values based on input from ui
  outVar_seperated = reactive({
    if (input$dataset_seperated == 'Genes'){mydata=row.names(genes)}
    else if (input$dataset_seperated == 'Numeric Metadata') {mydata=meta_nums}
    else if (input$dataset_seperated == 'PCs') {mydata=pcs}
    mydata
  })

  getResChoices = reactive({
    mydata = levels(eval(call("$", aggregate, input$identity_table)))
    mydata
  })
    
  # Reduction Type for the Single Marker Plot
  observe({
    updateSelectInput(session, "reduction_single",
                      choices = reductions
    )})
  
  # Reduction Type for the Double Marker Plot
  observe({
    updateSelectInput(session, "reduction_double",
                      choices = reductions
    )})
  
  # Primary numeric value in the double marker plot
  observe({
    updateSelectInput(session, "numeric",
                      choices = outVar_double()
    )})
  
  # Secondary numeric value in the double marker plot
  observe({
    updateSelectInput(session, "numeric2",
                      choices = outVar_double()
    )})
  
  # Numeric input list for the marker set (multiple =TRUE)
  observe({
    updateSelectInput(session, "numeric_b",
                      choices = row.names(genes)
    )})
  
  # Only numeric input for the single marker plot
  observe({
    updateSelectInput(session, "numeric_single",
                      choices = outVar_single()
    )})
  
  # Cluster Tree identity
  observe({
    updateSelectInput(session, "identity_tree",
                      choices = meta_cats
    )})
  
  # Seperated Identity
  observe({
    updateSelectInput(session, "identity_seperated",
                      choices = meta_cats
    )})
  
  # Seperated Numeric
  observe({
    updateSelectInput(session, "numeric_seperated",
                      choices = outVar_seperated()
    )})
  
  # Seperated Reduction
  observe({
    updateSelectInput(session, "reduction_seperated",
                      choices = reductions
    )})
  
  
  # Table Identity
  observe({
    updateSelectInput(session, "identity_table",
                      choices = meta_cats
    
    )})
  
  
  # Table Marker
  observe({
    updateSelectInput(session, "compare_table",
                      choices =getResChoices()
    )})
  
  # Table Compare
  observe({
    updateSelectInput(session, "markers_table",
                      choices = getResChoices()
    )})
  

  
  
  # Marker Plot Double
  output$MarkerGenePlot <- renderPlot({
    FeaturePlot(
      aggregate,
      c(input$numeric, input$numeric2),
      blend=TRUE,
      reduction=input$reduction_double
    )
  })
  
  # Marker Plot Single
  output$MarkerGenePlotSingle <- renderPlot({
    FeaturePlot(
      aggregate,
      c(input$numeric_single),
      reduction=input$reduction_single
    )
  })
  
  # Double Feature Categorical Feature Plot
  output$CategoricalPlot <- renderPlot({
    Idents(aggregate) <- input$categorical
    order <- sort(levels(aggregate))
    levels(aggregate) <- order
    DimPlot(object = aggregate, pt.size=0.5, reduction = input$reduction_double, label = T)
  })
  
  # Single Feature Categorical Feature Plot
  output$CategoricalPlotSingle <- renderPlot({
    Idents(aggregate) <- input$categorical_single
    order <- sort(levels(aggregate))
    levels(aggregate) <- order
    DimPlot(object = aggregate, group.by=input$categorical_single, pt.size=0.5, reduction = input$reduction_single, label = T)
  })
  
  # Double Feature Violin Plot
  output$ViolinPlot <- renderPlot({
    Idents(aggregate) <- input$categorical
    order <- sort(levels(aggregate))
    levels(aggregate) <- order
    VlnPlot(object =  aggregate, features = c(input$numeric, input$numeric2), pt.size = 0.05)
  })
  
  # Single Feature Violin Plot
  output$ViolinPlotSingle <- renderPlot({
    Idents(aggregate) <- input$categorical_single
    order <- sort(levels(aggregate))
    levels(aggregate) <- order
    VlnPlot(object =  aggregate, features = c(input$numeric_single), pt.size = 0.05)
  })
  
  # Cluster Tree Plot
  output$ClusterTree <- renderPlot({
    Idents(aggregate) <- input$identity_tree
    aggregate <- BuildClusterTree(
      aggregate, dims = use.pcs)
    PlotClusterTree(aggregate)
  })
  
  
  
  
  # Seperated Feature Plot
  output$SeperatedFeature <- renderPlot({
    Idents(aggregate) <- input$identity_seperated
    order <- sort(levels(aggregate))
    levels(aggregate) <- order
    FeaturePlot(aggregate, c(input$numeric_seperated), reduction=input$reduction_seperated,
      split.by = "orig.ident"
    )
  })
  
  # Seperated Dim Plot
  output$SeperatedDim <- renderPlot({
    Idents(aggregate) <- input$identity_seperated
    order <- sort(levels(aggregate))
    levels(aggregate) <- order
    DimPlot(aggregate, reduction=input$reduction_seperated,
                split.by = "orig.ident"
    )
  })
  
  
  
  # Seperated Counts table
  output$SeperatedCounts <- renderTable({
    lab_list = c()
    identities = c()
    num_list = c()
    
    marker = c(input$numeric_seperated)
    Idents(aggregate) <- input$identity_seperated
    
    if(input$dataset_seperated == 'Numeric Metadata'){
      nm <- data.frame(matrix(unlist(eval(call('$', aggregate, marker[1]))), nrow=length(eval(call('$', aggregate, marker[1]))), byrow=T))
      colnames(nm) = marker
      rownames(nm) = labels(eval(call('$', aggregate, marker[1])))
      widedat <- nm
    }
    else{widedat <- FetchData(aggregate, marker)}
    
    widedat$Cluster <- Idents(aggregate)
    widedat$orig.ident=sapply(strsplit(rownames(widedat), "-"), "[", 2)
    widedat$final <- paste(widedat$orig.ident, widedat$Cluster, sep="_")
    final_object = (aggregate(widedat[, 1:2], list(widedat$final), mean)[1:2])
    rownames(final_object) = final_object$Group.1
    tab =t(final_object[-c(1)])
    #tab = final_object#suppressMessages(AverageExpression(aggregate, add.ident = input$identity_seperated))
    
    for (i in names(tab[input$numeric_seperated,])){lab_list = c(lab_list,strsplit(i,'_')[[1]][1])}
    for (i in names(tab[input$numeric_seperated,])){identities = c(identities,strsplit(i,'_')[[1]][2])}
    for (i in tab[input$numeric_seperated,]){num_list = c(num_list,i)}
    df = as.data.frame(pivot_wider(data.frame(identities, num_list, lab_list), names_from = 'lab_list', values_from = 'num_list'))
    df[is.na(df)] <- 0
    
    df_p = as.data.frame.matrix(prop.table((table(eval(call("$", aggregate, input$identity_seperated)), eval(call("$", aggregate, "orig.ident")))),2))
    #df_p$identities = rownames(df_p)
    #df_p=df_p[, !(names(df_p)) %in% c('identities')]
    df_p=df_p/colSums(df_p)

    
    merged_final = as.data.frame.matrix(merge(df, df_p, by.x = 'identities', by.y = 'row.names', suffixes = c(".AvgExpression",".Proportion")))
    merged_final
    #as.data.frame(c(1,2))
    #as.data.frame.matrix(table(eval(call("$", aggregate, input$identity_seperated)), eval(call("$", aggregate, "orig.ident"))))
  }, width = "100%", colnames=TRUE, rownames=TRUE, digits=4)
  
  
  # Marker Table
  output$markers <- renderTable({
    Idents(aggregate) <- input$identity_table
    if (as.logical(length(c(input$compare_table)))){FindMarkers(aggregate, ident.1=input$markers_table, ident.2=input$compare_table)}
    else {FindMarkers(aggregate, ident.1=input$markers_table)}
  }, rownames = TRUE, colnames = TRUE, width = "100%", digits=-5)
  
  # Marker Set Plot
  output$MarkerSet <- renderPlot({
    Idents(aggregate) <- input$categorical_b
    markers = input$numeric_b
    expr.cutoff = 3
    widedat <- FetchData(aggregate, markers)
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
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("dataset", "Numeric Analysis Type:",
                                            c('Genes', 'Numeric Metadata','PCs'))),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("reduction_double", "Reduction:",
                                            c(reductions))),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("categorical", "Identity:",
                                            c(meta_cats))),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("numeric", "Primary Numeric:", "")),
                            
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
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
                                selectInput("reduction_single", "Reduction:",
                                            c(reductions))),
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
                            selectizeInput("numeric_b", "Primary Numeric (csv format works here if pasted in):", "", 
                                           options = list(
                                             delimiter = ',',
                                             create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                           selected = NULL, multiple = TRUE), ## and switch multiple to True,
                            mainPanel(width = 12,
                                      br(),
                                      br(),
                                      plotOutput("MarkerSet")
                            )
                   ),
                   
                   tabPanel("Cluster Tree", value=5,
                            br(),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("identity_tree", "Identity:",
                                            c(meta_cats))),
                            mainPanel(width = 12,
                                      br(),
                                      br(),
                                      #h3(textOutput("caption")),
                                      plotOutput("ClusterTree"),
                            )
                   ),
                   tabPanel("Seperated Feature", value=6,
                            br(),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("dataset_seperated", "Numeric Analysis Type:",
                                            c('Genes', 'Numeric Metadata','PCs'))),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("reduction_seperated", "Reduction:",
                                            c(reductions))),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("identity_seperated", "Identity:",
                                            c(meta_cats))),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("numeric_seperated", "Primary Numeric:", "")),

                            mainPanel(width = 12,
                                      br(),
                                      br(),
                                      #h3(textOutput("caption")),
                                      plotOutput("SeperatedFeature"),
                                      plotOutput("SeperatedDim"),
                                      tableOutput("SeperatedCounts")
                                      
                            )
                   ),
                   tabPanel("Marker Table", value=7,
                            br(),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("identity_table", "Identity:",
                                            c(meta_cats))),
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("markers_table", "Get markers for:", "", multiple = TRUE)),
                            
                            div(style="display: inline-block;vertical-align:top; width: 24%;",
                                selectInput("compare_table", "Compare to (blank is all other groups):", "", multiple = TRUE)),
                            
                            mainPanel(width = 12,
                                      br(),
                                      br(),
                                      #h3(textOutput("caption")),
                                      tableOutput("markers")
                            )
                   ),
                   id = "tabselected"
                 )
    ),
    mainPanel(width = 12)
  )
)

# Potential to do, add DimPlot or HeatMap


shinyApp(ui, server)

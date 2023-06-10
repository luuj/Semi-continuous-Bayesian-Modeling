library(shiny)
library(ggplot2 )
source("../gendata.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Clustered data generation"),
    
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(6,numericInput("bSigma", "Subject sigma", value=1, width=100, min=0, max=15, step=0.5)),
          column(6,numericInput("vSigma", "Cluster sigma", value=1.5, width=100, min=0, max=15, step=0.5))
        ),
        
        # Spacing
        br(),
        
        fluidRow(
          column(6,numericInput("beta1", "beta1", value=1.3, width=100, step=0.05)),
          column(6,numericInput("alpha1", "alpha1", value=-3, width=100, step=0.05))
        ),
        
        fluidRow(
          column(6,numericInput("beta2", "beta2", value=-0.5, width=100, step=0.05)),
          column(6,numericInput("alpha2", "alpha2", value=-3.5, width=100, step=0.05))
        ),
        
        fluidRow(
          column(6,numericInput("betat", "betat", value=1, width=100, step=0.05)),
          column(6,numericInput("alphat", "alphat", value=0.5, width=100, step=0.05))
        ),
        
        fluidRow(
          column(6,numericInput("betam", "betam", value=3, width=100, step=0.05)),
          column(6,numericInput("alpham", "alpham", value=1, width=100, step=0.05))
        ),
        
        fluidRow(
          column(6,numericInput("beta1y", "beta1y", value=0.8, width=100, step=0.05)),
          column(6,numericInput("beta2y", "beta2y", value=-1.3, width=100, step=0.05))
        ),
        
        submitButton("Update View", icon("refresh"))
      ),
      
      mainPanel(
        tabsetPanel(type = "tabs",
          tabPanel("Spaghetti Plot", plotOutput("distPlot")),
          tabPanel("Histograms", plotOutput("histDth"), plotOutput("histCost")),
          tabPanel("Scatter Plot", plotOutput("scatPlot"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    # Update dataset when parameters are changed
    dat <- reactive({
      set.seed(10)
      # gen_data(n_k=70, n_t=12, k=250,
      #   bSigma=(diag(2)*input$bSigma),
      #   vSigma=(diag(3)*input$vSigma),
      #   coefs=c(input$beta1, input$beta2, input$betat, input$beta1y, input$beta2y, input$betam,
      #           input$alpha1, input$alpha2, input$alphat, input$alpham))
      gen_data(n_k=70, n_t=12, k=250)
    })
    
    output$distPlot <- renderPlot({
      plot.c1 <- dat() %>% filter(clst==1) %>% ggplot(aes(t,cumcost,group=id,color=factor(trt))) + geom_line() +
          xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost over time C1") + theme(legend.position = "none")
      plot.c2 <- dat() %>% filter(clst==2) %>% ggplot(aes(t,cumcost,group=id,color=factor(trt))) + geom_line() +
          xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost over time C2") + theme(legend.position = "none")
      plot.c3 <- dat() %>% filter(clst==3) %>% ggplot(aes(t,cumcost,group=id,color=factor(trt))) + geom_line() +
          xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost over time C3") + theme(legend.position = "none")
      plot.c4 <- dat() %>% filter(clst==4) %>% ggplot(aes(t,cumcost,group=id,color=factor(trt))) + geom_line() +
          xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost over time C4") + theme(legend.position = "none")
      plot.c5 <- dat() %>% filter(clst==5) %>% ggplot(aes(t,cumcost,group=id,color=factor(trt))) + geom_line() +
        xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost over time C5") + theme(legend.position = "none")
      plot.c6 <- dat() %>% filter(clst==6) %>% ggplot(aes(t,cumcost,group=id,color=factor(trt))) + geom_line() +
        xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost over time C6") + theme(legend.position = "none")
      
      gridExtra::grid.arrange(plot.c1, plot.c2, plot.c3, plot.c4, plot.c5, plot.c6, ncol=2)
    })
      
    output$histDth <- renderPlot({
      dthPlot <- clusterSummary(t, mean, dat())
      dthPlot$trt<- factor(dthPlot$trt, labels=c("Trt 0", "Trt 1"))
      dthPlot %>% ggplot(aes(max)) + geom_histogram(color="black", fill="white", bins = 11) +
        ylab("Number of clusters") + ggtitle("Average time to death distribution by cluster") + 
        facet_wrap(~trt) + xlab("Time to death (months)")
    })
    
    output$histCost <- renderPlot({
      costPlot <- clusterSummary(cumcost, mean, dat())
      costPlot$trt<- factor(costPlot$trt, labels=c("Trt 0", "Trt 1"))
      costPlot %>% ggplot(aes(max)) + geom_histogram(color="black", fill="white", bins = 100) +
        ylab("Number of clusters") + ggtitle("Average total cost accrued distribution by cluster") +
        facet_wrap(~trt,  scales = "free") + xlab("Total cost accrued (dollars)")
    })
    
    output$scatPlot <- renderPlot({
      dthprop <- clusterSummary(y2, mean, dat())
      costprop <- clusterSummary(y1, mean, dat())
      scatplot <- left_join(dthprop, costprop, by=c("clst","trt"))
      colnames(scatplot) <- c("Cluster", "Trt", "Death", "Cost")
      scatplot$Trt <- factor(scatplot$Trt)
      
      # Scatter plot for proportion of at least one cost and death for each cluster
      scatplot %>% ggplot(aes(Death, Cost, color=Trt)) + geom_point() + ylab("Proportion with any cost accrued") +
        xlab("Proportion of deaths") + ggtitle("Cost vs. death proportions by cluster")
    })
}

# Return dataset with maximum value of a covar for each id by cluster
clusterSummary <- function(covar, fn, dat, includeTrt = T){
  pars <- as.list(match.call()[-1]) # Convert covar into variable
  out <- aggregate(eval(pars$covar) ~ id+clst+trt, dat, max)
  colnames(out) <- c("id", "clst", "trt", "max")
  if (includeTrt)
    return(aggregate(max~clst+trt, out, fn))
  else
    return(aggregate(max~clst, out, fn))
}

# Run the application 
shinyApp(ui = ui, server = server)


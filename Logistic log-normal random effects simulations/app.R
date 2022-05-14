library(shiny)
source("C:\\Users\\LuuJonathan\\Dropbox\\Jonathan-Sebastian\\Dissertation\\Code\\Logistic log-normal random effects simulations\\TPM_Simulations.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Two-part model data generation"),
    
    fluidRow(
        column(1,numericInput("n", "# of samples", value=1000,width=100, min=1, max=30000)),
        column(1,numericInput("ymax", "Ymax", value=3000,width=100, min=1, max=30000))
    ),
    
    fluidRow(
        column(1,numericInput("re_sigma1", "RE_Sigma1", value=1,width=100,min=0.01, step=0.05)),
        column(1,numericInput("re_sigma2", "RE_Sigma2", value=2,width=100,min=0.01, step=0.05)),
        column(1,numericInput("re_rho", "RE_rho", value=0.2,width=100,min=-1, max=1,step=0.05))
    ),
    
    fluidRow(
        column(1,numericInput("ln_sigma", "LN_sigma", value=1,width=100,min=0.01, step=0.05))
    ),
    
    fluidRow(
        column(1,numericInput("beta1", "Beta1", value=0,width=100, step=0.25)),
        column(1,numericInput("beta2", "Beta2", value=0,width=100, step=0.25)),
        column(1,numericInput("beta3", "Beta3", value=1,width=100, step=0.25)),
        column(1,numericInput("gamma1", "Gamma1", value=0,width=100, step=0.25)),
        column(1,numericInput("gamma2", "Gamma2", value=0,width=100, step=0.25)),
        column(1,numericInput("gamma3", "Gamma3", value=1.5,width=100, step=0.25))
    ),

    mainPanel(
       plotOutput("distPlot"),
    )
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$distPlot <- renderPlot({
        dat <- genData(beta = c(input$beta1, input$beta2, input$beta3), 
                       gamma=c(input$gamma1, input$gamma2, input$gamma3),
                       re_sigma1 = input$re_sigma1, re_sigma2 = input$re_sigma2, 
                       re_rho = input$re_rho,
                       ln_sigma = input$ln_sigma)
        
        sub_dat <- dat %>% filter(id %in% sample(dat$id, input$n))
        cumu_p <- ggplot(data=sub_dat, aes(x=time, y=csum, group=id))
        cumu_p + geom_line() + xlab("Days") + ylab("Cumulative sum") + 
            ggtitle("Cumulative cost sum over 3 days", 
                    subtitle = paste0("Alpha={", input$beta1, ",", input$beta2, ",", input$beta3, "}", " Beta={", input$gamma1, ",", input$gamma2, ",", input$gamma3, "}")) + ylim(0,input$ymax)

    })
}

# Run the application 
shinyApp(ui = ui, server = server)



# Visualize the effect of gender - two plots side by side
# y axis on the side to see marginal probability of zero
# look into analysis procedure for tpm - see Neelon
    




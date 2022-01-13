library(shiny)
library(tidyverse)
library(lavaan)
library(shinycssloaders)
library(here)
source(here("Simulator_PokropekEtAl.R"))
source(here("Detector_ByrneVandeVijer.R"))
source(here("Detector_MInd.R"))
source(here("Detector_Rieger_v3.R"))
source(here("runsim.R"))


# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Performance of Detection Methods for Non-invariant Items in CFA."),
  strong("Pit Rieger - Jan 2022"),
  p("This simple shinyapp is an interactive version of the simulation study in my master thesis. It allows the user to specify several
  parameters for generating data from a single-factor CFA model with item-level violations of measurement invariance. Through replications, the app provides
    estimates for the sensitivity and specificity of four different detection methods. Further details can be found on the", 
    a("GitHub Repo", href = "https://github.com/pitrieger/masterthesis"), 
    "and in the thesis itself which is available on",
    a("Overleaf", href = "https://www.overleaf.com/read/ggznrkmtkwxp"), 
    ". Note that this app is rather slow, so please be patient when hitting 'Simulate'."),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      sliderInput(inputId = "n",
                  label = "Number of observations per group",
                  min = 50,
                  max = 500,
                  step = 10,
                  value = 250),
      sliderInput(inputId = "g",
                  label = "Number of groups",
                  min = 2,
                  max = 48,
                  step = 2,
                  value = 2),
      sliderInput(inputId = "p",
                  label = "Number of items",
                  min = 3,
                  max = 10,
                  value = 5),
      sliderInput(inputId = "k",
                  label = "Number of non-invariant items",
                  min = 1,
                  max = 10,
                  value = 2),
      sliderInput(inputId = "loadingbias",
                  label = "Bias on loadings",
                  min = 0,
                  max = 2,
                  step = 0.05,
                  value = 0.2),
      sliderInput(inputId = "interceptbias",
                  label = "Bias on intercepts",
                  min = 0,
                  max = 2,
                  step = 0.05,
                  value = 0.2),
      sliderInput(inputId = "nsim",
                  label = "Number of simulations (may take substantial amount of time for large values!)",
                  min = 3,
                  max = 100,
                  value = 20),
      actionButton("go", "Simulate")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      strong("Dot Plot"),
      withSpinner(
        plotOutput(outputId = "dotPlot")
      ),
      strong("Raw Confusion Matrix Entries"),
      withSpinner(
        tableOutput(outputId = "confusionmat")
      )
    )
  )
)

server <- function(input, output) {
  pars = reactiveValues(n = NULL,
                        g = NULL,
                        p = NULL, 
                        k = NULL, 
                        loadingbias = NULL, 
                        interceptbias = NULL, 
                        nsim = NULL)
  confusion = reactiveValues(sim_out_df = NULL)

  observeEvent(input$go,{
    pars$n = input$n
    pars$g = input$g
    pars$p = input$p
    pars$k = input$k
    pars$loadingbias = input$loadingbias
    pars$interceptbias = input$interceptbias
    pars$nsim = input$nsim
    if (pars$k > pars$p) stop("k cannot be larger than p")
    
    sim_out = replicate(pars$nsim, run_sim(pars))
    confusion$sim_out_sum = apply(sim_out, 1:2, function(x) sum(x, na.rm = T)) 
  })

  #output$dotPlot <- renderPlot({
  #  if (is.null(pars$n)) return()
  #  hist(rnorm(pars$nsim))
  #})  
  
  output$dotPlot <- renderPlot({
    if (is.null(pars$n)) return()
    plot_df = rbind(cbind(get_sensspec(confusion$sim_out_sum), type = "Sensitivity"), 
                    cbind(get_sensspec(confusion$sim_out_sum, type = "specificity"), type = "Specificity"))
    plot_df$method = rep(c("MInd", "BV", "R1", "R2"), 2)
    plot_df$method = factor(plot_df$method, levels = rev(c("MInd", "BV", "R1", "R2")))
    
    ggplot(plot_df, aes(xmin = lowerCI, xmax = upperCI, x = est, y = method)) +
      geom_pointrange() + 
      facet_wrap(~type, nrow = 2) + 
      scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
      theme_bw() + 
      labs(x = "Estimate") + 
      theme(axis.title.y = element_blank(),
            strip.background = element_rect(color = "black", fill = "white"),
            strip.text = element_text(face = "bold"))
  })
  
  output$confusionmat <- renderTable({
    if (is.null(pars$n)) return()
    tab = cbind(Method = c("MInd", "BV", "R1", "R2"), confusion$sim_out_sum)
    tab
  })
  
}

shinyApp(ui = ui, server = server)

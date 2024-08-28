################################################################################
# An R shiny to illustrate the shapes of the harmonic oscillator hazard function
################################################################################

# Required packages
library(shiny)

# Routines

source("routinesHO.R")



# minimum time for plot
#tmin = -100
# maximum time for plot
tmax = 100
# Minimum value of parameter 1
p1min = 0
# Minimum value of parameter 2
p2min = 0
# Minimum value of parameter 3
p3min = 0
# Minimum value of parameter 4
p4min = 0
# Minimum value of parameter 5
p5min = -10

# Maximum value of parameter 1
p1max = 10
# Maximum value of parameter 2
p2max = 10
# Maximum value of parameter 3
p3max = 10
# Maximum value of parameter 4
p4max = 10
# Maximum value of parameter 5
p5max = 10


uiHO <- fluidPage(
  # Title
  titlePanel("Harmonic Oscillator Hazard"),
  
  # Sidebar for parameter selection
  sidebarLayout(
    sidebarPanel(
      sliderInput("t", "t", min = 0, max = tmax, value = 10),
      sliderInput("par1", "Parameter 1 (eta)", min = p1min, max = p1max, value = 1, step = 0.1),
      sliderInput("par2", "Parameter 2 (w0)", min = p2min, max = p2max, value = 1, step = 0.1),
      sliderInput("par3", "Parameter 3 (hb)", min = p3min, max = p3max, value = 1, step = 0.1),
      sliderInput("par4", "Parameter 4 (h0)", min = p4min, max = p4max, value = 1, step = 0.1),
      sliderInput("par5", "Parameter 5 (r0)", min = p5min, max = p5max, value = 0, step = 0.1)
    ),
    
    # Main panel for plot
    mainPanel(
      plotOutput("functionPlot")
    )
  )
)

serverHO <- function(input, output) {
  output$functionPlot <- renderPlot({
    tt = seq(from = 0, to = input$t, length = 1000)
    # Calculate function values for current parameters
    y1 <- hHO(tt, input$par1, input$par2, input$par3, input$par4, input$par5)
    
     # Plot the function
    plot(tt, y1, type = "l", 
         main = paste("h(t, eta = ", input$par1, 
                      ", w0 = ", input$par2, 
                      ", hb = ", input$par3,
                      ", h0 = ", input$par4,
                      ", r0 = ", input$par5, "),", sep=""),
         xlab = "t", ylab = "h(t)",
         cex.lab = 1.5, cex.axis = 1.5, lwd = 2)
    abline(h = input$par3, col = "blue", lty = 2, lwd = 2)
    abline(h = 0, col = "red", lty = 1, lwd = 2)
  })
}

# Run the app
shinyApp(ui = uiHO, server = serverHO)



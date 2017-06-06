library(shiny)
library(plotly)
# devtools::install_github("ropensci/plotly#1035", force = T)
ui <- fluidPage(
  plotlyOutput("plot")
)

server <- function(input, output) {
  output$plot <- renderPlotly({
    plot_ly(x = 1:10, y = 10:1, z = 1:10)
  })
}

shinyApp(ui, server)

# sessionInfo()

# require(devtools)
# install_version("plotly", version = "4.5.6")
# packageVersion('shiny')
# install_version("shiny", version = "1.0.0")
# packageVersion('plotly')
# packageVersion('ggplot2')
# install_version("ggplot2", version = "2.2.1")
#
#
# devtools::install_github("ropensci/plotly#1035")

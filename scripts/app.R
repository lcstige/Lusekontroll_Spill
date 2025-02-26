library(shiny)
library(ggplot2)
## Load initials
source('global.R')
# Update model ------------------------------------------------------------
source('update.R')
summarised_data <- summarise_data()
ui <- fluidPage(
    
    # Application title
    titlePanel("Lusespill"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
         sidebarPanel(
        #     sliderInput("TreatmentThreshold",
        #                 "Treatment threshold:",
        #                 min = 0.1,
        #                 max = 1.0,
        #                 value = 0.5),
        #     actionButton("Go","Test")
        # ),
         radioButtons("Treatmenttype", "Type behandling:",
                          c("Hydrogenperoksid" = "HPcht", 
                            "Deltametrin" = "DMcht", 
                            "Azametifos" = "AZcht", 
                            "Emamectin" = "EMcht", 
                            "Diflubenzuron" = "DBcht",
                            "Termisk" = "therm", 
                            "Ferskvann" = "freshw", 
                            "Mekanisk" = "mech")
         ),
         
        radioButtons("TreatmentThreshold", "Behandle?:",
                     c("Ja" = 1,
                       "Nei" = 0)
                     ),
                     
           actionButton("Go","Test"),
        
        radioButtons("Ncages", "Antall merder:",
                     c("1" = 1,
                       "2" = 2,
                       "3" = 3,
                       "4" = 4,
                       "5" = 5,
                       "6" = 6)
                     ),
        ),
        
        mainPanel(
            plotOutput("sim_plot"),
            tableOutput('table'), 
            textOutput('text')
            
        )
    )
)


server <- function(input, output){
  ymax.plot <- 20
  logoffset <- .01
  labels <- c(0,.1,.5,2,5,20,50,100,1000,1e4,1e5,1e6,1e7,1e8,1e9,1e19)
  yat    <- log10(labels+logoffset)
  #Ncages <- input$Ncages

    output$sim_plot <- renderPlot({
      ggplot(summarised_data, aes(day, Y.CH)) +
        scale_y_continuous(name="Lus pr laks", breaks = yat, labels = labels, limits=c(min(yat), log10(ymax.plot+logoffset))) +
        geom_point(shape = 1, colour = "blue") +
        geom_point(aes(day, Y.AF),shape = 1, colour = "red") +
        facet_wrap(~cage,  ncol=1) 
    })
    output$table <- renderTable(summarised_data %>% head(10))
    observeEvent(input$Go,{
        SV_T <- update_SV( do_treat = input$TreatmentThreshold, trt.type = input$Treatmenttype)
        t <<- SV_T$t
        SV <<- SV_T$SV
        summarised_data <<- summarise_data()

        output$sim_plot <- renderPlot({
            ggplot(summarised_data, aes(day, Y.OM)) +
                scale_y_continuous(name="Lus pr laks", breaks = yat, labels = labels, limits=c(min(yat), log10(ymax.plot+logoffset))) +
                geom_point(shape = 1, colour = "blue") +
                geom_point(aes(day, Y.AF),shape = 1, colour = "red") +
                geom_vline(xintercept = summarised_data$treatment +1, colour = "rosybrown") +
                facet_wrap(~cage,  ncol=1) 
          
        })
         output$table <- renderTable(summarised_data %>% na.omit %>% tail)
         output$tekst <- renderText(input$Treatmenttype)
        
    })



}

shinyApp(ui=ui, server=server)
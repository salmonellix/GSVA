dashboardPage(
  title = "GSVA Shiny Application",
  
  dashboardHeader(
    tags$li(class = "dropdown",
            tags$div(id = "app_title", "GSVA Shiny Application")
    ),
    title = tags$img(src="GSVA.png", height=75, width=75)
  ),
  
  dashboardSidebar(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    fluidRow(
      column(
        width = 12,
        align = "center",
        h3("Data Input", style="font-weight: bold")
      )
    ),
    # h3("Data input"),
    matrixUI("matrix1"),
    br(),
    geneSetsUI("genes1"),
    br(),
    radioButtons("arg", "Change default settings:",
                 c("No" = "no",
                   "Yes" = "yes")),
    br(),
    fluidRow(
      actionButton("button", "Run"),
      downloadUI("download"),
      closeBtnUI("close"),
    )
  ),
  
  dashboardBody(
    shinyjs::useShinyjs(),
    add_busy_spinner(spin = "double-bounce", position = "bottom-right",
                     height = "100px", width = "100px"),
    fluidRow(
      box(
        width = 9,
        tabsetPanel(id = "Panels", type="tabs",
                    tabPanel("Samples",
                             textOutput("errorsGsva"),
                             htmlOutput("text1"),
                             plot1_UI("plot1"),
                             tableOutput("result")
                    ),
                    tabPanel("GeneSets",
                             uiOutput("text2"),
                             htmlOutput("text3"),
                             plot2_UI("plot2"),
                             plot3_UI("plot3")
                    ),
                    tabPanel("Session Info",
                             verbatimTextOutput("sessionInfo"))
        )
      ),
      box(
        width = 3,
        argumentsDataUI("argumentsInput")
      )
    )
  )
  
)

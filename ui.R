
library(shiny)
library(shinythemes)
library(shinydashboard)
library(zoo)
library(xts)
library(quantmod)
library(quadprog)
library(PortfolioAnalytics)
library(PerformanceAnalytics)
library(ROI)
library(ggplot2)
library(plotly)
library(dplyr)
library(readxl)


# "Portfolio Optimization"

header = dashboardHeader(title = tags$img(src = "Logo.png", height = "40", width = "140"), 
    dropdownMenu(type = "messages", 
                 messageItem(
                     from = "Andrés Felipe Cortés Bello",
                     message = "Click here to go to LinkedIn!",
                     href = "https://www.linkedin.com/in/felipe-cortes-bello/",
                     icon = icon("linkedin", "fa-2x")),
                 messageItem(
                     from = "Andrés Camilo Chacón Briceño",
                     message = "Click here to go to LinkedIn!",
                     href = "https://www.linkedin.com/in/andrescachaconb/",
                     icon = icon("linkedin", "fa-2x")),
                 icon = icon("linkedin", "fa-2x")
                 
                 ),
    
    dropdownMenu(type = "messages", messageItem(
        from = "felipecortesbusiness@gmail.com",
        message =  "Feedback and suggestions",
        icon = icon("envelope")),
        messageItem(
            from = "andreschabri@gmail.com",
            message =  "Feedback and suggestions",
            icon = icon("envelope"))
        , icon = icon("book", "fa-2x")
    )
    )



sidebar = dashboardSidebar(sidebarMenu(
    menuItem("Main Panel/ About us",
             tabName = "Main",
             icon = icon("archive")),
    menuItem("Optimization",
             tabName = "Optim",
             icon= icon("chart-pie"),
             menuSubItem("Analysis", tabName = "Ana"),
             menuSubItem("Efficient Frontier", tabName = "Effi")),
    menuItem("Performance",
             tabName = "Perfor",
             icon = icon("chart-line")),
    menuItem("Black- Litterman",
             tabName = "Litter",
             icon= icon("hand-holding-usd")),
    menuItem("Random Portfolio",
             tabName = "Random",
             icon = icon("balance-scale-right"))
))



body = dashboardBody(tags$script(HTML("$('body').addClass('fixed');")),
    tabItems(
        tabItem(tabName = "Main", h1("Portfolio Optimization", align = "center"), h2("Welcome!"), br(), h4(
            "This app summarizes
            different portfolio methodologies based on several
            metrics. You will find portfolios such as the tangency
            portfolio, Sortino portfolio, Treynor portfolio, Omega
            portfolio and minimum variance portfolio (MVP). 
            Considering that this project is based on the nature of Shiny
            reactivity, now we will provide a couple of dynamic options that
            you can adjust according to your preferences' profile:"), br(), h4(strong("
                Risk Free")),
                sliderInput(inputId = "rf", "Select the risk free rate", min = 0, max = 0.02, value = 0.011),
                h4(strong("
            Data range")), h6("Warning: You must select a data range between
            2008-01-01 and 2020-10-30"), dateRangeInput("dates", h6(strong("Choose date range:")),
                                                        start = "2008-01-01",
                                                        end = "2020-10-30"), br(), h3("
            Stock picking"), br(), h4("
            We selected 20 assets from de Nasdaq 100 index. We took this index
            because of its particular behavior during the last Covid- 19 crisis
            in which the Nasdaq broke historical records. The selection of the 20
            stocks enjoy the principle of diversification given these industries: 
            IT, Healthcare, Industrials, Consumption, Essentials and Utilities. With
            the first filter over the table, then we decided to select the top perform
            companies based on 2 elements: high ROE and historical performance."), br(), h4("
            The stocks are: Citrix Systems, Qualcomm, Apple Inc., IDEXX Laboratories Inc., 
            Align Technology, Amgen, Copart, Fastenal, Cintas Corporation, O’Reilly Auto Parts, 
            eBay Inc., Booking Holdings Inc., Netflix, Electronic Arts, Alphabet Inc Class A,
            PepsiCo, Monster Beverage, Costco, Xcel Energy Inc., Exelon Corporation.", br(), h2("Authors", align = "center"),
                                                                                            
                                                                                            HTML('<center><img src="Back.png" width="1000"></center>')
                                                                                            
                                                                                            
                                                                                            
                                                                                            )
            
        ),
        tabItem(tabName = "Ana", h1("Portfolio Optimization", align = "center"), tabBox(width = 12, 
            title = h5("Select Portfolio Methodology"),
            tabPanel("Tangency", h3("Tangency Portfolio"), 
                     
                    
                     fluidRow(
                         column(width = 4,
                                box(width = NULL,
                                    infoBoxOutput("RPT", width = "100%"),
                                    infoBoxOutput("SIGMAPT", width = "100%"),
                                    infoBoxOutput("Sharpe", width = "100%"),
                                    infoBoxOutput("MDS", width = "100%")  
                                         
                                    )
                         ),
                         
                         column(width = 8,
                                box(width = NULL, h1("Weights Tangency Portfolio", align = "center"),
                                    plotlyOutput("PIET")        
                                )
                         )
                     )
                     
                     
                     ),
            tabPanel("Treynor", h3("Treynor Portfolio"),
                     
                     
                     fluidRow(
                         column(width = 4,
                                box(width = NULL,
                                    infoBoxOutput("RTR", width = "100%"),
                                    infoBoxOutput("SIGMATR", width = "100%"),
                                    infoBoxOutput("SharpeTR", width = "100%"),
                                    infoBoxOutput("MDSTR", width = "100%")   
                                    
                                )
                         ),
                         
                         column(width = 8,
                                box(width = NULL, h1("Weights Treynor Portfolio", align = "center"),
                                    plotlyOutput("PIETR")         
                                )
                         )
                     )
                     
                     ),
            tabPanel("Sortino", h3("Sortino Portfolio"),
                     
                     
                     fluidRow(
                         column(width = 4,
                                box(width = NULL,
                                    infoBoxOutput("RSOR", width = "100%"),
                                    infoBoxOutput("SIGMASOR", width = "100%"),
                                    infoBoxOutput("SharpeSOR", width = "100%"),
                                    infoBoxOutput("MDSOR", width = "100%")  
                                    
                                )
                         ),
                         
                         column(width = 8,
                                box(width = NULL, h1("Weights Sortino Portfolio", align = "center"),
                                    plotlyOutput("PIESOR")        
                                )
                         )
                     )
                    
                     
                     ),
            tabPanel("Omega", h3("Omega Portfolio"),
                     
                     fluidRow(
                         column(width = 4,
                                box(width = NULL,
                                    infoBoxOutput("ROME", width = "100%"),
                                    infoBoxOutput("SIGMAOME", width = "100%"),
                                    infoBoxOutput("SharpeOME", width = "100%"),
                                    infoBoxOutput("MDOMEE", width = "100%")  
                                    
                                )
                         ),
                         
                         column(width = 8,
                                box(width = NULL, h1("Weights Omega Portfolio", align = "center"),
                                    plotlyOutput("PIEOME")        
                                )
                         )
                     )
                     
                     )
            
            
            
            
        )),
        tabItem(tabName = "Effi", h1("Efficient Frontier", align = "center"), 
                                                                   
                                                                  fluidRow(h4("Here you will find two efficient frontiers. We decided to separte them so
                                                                             as to make possible for you to analize them individually. The first graph 
                                                                              shows the 5 main portfolios (Tangecy Portfolio, Sortino Portfolio, Treynor
                                                                              Portfolio, Omega Portfolio and Minimum Variance Portfolio). The second graph has the 20 selected assets
                                                                              within a risk- reward plane. Feel free to interact with the graphs!"),
                                                                       box(width = 6,
                                                                           plotlyOutput("EFP")
                                                                           
                                                                       ),
                                                                       box(width = 6,
                                                                         
                                                                           plotlyOutput("EFA")
                                                                       )
                                                                       
                                                                       )

        ),
        tabItem(tabName = "Perfor", h1("Portfolio Performance", align = "center") ,tabBox(width = 12,
                                           title = "Select sample",
                                           tabPanel("In Sample", h2("Performance in-sample"),
                                                    
                                                    sidebarLayout(
                                                        sidebarPanel(
                                                            radioButtons(inputId = "Bench", label = "Select the benchmark you prefer:", choices = c("Nasdaq" = "^IXIC","S&P500" = "^GSPC", "Dow Jones" = "^DJI")),
                                                            infoBox(title = "Omega Portfolio", subtitle = "Best Performance", icon = icon("trophy"), width = "100%", fill = T, color = "green"),
                                                            infoBox(title = "Benchmark", subtitle = "Worst Performance", icon = icon("battery-quarter"), width = "100%", fill = T, color = "red")                                                
                                                        ),
                                                        
                                                        mainPanel(plotlyOutput("PERIN"))
                                                        
                                                    )
                                                    
                                                    
                                                    ),
                                           tabPanel("Out sample", h2("Performance out-sample"),
                                                        
                                                        plotlyOutput("PEROUT")
                                                    
                                                    ))),
                                          
        tabItem(tabName = "Litter", h1("Black- Litterman model", align = "center"),
                
                h4("Showing up next, you will find 3 elements. First and below, there are the weights of the traditional
                   mean- variance model and also the weights of the Black Litterman's optimization. Then you will find the comparison in terms
                   of performance between the same two models along with the possibility to select the initial investment and the views' explanation. 
                   Finally, in case you are interested, you will have access to the web page of investors relations of each one of the companies that we 
                   selected for view purposes. Just click on the box at the end of the page!"),
                
                fluidRow(
                    column(width = 12,
                           box(width = 6, h3("Weights Equilibrium Model", align = "center"),
                               plotlyOutput("EM")        
                           ),
                           box(width = 6, h3("Weights Black Litterman", align = "center"),
                               plotlyOutput("BLM")        
                           )
                           
                    )
                ),
                fluidRow(
                    
                    column(width = 12,
                           box(width = 4,
                               selectInput(inputId = "value", label = "Select first value index:", 
                                           choices = c(1, 10, 100, 1000)),
                               h3("Views explanation", asign = "center"), br(), p("- Apple, expected view: 4.4%. iPhone Availability Tracker: Lead Times Moderate for iPhone 12; But Stable for iPhone 12 Pro."), 
p("- Xcel Energy, expected view: -6.6%. Meet the New Boss, Same as the Old Boss: Strong 3Q20 Results, No Surprises on Five-Year Plan Roll-Forward"), p("- O'Reilly Automotive, expected view: 12.1%. Strong Results But 2021 Remains a Show-Me; Management Takes."),
p("- Cintas, expected view: -7.6%. Strong 1QF21 (Aug) with narrowing rev declines, though 2Q expectations are a bit fuzzy"), p(strong("Source: Bloomberg Terminal"))
                           ),
                           
                           box(width = 8, h3("Performance Black Litterman Model", align = "center"),
                               plotlyOutput("PerBL")        
                           )
                    )
                ),

fluidRow(
    box(width = 12,
        infoBox(title = "Apple", subtitle = "Click here to access to investor relations", width = 3, fill = T, col = "navy", icon = icon("apple"), href = "https://investor.apple.com/investor-relations/default.aspx
"),
        infoBox(title = "XCEL Energy", subtitle = "Click here to access to investor relations", width = 3, fill = T, col = "lime", icon = icon("broadcast-tower"), href = "https://www.xcelenergy.com/stateselector?stateSelected=true&goto=%2Fcompany%2Finvestor_relations
"),
        infoBox(title = "O'Reilly", subtitle = "Click here to access to investor relations", width = 3, fill = T, col = "aqua", icon = icon("car"), href = "https://corporate.oreillyauto.com/investor-relations-general-info
"),
        infoBox(title = "CINTAS", subtitle = "Click here to access to investor relations", width = 3, fill = T, col = "red", icon = icon("tshirt"), href = "https://www.cintas.com/investors/
")
    )
)

                
                ),
        tabItem(tabName = "Random", h1("Random Portfolios", align = "center"),
                
                h4("In this section we present every single thing related to random portfolios and its optimization.
                   There are two graphs, one of them is the Risk- Reward plane with the optimum portfolio in red. There are also some boxes set aside that
                   provide useful information. The other graph includes 
                   the weights of the optimum portfolio. Enjoy!"),
                
                fluidRow(
                    
                    column(width = 12,
                           box(width = 8, h3("Efficient Frontier Random Portfolios", align = "center"),
                               plotlyOutput("EFRP")
                                ),
                           
                           box(width = 4, h3("Measures", align = "center"),
                               
                               infoBox(title = "Return O.Portfolio", value = paste(round(0.07350797*100, 3), "%"), icon = icon("chart-line"), fill = T, col = "lime", width = "100%"),
                               infoBox(title = "Risk O.Portfolio", value = paste(round(0.1317854*100, 3), "%"), icon = icon("wave-square"), fill = T, col = "blue", width = "100%"),
                               infoBoxOutput("SRPA", width = "100%")
                           )
                    )
                ),
                fluidRow(
                    
                    column(width = 12,
                           box(width = 4, h3("Most Weighted Assets", align = "center"),
                               infoBox(title = "PepsiCo", value = paste(48.4, "%"), icon = icon("cocktail"), fill = T, col = "navy", width = "100%"),
                               infoBox(title = "Exelon", value = paste(23.6, "%"), icon = icon("broadcast-tower"), fill = T, col = "red", width = "100%"),
                               infoBox(title = "Electronic Arts", value = paste(8.2, "%"), icon = icon("xbox"), fill = T, col = "green", width = "100%")
                               
                           ),
                           
                           box(width = 8, h3("Weights Optimum Random Portfolio", align = "center"),
                               plotlyOutput("WRP")
                           )
                    )
                )
                
                
                
                )
    ))



ui <- dashboardPage(header = header, sidebar = sidebar, body = body, skin = "green")

    



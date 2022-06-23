library(shiny)
library(shinyMatrix)
library(xtable)
library(shinyBS)
library(knitr)
library(DeltaMAN)




generate.matrix <- function(k) {

  n = 100
  a = runif(k)
  prob = a/sum(a);prob
  gold = sample(1:k, n, replace = T, prob = prob)
  
  
  # seleccionar elementos que el observador C reconoce (e identifica correctamente)
  d = 68
  id = sort(sample(1:n, d))
  
  random_id = c(1:n)[-id]
  responses_random = sample(1:k, length(random_id), replace = T)
  
  
  responses = gold
  responses[random_id] = responses_random
  A = table(gold, responses)
  suppressMessages({
    A = DeltaMAN:::getM1(A)
  })
  
  A = matrix(A, ncol = k);A
  
  return(A)
}

html_title <-
  '<span class="logo">
    <div style="display:inline-block;">
      <a href="https://www.ugr.es/~bioest/software/delta/index.php"><img src="DELTA.gif" height="35"/></a>
      <b style="font-size:22px;font-family:Tahoma; color: #4682B4" >Delta</b> 
      <div style="font-size:12px;font-family:Tahoma; color: #4682B4; margin-left: 40px;">Measure of nominal agreement between two raters </div>
    </div>
  </span>'

# Define UI for application that draws a histogram
ui <- navbarPage(
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
    # Application title
    title = HTML(html_title),

    header = tags$style(
      ".navbar-static-top {height:75px}"
      ),

    footer = includeHTML("./www/footer.html"),
    windowTitle= tags$head(
      tags$link(rel = "icon", type = "image/png", href = "DELTA.gif"),
      tags$title("DELTA")
    ),

  #   tags$head(
  #     tags$style(HTML("
  #                   div.MathJax_Display{
  #                   text-align: left !important;
  #                   }
  # "))
  #   ),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(

          fluidRow(
            # column(3, HTML('<b>A1</b>')),
            column(4, numericInput("bins","Matrix order", min = 2, max = NA, value = 3, step=1))
          ),
            matrixInput("matrizA", label="Classification matrix",
                        value =  matrix("", ncol = 3, nrow = 3), 
                        rows = list(names = F),
                        cols = list(names = F), 
                        class = 'numeric'
            ),
          actionButton("random", "Get example"),            
          actionButton("delete", "Delete data"),
          br(),
          br(),
          bsCollapse(id = "collapsePanels", open = "Sampling",
                     bsCollapsePanel("Sampling",
                                     checkboxInput('standard', "The rater on the rows is a gold standard"),
                                     checkboxInput('fixedRows', "Row marginals have been set beforehand")
                                     , style = "info"),
                     bsCollapsePanel("Numeric procedure",
                                     fluidRow(
                                       column(width = 4,
                                         div(style = "white-space: nowrap;", 
                                             h5('Tolerance: 1.0E-',style="display:inline-block"),
                                             div(style="display: inline-block; width: 100%;",
                                                 numericInput("tol",label = NULL, min = 1, max = NA, value = 7, step=1))))
                                       ),
                                     fluidRow(
                                       column(width = 4,
                                         div(style = "white-space: nowrap;", 
                                             h5('Max.  iterations :',style="display:inline-block"),
                                             div(style="display: inline-block; width: 100%;",
                                                 numericInput("mxits",label = NULL, min = 0, max = NA, value = 100, step=1))))
                                       )
                                     , style = "primary"),
                     bsCollapsePanel("Output",
                                     fluidRow(
                                       column(width = 3,
                                              div(style = "white-space: nowrap;", 
                                                  h5('Significant digits:',style="display:inline-block"),
                                                  div(style="display: inline-block; width: 100%;",
                                                      numericInput("digits",label = NULL, min = 0, max = NA, value = 3, step=1))))
                                     ),
                                     checkboxInput('fullReport', "Full report"),
                                     radioButtons('reportFormat', 'Format', choices = c('.pdf', '.tex'), inline = T),
                                     downloadButton("report", "Download report")
                                     , style = "primary")
          ),

        ),

        # Show the report
        mainPanel(
          uiOutput("outputMatrix"),
          uiOutput("outputSampling"),
          uiOutput("ND"),
          uiOutput("AD")


        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  
  output$ND<- renderUI({
    
    bsCollapse(#open = "Standard analysis",
               bsCollapsePanel("Standard analysis", style = "warning",
                               bsCollapse(open = "Summary: Goodness of fit, Kappa and Delta", multiple = T,
                                          bsCollapsePanel("Summary: Goodness of fit, Kappa and Delta", style = "primary",
                                                          uiOutput("Summary_DN"),
                                                          bsCollapse(
                                                            bsCollapsePanel("Expected frequencies", style = "info",
                                                                            uiOutput("expected_freq")
                                                            )
                                                          )
                                          ),
                                          bsCollapsePanel("Parameters and measures (SE) of concordance", style = "primary",
                                                          verbatimTextOutput("validMeasures_DN"),
                                                          bsCollapse(
                                                            bsCollapsePanel("All parameters and measures (SE) of concordance", style = "info",
                                                                            verbatimTextOutput("allMeasures_DN")
                                                            )
                                                          )
                                          ),
                                          bsCollapsePanel("Full report", style = "primary",
                                                          bsCollapse(
                                                            # bsCollapsePanel("All parameters and measures of concordance", style = "info",
                                                            #                 verbatimTextOutput("allMeasures_DN")
                                                            # ),
                                                            bsCollapsePanel("Determination of B0 and B", style = "info",
                                                                            verbatimTextOutput("B0")
                                                            ),
                                                            bsCollapsePanel("Variances and covariances", style = "info",
                                                                            verbatimTextOutput("covs")
                                                            )
                                                          )
                                                          
                                                          
                                          )
                               )
               )
    )

    
  })
  
  
  output$AD<- renderUI({
    if(orden()==2){
      bsCollapse(#open = "Asymptotic analysis",
                 
                 bsCollapsePanel("Asymptotic analysis",style = "warning",
                                 bsCollapse(open = c("ASYMPTOTIC SOLUTION I","ASYMPTOTIC SOLUTION II"), multiple = T,
                                            bsCollapsePanel("ASYMPTOTIC SOLUTION I", style = "primary",
                                                            helpText("Solution based on the original data."),
                                                            verbatimTextOutput("asymtotic1"),
                                                            bsCollapse(
                                                              bsCollapsePanel("All parameters and measures (SE) of concordance", style = "info",
                                                                              verbatimTextOutput("allMeasures_DA1")
                                                              )
                                                            )
                                            ),
                                            bsCollapsePanel("ASYMPTOTIC SOLUTION II", style = "primary",
                                                            helpText("Solution based on the original data incremented + 1."),
                                                            verbatimTextOutput("asymtotic2"),
                                                            bsCollapse(
                                                              bsCollapsePanel("All parameters and measures (SE) of concordance", style = "info",
                                                                              verbatimTextOutput("allMeasures_DA2")
                                                              )
                                                            )
                                            )
                                 )

                 )
      )
    }
    
  })
  
  
  # INPUTS ----
  
    orden = reactive({
      input$bins
    })
    
    tol = reactive({
      as.numeric(eval(paste0("1e-",input$tol)))
    })
    
    mxits = reactive({
      input$mxits
    })
    
    standard <- reactive({
      input$standard
    })
    
    fixedRows <- reactive({
      input$fixedRows
    })
    
    fullReport <- reactive({
      input$fullReport
    })
    
    digits <- reactive({
      input$digits
    })
  
    
  
    # CLASSIFICATION MATRIX ----
    matriz = eventReactive(input$matrizA, {
      M <- input$matrizA
      
      dimnames(M) = list(LETTERS[1:orden()], LETTERS[1:orden()])
      
      if(!any(is.na(M))) {
        return(M)
      }
      else {
        return(NULL)
      }
    })
    
    
    # UPDATE MATRIX to the new order ----
    observeEvent(input$bins, {
      A = matrix("", ncol = input$bins, nrow = orden())
      updateMatrixInput(session, inputId = "matrizA", value = A)
      
    }, ignoreInit = TRUE)
  
    
    

    # Get a random matrix of given order
    observeEvent(input$random, {
      A = generate.matrix(orden())
      updateMatrixInput(session, inputId = "matrizA", value = A)
      
    }, ignoreInit = TRUE)
  
    

    observeEvent(input$delete, {
      updateMatrixInput(session, inputId = "matrizA", value = matrix("",nrow=orden(),ncol=orden()))
      #M <- matriz() 
    }, ignoreInit = TRUE)
    

    # DELTA RESULTS -----
    delta_obj <- reactive({
      M <- matriz()
      validate(
        need(!is.null(M) && all(!is.na(M)), '')
      )
      delta(M, standard = standard(), fixedRows = fixedRows(), tol = tol(), mxits = mxits())
    })
    
    # OUTPUTS -----
    # .. Matrix ----
    output$outputMatrix <- renderUI({
      M <- matriz() 
      validate(
        need(!is.null(M) && all(!is.na(M)), 'Empty matrix')
      )

      M = cbind(M, rowSums(M))
      M = rbind(M, colSums(M))
      
      dimnames(M) = lapply(dimnames(M), function(x) ifelse(x == '', 'sum',x))
      
      if (!is.null(M) && !any(is.na(M))) {
        Mlatex <- print(xtable(M, align=c("c",rep("c", ncol(M)-1), "|c"), auto = TRUE),
                   floating=FALSE, tabular.environment="array", comment=FALSE, print.results=FALSE,
                   hline.after=c(-1, 0,nrow(M)-1,nrow(M)))
        list(
          withMathJax(HTML(Mlatex))
        )
      }
    })
    


    

    # .. Sampling ----
    output$outputSampling <- renderUI({
      
      obj = delta_obj()
      # cat('SELECTED MODEL CONDITIONS\n')
      
      if(obj$standard){
        txt1 = 'There is a gold standard classification by rows.'
        # txt2 = '   The rater by rows is a gold standard.\n'
      }else{
        txt1 = 'None of the raters is a gold standard.'
        # txt2 = '   None of the raters is a gold standard.\n'
      }
      
      if(obj$fixedRows){
        txt2 = 'Row marginal frequencies have been set beforehand (type II sampling).'
        # txt1 = '   Measures (SE) under type-II sampling.\n'
      }else{
        txt2 = 'Row marginal frequencies have not been set beforehand (type I sampling).'
        # txt1 = '   Measures (SE) under type-I sampling.\n'
      }
      txt = paste("Selected model conditions: <ul><li>",txt1, "</li><li>",txt2,"</li></ul>")
      HTML(txt)
    })
    
    # .. Standard analysis ----
    output$Summary_DN <- renderUI({
      obj = delta_obj()
      
      chi_val = DeltaMAN:::.round(obj$GOF$statistic, digits())
      chi_df = obj$GOF$df
      chi_pval = DeltaMAN:::.round(obj$GOF$pval, digits())

      Kappa_val = DeltaMAN:::.round(obj$Kappa$Estimate, digits())
      Kappa_SE = DeltaMAN:::.round(obj$Kappa$SE, (digits()+1))
      Delta_val = DeltaMAN:::.round(obj$Delta$Estimates$Delta, digits())
      Delta_SE = DeltaMAN:::.round(obj$Delta$SE$Delta, (digits()+1))

      if(attr(obj$GOF, 'validity')==T){
        chi2 = paste0("$$\\chi^2 =",chi_val,",\\text{ df =  }", chi_df,",\\text{ p-value =}",chi_pval,"$$")
      }else{
        x = gsub('>20%', '$>20\\\\%$', attr(obj$GOF, 'validity'))
        chi2 = paste0("$$\\chi^2 \\text{ ",tolower(x),"}","$$")
      }

      withMathJax(helpText(
        paste0(chi2),
        paste0("$$\\kappa~(SE) = ",Kappa_val,"~(",Kappa_SE,")$$"),
        paste0("$$\\Delta~(SE) = ",Delta_val,"~(",Delta_SE,")$$")
      ))
    })
    
    output$expected_freq <- renderUI({
      obj = delta_obj()
      Mexp = obj$GOF$Expected
      
      Mlatex <- print(xtable(Mexp), 
                      floating=FALSE, tabular.environment="array", comment=FALSE, print.results=FALSE)
      list(
        withMathJax(HTML(Mlatex))
      )
    })
    
    
    # output$outputMatrix <- renderUI({
    #   M <- matriz() 
    #   validate(
    #     need(!is.null(M) && all(!is.na(M)), 'Empty matrix')
    #   )
    #   
    #   M = cbind(M, rowSums(M))
    #   M = rbind(M, colSums(M))
    #   
    #   dimnames(M) = lapply(dimnames(M), function(x) ifelse(x == '', 'sum',x))
    #   
    #   if (!is.null(M) && !any(is.na(M))) {
    #     Mlatex <- print(xtable(M, align=c("c",rep("c", ncol(M)-1), "|c"), auto = TRUE),
    #                     floating=FALSE, tabular.environment="array", comment=FALSE, print.results=FALSE,
    #                     hline.after=c(-1, 0,nrow(M)-1,nrow(M)))
    #     list(
    #       withMathJax(HTML(Mlatex))
    #     )
    #   }
    # })
    
    output$validMeasures_DN <- renderPrint({
      obj = delta_obj()
      print(obj$Delta, digits = digits(), tex = F)
    })
    
    
    output$allMeasures_DN <- renderPrint({
      obj = delta_obj()
      print(obj$all.measures, digits = digits(), tex = F)
    })
    
    
    output$B0 <- renderPrint({
      obj = delta_obj()
      print(obj$problem.parameters, digits = digits(), tex = F)
    })
    
    
    output$covs <- renderPrint({
      obj = delta_obj()
      print(obj$cov, digits = digits(), tex = F)
    })
    

    # .. Asymptotic ----
    output$asymtotic1 <- renderPrint({
      
      obj = delta_obj()
      
      validDA = obj$asymptoticDelta$validMeasures

      print(validDA$DA_A1, digits = digits(), tex = F)

      
    })
    
    
    output$allMeasures_DA1 <- renderPrint({
      
      obj = delta_obj()
      
      allDA = obj$asymptoticDelta$allMeasures
      
      print(allDA$DA_A1, digits = digits(), tex = F)
      
    })
    
    
    output$asymtotic2 <- renderPrint({
      
      obj = delta_obj()
      
      validDA = obj$asymptoticDelta$validMeasures

      print(validDA$DA_A2, digits = digits(), tex = F)

    })
    
    output$allMeasures_DA2 <- renderPrint({
      
      obj = delta_obj()
      
      allDA = obj$asymptoticDelta$allMeasures
      
      print(allDA$DA_A2, digits = digits(), tex = F)
      
    })
    
    format = reactive({
      format = input$reportFormat
      ifelse(format == '.pdf', 'pdf', 'zip')
    })
    
    filename = reactive({
      format = input$reportFormat
      ifelse(format == '.pdf', 'myreport.pdf', 'myreport.zip')
    })
    # DOWNLOAD REPORT ----
    # output$report = downloadHandler(
    # 
    #   filename = paste0('myreport.',format()),
    #   
    #   content = function(file) {
    #     # out = knit2pdf('input.Rnw', clean = TRUE)
    #     Sweave2knitr("input.Rnw")
    #     knit('input-knitr.Rnw', output = paste('input.tex', sep = ''))
    #     
    #     if(format() == 'pdf'){
    #       tools::texi2pdf('input.tex')
    #       # file.rename(out, file) # move pdf to file for downloading
    #       file.copy('input.pdf', file, overwrite = TRUE)
    #     }else{
    #       zip::zipr(zipfile=file, files='input.tex')
    #     }
    #   },
    #   
    #   contentType =paste0('application/',format())
    # )
    
    
    output$report = downloadHandler(
      
      # filename = 'myreport.zip',
      filename = function(){filename()},
      content = function(file) {
        # out = knit2pdf('input.Rnw', clean = TRUE)
        Sweave2knitr("input.Rnw")
        knit('input-knitr.Rnw', output = paste('report.tex', sep = ''))
        if(format() == 'pdf'){
          tools::texi2pdf('report.tex')
          file.copy('report.pdf', file, overwrite = TRUE)
        }else{
          zip::zipr(zipfile=file, files='report.tex')
        }
        
      },
      
      contentType = NULL
    )
}

# Run the application 
shinyApp(ui = ui, server = server)


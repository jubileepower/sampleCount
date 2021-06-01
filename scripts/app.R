#############################################
### Draft app for sample size calculations in R:
### https://jubileepower.shinyapps.io/sampleCount/
### Taking the sample counts of provinces (2nd column),
### calculate and output the numbers of samples required for surveillance
### The unit entry is used for dividing counts.
### Eg. To calculate weekly average using counts from one month, unit entry should be 4.4 (31/7).
### author: "Julie Chih-yu Chen"
### date: "28/05/2021"
#############################################


library(shiny)
library(shinysky)

######
# sample size for detection
######
sampDetectFcn<-function(Ms, p, d=0.5, crit=1.96){
  n <- ceiling(crit^2*p*(1-p)/(p*d)^2)
  ceiling((n*Ms)/(n+Ms))
}

######
# sample size for difference
######
sampDiffFcn<-function(Ms, p1, p2, d=0.5, critBeta=0.87, critAlpha=1.64){
  n <- ceiling((critAlpha+critBeta)^2*(p1*(1-p1)+p2*(1-p2))/(p1-p2)^2)
  ceiling((n*Ms)/(n+Ms))
}


######
# Shiny
######
ui  <- shinyUI (       wellPanel(hotable("matrixTable"),hr(),
                                 numericInput("unit", label = h3("Unit"), value = 1),
                                 hr(),dataTableOutput("resultTable"))     )


server <- shinyServer (function(input, output) {
  A = cbind(rep(1000,13), rep(100/13,13)) # init - input matrix A
  colnames(A)=c("Count","Percent")
  prov=c("Alberta","British Columbia","Manitoba","New Brunswick","Newfoundland and Labrador","Northwest Territories","Nova Scotia","Nunavut","Ontario","Prince Edward Island","Quebec","Saskatchewan","Yukon")

  output$matrixTable <- renderHotable({data.frame(prov,A)}, readOnly = FALSE)
  unit <- renderText({ input$unit })

  R = matrix(rep(0,52), nrow=13) # init - result matrix R
  output$resultTable <- renderDataTable({data.frame(R)})

  observe({  # process matrix
    df <- hot.to.df(input$matrixTable)
    if(!is.null(df)) {    # ensure data frame from table exists
      print(df)
      B = data.matrix(df) # ensure its numeric

      ## taking the counts frmo provs
      Ms=ceiling(B[,2]/input$unit)

      ## sample size for detection
      ncorr25 <- sampDetectFcn(Ms, p=0.025)
      ncorr1 <- sampDetectFcn(Ms, p=0.01)

      ## sample size for detecting a difference
      ncorrDiff_2.5_5 <- sampDiffFcn(Ms, p1=0.025, p2=0.05 )
      ncorrDiff_1_3 <- sampDiffFcn(Ms, p1=0.01, p2=0.03 )


      R = data.frame(prov=df[,1],counts=Ms, toSamp_2.5=ncorr25,toSamp_2.5_5= ncorrDiff_2.5_5, toSamp_1=ncorr1, toSamp_1_3=ncorrDiff_1_3)
      output$resultTable <- renderDataTable({R})
    }
  }) # end of observe
}) # end of server

shinyApp(ui, server)

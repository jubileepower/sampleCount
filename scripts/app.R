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
ui  <- shinyUI (       wellPanel( p("Paste in the number of SARS-CoV-2 case counts below"),
                       hotable("matrixTable"),hr(),
                       numericInput("unit", label = h3("Unit"), value = 1),
                       p("If you pasted counts from 4 weeks above and prefer the calculation of weekly required samples, modify the unit to 4."),
                       hr(),h4("Calculations for representative sampling of SARS-CoV-2 cases for genomic monitoring from routine surveillance"),
                       tableOutput("resultTable"),hr(),
                       h4("The summary table of samples required in ranges of cases"),
                       tableOutput("summaryTable")))

server <- shinyServer (function(input, output) {
  #A = cbind(rep(1000,13), rep(100/13,13)) # init - input matrix A
  #colnames(A)=c("Count","Percent")
  A=rep(1000,13)# init - input matrix A
  prov=c("Alberta","British Columbia","Manitoba","New Brunswick","Newfoundland and Labrador","Northwest Territories","Nova Scotia","Nunavut","Ontario","Prince Edward Island","Quebec","Saskatchewan","Yukon")

  output$matrixTable <- renderHotable({data.frame(Province=prov,Count=A)}, readOnly = FALSE)
  unit <- renderText({ input$unit })

  R = matrix(rep(0,13), nrow=13) # init - result matrix R
  output$resultTable <- renderTable({data.frame(R)})

  observe({  # process matrix
    df <- hot.to.df(input$matrixTable)
    if(!is.null(df)) {    # ensure data frame from table exists
      print(df)
      B = data.matrix(df) # ensure its numeric
      print(B)
      ## taking the counts frmo provs
      Ms=ceiling(B[,2]/input$unit)

      ## sample size for detection
      ncorr25 <- sampDetectFcn(Ms, p=0.025)
      ncorr1 <- sampDetectFcn(Ms, p=0.01)

      ## sample size for detecting a difference
      ncorrDiff_2.5_5 <- sampDiffFcn(Ms, p1=0.025, p2=0.05 )
      ncorrDiff_1_3 <- sampDiffFcn(Ms, p1=0.01, p2=0.03 )


      R <- data.frame(Province=df[,1],Count=Ms, toSamp_2.5=ncorr25, toSamp_1=ncorr1, toSamp_2.5_5= ncorrDiff_2.5_5, toSamp_1_3=ncorrDiff_1_3)
      R[] <- lapply(R, as.character)
      colnames(R)=c("Province","Counts","Detect_@_2.5%", "Detect_@_1%", "Detect_Δ2.5%–5%", "Detect_Δ1%–3%")

      output$resultTable <- renderTable({R})

      ###summary table
      maxCase=max(Ms) ### max count of COVID cases this week

      MssAll<-c(100000, 50000, 25000, 10000, 5000, 2500, 1000, 500 )
      MssAllName<-c("50001-100000", "25001-50000", "10001-25000", "5001-10000", "2501-5000", "1001-2500", "501-1000", "<500" )

      ### the rows to return according to the max case count
      Mss <- MssAll[max(which(MssAll>maxCase)):length(MssAll)]
      MssName <- MssAllName[max(which(MssAll>maxCase)):length(MssAll)]

      summaryTable<-data.frame(rangeTop=MssName, at2.5=sampDetectFcn(Mss, p=0.025), at1=sampDetectFcn(Mss, p=0.01), from2.5to5=sampDiffFcn(Mss, p1=0.025, p2=0.05), from1to3=sampDiffFcn(Mss, p1=0.01, p2=0.03))
      summaryTable[] <- lapply(summaryTable, as.character)
      colnames(summaryTable)=c("Number of positive SARS-CoV-2 cases","Detect @ 2.5%","Detect @ 1%","Detect Δ 2.5% – 5%","Detect Δ 1% – 3%")
      output$summaryTable <- renderTable({summaryTable})

    }
  }) # end of observe
}) # end of server

shinyApp(ui, server)

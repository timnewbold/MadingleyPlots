
.gridArrange<-function(nCells){
  arrange<-list('1'=c(1,1),
                '2'=c(1,2),
                '3'=c(1,3),
                '4'=c(2,2),
                '5'=c(2,3),
                '6'=c(2,3),
                '7'=c(3,3),
                '8'=c(3,3),
                '9'=c(3,3),
                '10'=c(4,3),
                '11'=c(4,3),
                '12'=c(4,3)
  )
  
  stopifnot(nCells %in% names(arrange))
  
  return(arrange[nCells][[1]])
}

.plotDims<-function(nCells){
  widths<-list('1'=8.5/2.54,
               '2'=12.5/2.54,
               '3'=17.5/2.54,
               '4'=12.5/2.54,
               '5'=17.5/2.54,
               '6'=17.5/2.54,
               '7'=17.5/2.54,
               '8'=17.5/2.54,
               '9'=17.5/2.54,
               '10'=17.5/2.54,
               '11'=17.5/2.54,
               '12'=17.5/2.54
  )
  heights<-list('1'=6/2.54,
               '2'=6/2.54,
               '3'=6/2.54,
               '4'=12/2.54,
               '5'=12/2.54,
               '6'=12/2.54,
               '7'=18/2.54,
               '8'=18/2.54,
               '9'=18/2.54,
               '10'=24/2.54,
               '11'=24/2.54,
               '12'=24/2.54
  )
  
  stopifnot(nCells %in% names(widths))
  stopifnot(nCells %in% names(heights))
  
  return(list(width=widths[nCells][[1]],height=heights[nCells][[1]]))
  
}
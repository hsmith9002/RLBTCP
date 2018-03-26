#' Plot Transcript Content
#' 
#' This function takes the output from the content.for.plot() function and plots the transcript content that was summarized. It generates 4 barcharts. 
#' A low and high resolution for the "At least one" criteria, and a low and high resolution for the "in all" criteria.
#' @param data The object created by the content.for.plot()
#' @keywords barchart
#' @export
#' @examples 
#' transcriptBarChart(content.summary)

transcriptBarChart <- function(content.plot.data, directory, filename) {
  # plot everything
  ALO.plot <- ggplot(content.plot.data$ALO.ForPlot, aes(Names, value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
    ggtitle("Transcript Content Using At Least 1 Criteria") + labs(y=paste("# of Transcripts \n Total: ", sum(content.plot.data$ALO.ForPlot$value), sep = ""), x = "Transcript Type")
  
  ALOSS.plot <- ggplot(content.plot.data$ALO.ForPlotSS, aes(Names, value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
    ggtitle("Transcript Content Using At Least 1 Criteria") + labs(y=paste("# of Transcripts \n Total: ", sum(content.plot.data$ALO.ForPlotSS$value), sep = ""), x = "Transcript Type")
  
  # plot everything
  ALL.plot <- ggplot(content.plot.data$ALL.ForPlot, aes(Names, value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
    ggtitle("Transcript Content Using In All Criteria") + labs(y=paste("# of Transcripts \n Total: ", sum(content.plot.data$ALL.ForPlot$value), sep = ""), x = "Transcript Type")
  
  ALLSS.plot <- ggplot(content.plot.data$ALL.ForPlotSS, aes(Names, value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
    ggtitle("Transcript Content Using In All Criteria") + labs(y=paste("# of Transcripts \n Total: ", sum(content.plot.data$ALL.ForPlotSS$value), sep = ""), x = "Transcript Type")
  pdf(file = paste(directory, "/", filename, ".pdf", sep = ""), height = 8, width = 10)
  multiplot(ALO.plot, ALOSS.plot, ALL.plot, ALLSS.plot, cols = 2)
  dev.off()
}
# --------------------------------------------------
# File Name: run.r
# Purpose:
# Creation Date: Jan 2016
# Last Modified: Wed Jan 27 12:40:45 2016
# Author(s): The DeepSEQ Team, University of Nottingham UK
# Copyright 2016 The Author(s) All Rights Reserved
# Credits: 
# --------------------------------------------------


# Config ...
graphics.off()
.SavedPlots <- NULL # Clear plots history
windows(record = TRUE) 

pdfout = 1

fname = 'out'

if (pdfout == 1) 
	pdf(file = 'plot_bw.pdf'
		, paper = "a4"
		, pointsize = 8
		)

# plot.new()
par(mfrow = c(4, 2))

#-------------------------------------------------------------------------------- 
# Plot labels

titles = c( 'Lambda Synthetic Reads', 'Lambda Amp3 Reads RU21')

normText = c("Un-normalised", "Normalised")

letters = c("A", "B", "C", "D", "E", "F", "G", "H")


#-------------------------------------------------------------------------------- 
# Read the data 

a<-read.table(file = fname, sep = ',', header = F)


colnames(a) <- c(
	"ReadFolder", "Normalise", "Chapter", "File", "Amplicon"
		, "winSz", "Len", "qrySz", "Offset"
		, "Read", "Ref", "Time_t", "trg", "Pos", "outBy", "read_t_Success"
		, "Read_", "Ref_", "Time_c", "trg_", "Pos_", "outBy"
		, "read_c_Success_", "dim", "quasi2d"
		)

# Total time taken on template and complement ...
a$Time <- a$Time_t + a$Time_c


#a <- a[ a$Pos ! = -1 , ]
#a$ReadFolder <- factor(a$ReadFolder)
#a$Normalise <- factor(a$Normalise)

summary(a)

#-------------------------------------------------------------------------------- 


# Loop through the folders and variations (normalised / un-normalised) ...
 
folders = rev(levels(a$ReadFolder))
print(folders)
 
normalisations = levels(a$Normalise)
print(normalisations)

 
for (plotType in c('Unstacked', 'Stacked')) {
	# Unstacked -- Plot only successful reads .. 
	# Stacked -- Plot both success and failed reads
	print(plotType)
 
 
	# NB R counts list elements starting from 1 ..
	j = 1 # Index into list of letters for labeling plots 
 
	for (n in normalisations) {
		for (f in folders) {
		
			if (f == "synthetic") numReads = 100 else numReads = 48
			print(numReads)
		
			folderNum = which(folders == f)
			print(f)
			normNum = which(normalisations == n)
			print(n)
		
			# Get relevant subset of data ...
			b <- a[ a$ReadFolder == f 
				& a$Normalise == n
				, ] 
		
			if (plotType == 'Unstacked') 
				b <- b[ b$read_t_Success == ' True', ] # NB space ..
		
			attach(b)
			summary(b)
		
			# Cross tabulate the data ....
			x0 <- xtabs(~winSz + Read + Ref + winSz)
			print(ftable(x0, col.vars = c(2, 3))) 
		
			#-------------------------------------------
			# Plot "Read Counts" ..
		
			# Set it up ...
			txt <- paste( letters[j], ") "
					, titles[folderNum] , ", "
					, normText[normNum]
					)
			print(txt)
		
		
			if (plotType == 'Unstacked') {
				xc <- xtabs(~winSz) 
				m <- (as.matrix(t(xc)) / numReads) * 100 # As % ..
				pal = 'grey'
				yLim = c(0, 110)
				yLab = 'Aligned Reads (%)'
			}
			else {
				xc <- xtabs(~winSz + read_t_Success) 
				m <- t(as.matrix(xc))
				m <- prop.table(m, margin = 2)
				pal = heat.colors(length(rownames(m)))
				yLim = c(0, 1)
				yLab = 'Proportion Reads'
			}
		
			# Do plot ...
			barplot(m
				, beside = F
				, col = pal
				, ylim = yLim
				, ylab = yLab
				, xlab = 'Window Size'
				, main = txt
				)
			grid()
			j = j + 1
			
			#-------------------------------------------
			# Plot "Read Times" ..
		
			# Set it up ...
			txt <- paste( letters[j], ") "
					, titles[folderNum] , ", "
					, normText[normNum]
					)
			print(txt)
		
			if (plotType == 'Unstacked') {
				xt <- xtabs(Time~winSz)
				m <- t(as.matrix(xt))
				pal = 'grey'
			}
			else {
				xt <- xtabs(Time~winSz + read_t_Success) 
				m <- as.matrix(t(xt))
				pal = heat.colors(length(rownames(m)))
			}
		
			# Do plot ...
			barplot(m
				, beside = F, 
				, col = pal
				, ylab = 'Time (seconds)'
				, xlab = 'Window Size'
				, main = txt
				)
			grid()
			j = j + 1
		
			#-------------------------------------------
		
		}
	}
}
		
#-------------------------------------------------------------------------------- 
# End ...
		
if (pdfout == 1) 	
	dev.off()
		

summary(a)	
warnings() 
